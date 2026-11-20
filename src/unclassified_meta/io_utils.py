"""Small FASTA/FASTQ helpers (pure Python, no third-party dependencies)."""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator


def open_text(path: str | Path, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        if "w" in mode:
            return gzip.open(path, mode, compresslevel=1)
        return gzip.open(path, mode)
    return open(path, mode.replace("t", "") if "b" in mode else mode)


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield (id, sequence); id is the first word of the header."""
    name, chunks = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0] if line[1:].strip() else "", []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def write_fasta(records, path: str | Path, width: int = 80) -> int:
    n = 0
    with open_text(path, "wt") as fo:
        for name, seq in records:
            fo.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fo.write(seq[i:i + width] + "\n")
            n += 1
    return n


def count_fastq(path: str | Path) -> int:
    """Number of records in a FASTQ (4 lines each)."""
    n = 0
    with open_text(path) as fh:
        for n, _ in enumerate(fh, 1):
            pass
    if n % 4:
        raise ValueError(f"{path}: line count {n} is not a multiple of 4")
    return n // 4


def normalise_read_id(read_id: str) -> str:
    """Strip a leading '@', anything after whitespace and a trailing /1 or /2 mate suffix."""
    rid = read_id.lstrip("@").split()[0] if read_id.strip() else ""
    if len(rid) > 2 and rid[-2] == "/" and rid[-1] in "12":
        rid = rid[:-2]
    return rid


# ---------------------------------------------------------------- output safety (used by every command)

def same_file(a: str | Path, b: str | Path) -> bool:
    """True when two paths name the same file: same resolved path (symbolic links followed), or the same
    device and inode (hard links)."""
    import os
    pa, pb = Path(a), Path(b)
    if os.path.realpath(pa) == os.path.realpath(pb):
        return True
    try:
        sa, sb = pa.stat(), pb.stat()
    except OSError:
        return False
    return (sa.st_dev, sa.st_ino) == (sb.st_dev, sb.st_ino)


def tmp_for(dst: str | Path) -> Path:
    """Pattern of the temporary names used for `dst` (unique names are made by `staged`; this is the stem only,
    used by overlap checks)."""
    dst = Path(dst)
    name = dst.name
    return dst.with_name(("." + name[:-3] + ".part.gz") if name.endswith(".gz") else ("." + name + ".part"))


def refuse_overlap(outputs, inputs, command: str = "this command") -> None:
    """Refuse when any output is any input, through the same path, a symbolic link or a hard link (same device
    and inode). Call it with ALL outputs and ALL inputs of a command, across samples and mates."""
    ins = [Path(i) for i in inputs if i is not None]
    for o in outputs:
        if o is None:
            continue
        for i in ins:
            if same_file(o, i):
                raise ValueError(f"{command}: output {o} is the input file {i} (same path or a link to it); "
                                 "write the results to a folder that holds no input files")


def preflight(command: str, inputs, owned_dirs=(), files=()) -> None:
    """The check every command runs before its FIRST write (log files and sample sheets included).

    owned_dirs  folders the command writes into, cleans or deletes (output, scratch and staging folders)
    files       single output files outside those folders
    Refused: an input that lies inside an owned folder (it could be overwritten or deleted with it), an owned
    folder that is or lies inside an input, any symbolic link (to a file or a folder) and any file with more
    than one hard link inside an owned folder, and any output file that is an input. The owned folder itself
    may be a symbolic link (e.g. to a scratch file system); its real path is used."""
    import os
    ins = [Path(i) for i in inputs if i is not None]
    real_in = [Path(os.path.realpath(i)) for i in ins]
    ids = set()
    for i in ins:
        try:
            st = i.stat()
            ids.add((st.st_dev, st.st_ino))
        except OSError:
            pass
    for d in owned_dirs:
        if d is None:
            continue
        rd = Path(os.path.realpath(d))
        for i, ri in zip(ins, real_in):
            if ri == rd or rd in ri.parents:
                raise ValueError(f"{command}: input {i} lies inside the output folder {d}; inputs must stay outside "
                                 "the folders a command writes, cleans or deletes")
            if ri in rd.parents:
                raise ValueError(f"{command}: output folder {d} lies inside the input {i}")
        if rd.is_dir():
            # nothing inside an owned folder may be a link: a symbolic link (to a file or a folder) or a file with
            # a second hard link would let a write or a clean-up reach a file outside the folder
            for root, dirs, fnames in os.walk(rd, followlinks=False):
                for fn in fnames + dirs:
                    f = Path(root) / fn
                    try:
                        st = os.lstat(f)
                    except OSError:
                        continue
                    if (st.st_dev, st.st_ino) in ids or (os.path.islink(f) and f.resolve() in real_in):
                        raise ValueError(f"{command}: {f} in the output folder {d} is (a link to) an input file")
                    if os.path.islink(f):
                        raise ValueError(f"{command}: {f} in the output folder {d} is a symbolic link; output "
                                         "folders must hold no links (remove it, or choose another output folder)")
                    if not os.path.isdir(f) and st.st_nlink > 1:
                        raise ValueError(f"{command}: {f} in the output folder {d} has another hard link "
                                         f"({st.st_nlink} names); output folders must hold no links")
    refuse_overlap(files, ins, command)


class folder_lock:
    """Exclusive lock of an output folder for one command (fcntl lock on a lock file in it). A second command on
    the same folder fails at once instead of mixing its files with the first one's."""

    def __init__(self, folder: str | Path, command: str):
        self.path = Path(folder) / ".unclassified_meta.lock"
        self.command = command
        self.fh = None

    def __enter__(self):
        import fcntl
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = open(self.path, "a+")
        try:
            fcntl.flock(self.fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            self.fh.close()
            raise RuntimeError(f"another run is using this folder ({self.path.parent}); wait for it to finish or "
                               "choose another output folder") from None
        self.fh.seek(0)
        self.fh.truncate()
        self.fh.write(f"{self.command} pid {__import__('os').getpid()}\n")
        self.fh.flush()
        return self

    def __exit__(self, *exc):
        import fcntl
        fcntl.flock(self.fh, fcntl.LOCK_UN)
        self.fh.close()
        return False


def current_umask() -> int:
    """The process umask (read by setting and restoring it; there is no read-only call)."""
    import os
    mask = os.umask(0o022)
    os.umask(mask)
    return mask


class staged:
    """Context manager: write to a UNIQUE temporary file in the destination folder (created exclusively, so two
    writers never share one), then os.replace it into place. On an error, or when the rename itself fails, only
    that temporary file is removed; nothing else is ever deleted. Replacing a symbolic link at the destination
    replaces the link itself, never the file it points to."""

    def __init__(self, dst: str | Path):
        self.dst = Path(dst)
        self.tmp = None

    def __enter__(self) -> Path:
        import os
        import tempfile
        self.dst.parent.mkdir(parents=True, exist_ok=True)
        name = self.dst.name
        suffix = ".part.gz" if name.endswith(".gz") else ".part"
        fd, path = tempfile.mkstemp(prefix=f".{name}.", suffix=suffix, dir=self.dst.parent)
        try:
            # mkstemp creates 0600; give the output the mode any other new file would get under the umask
            os.fchmod(fd, 0o666 & ~current_umask())
        finally:
            os.close(fd)
        self.tmp = Path(path)
        return self.tmp

    def __exit__(self, exc_type, exc, tb) -> bool:
        import os
        if exc_type is None:
            try:
                os.replace(self.tmp, self.dst)
            except BaseException:
                self.tmp.unlink(missing_ok=True)
                raise
        else:
            self.tmp.unlink(missing_ok=True)
        return False


def make_staging_dir(parent: str | Path, owner: str) -> Path:
    """A fresh, uniquely named staging folder with an ownership marker; only such folders may be removed by
    `remove_staging_dir`."""
    import tempfile
    d = Path(tempfile.mkdtemp(prefix=".staging-", dir=parent))
    (d / ".owned_by_unclassified_meta").write_text(owner + "\n")
    return d


def remove_staging_dir(d: str | Path) -> None:
    """Remove a staging folder made by `make_staging_dir` (never a symbolic link, never a folder without the
    ownership marker)."""
    import shutil
    d = Path(d)
    if d.is_symlink() or not (d / ".owned_by_unclassified_meta").is_file():
        raise RuntimeError(f"refusing to delete {d}: not a staging folder made by this program")
    shutil.rmtree(d)
