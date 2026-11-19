"""Run external commands with logging, and mark finished steps so a rerun resumes."""
from __future__ import annotations

import hashlib
import json
import logging
import os
import shlex
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger("unclassified_meta")


class StepError(RuntimeError):
    pass


def open_log(logfile: Path):
    """Open a log file for appending without following a symbolic link (O_NOFOLLOW) and refuse a file with a
    second hard link: the log can never write into a file outside the output folder."""
    import os
    logfile.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_APPEND | getattr(os, "O_NOFOLLOW", 0)
    if not hasattr(os, "O_NOFOLLOW") and logfile.is_symlink():
        raise OSError(f"{logfile} is a symbolic link; refusing to write the log through it")
    try:
        fd = os.open(logfile, flags, 0o644)
    except OSError as e:
        raise OSError(f"cannot open the log {logfile} safely ({e.strerror}; a symbolic link is refused)") from None
    if os.fstat(fd).st_nlink > 1:
        os.close(fd)
        raise OSError(f"{logfile} has another hard link; refusing to write the log into it")
    return os.fdopen(fd, "a")


class _LogHandler(logging.StreamHandler):
    """File handler on a stream opened by `open_log`."""

    def __init__(self, logfile: Path):
        super().__init__(open_log(logfile))
        self.baseFilename = str(logfile.resolve())

    def close(self):
        try:
            self.flush()
            self.stream.close()
        finally:
            super().close()


def setup_logging(logfile: Path) -> None:
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S")
    root = logging.getLogger("unclassified_meta")
    root.setLevel(logging.INFO)
    if not any(getattr(h, "baseFilename", None) == str(logfile.resolve()) for h in root.handlers):
        fh = _LogHandler(logfile)
        fh.setFormatter(fmt)
        root.addHandler(fh)
    if not any(type(h) is logging.StreamHandler for h in root.handlers):
        sh = logging.StreamHandler()
        sh.setFormatter(fmt)
        root.addHandler(sh)


def require(*tools: str) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise StepError(f"required tool(s) not on PATH: {', '.join(missing)} "
                        "(activate the micromamba environment from environment.yml)")


def sh(cmd: str | list[str], stdout: Path | None = None, cwd: Path | None = None) -> None:
    """Run a command (list = exec, str = bash -o pipefail). Raise StepError with the stderr tail on failure."""
    text = cmd if isinstance(cmd, str) else " ".join(shlex.quote(c) for c in cmd)
    log.info("RUN %s%s", text, f" > {stdout}" if stdout else "")
    args = ["bash", "-o", "pipefail", "-c", cmd] if isinstance(cmd, str) else cmd
    out = open(stdout, "w") if stdout else subprocess.PIPE
    try:
        proc = subprocess.run(args, stdout=out, stderr=subprocess.PIPE, text=True, cwd=cwd)
    finally:
        if stdout:
            out.close()
    if proc.returncode != 0:
        tail = "\n".join((proc.stderr or "").strip().splitlines()[-20:])
        raise StepError(f"command failed (exit {proc.returncode}): {text}\n{tail}")


# Part of every step key. Bump it whenever a counting rule, a file layout or the meaning of a cached result
# changes, so results of older versions are recomputed instead of reused.
CACHE_VERSION = "counting-policy-2026-11-18"
SAMPLE_BYTES = 1 << 20          # file_id: bytes read from the start and from the end of a file


def params_key(*parts) -> str:
    """Short hash of the cache version and the parameters and inputs a step depends on (order-sensitive)."""
    return hashlib.sha256(json.dumps([CACHE_VERSION, *parts], sort_keys=True, default=str).encode()).hexdigest()[:16]


def _digest(p: Path, size: int) -> str:
    """sha256 of the first and the last SAMPLE_BYTES of a file (of all of it when it is small, or when the
    environment variable UNCLASSIFIED_META_FULL_HASH=1 asks for full-file hashing)."""
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        if size <= 2 * SAMPLE_BYTES or os.environ.get("UNCLASSIFIED_META_FULL_HASH") == "1":
            for chunk in iter(lambda: fh.read(1 << 22), b""):
                h.update(chunk)
        else:
            h.update(fh.read(SAMPLE_BYTES))
            fh.seek(size - SAMPLE_BYTES)
            h.update(fh.read(SAMPLE_BYTES))
    return h.hexdigest()[:16]


def file_id(path) -> list:
    """Identity of an input file or folder: absolute path, size, modification time in nanoseconds and a digest of
    its first and last megabyte (whole file below 2 MB; set UNCLASSIFIED_META_FULL_HASH=1 to hash whole files).
    Two different files therefore never share an identity because of a matching name, size and second-rounded time,
    and a file replaced in place is noticed. Hashing whole FASTQ files would cost minutes per sample, hence the
    sampled digest. A folder (e.g. a Kraken2 database) is the list of the identities of the files directly in it.
    Moving or copying inputs changes the identity, so finished steps are recomputed."""
    if path is None:
        return [None]
    p = Path(path)
    if not p.exists():
        return [str(p.absolute()), None]
    p = p.resolve()
    if p.is_dir():
        return [str(p), "dir", [file_id(c) for c in sorted(p.iterdir()) if c.is_file()]]
    st = p.stat()
    return [str(p), st.st_size, st.st_mtime_ns, _digest(p, st.st_size)]


class Steps:
    """Done-markers under <workdir>/.done so finished steps are skipped on a rerun.

    A marker stores the key of the parameters and inputs the step ran with (keys include CACHE_VERSION).
    `done(name, key)` is true only when the stored key equals `key`, so a step is recomputed when, for example,
    --subsample, an input file or the counting policy changes. Markers written by versions before keys existed
    are empty; the parameters they ran with are unknown, so they are not trusted and the step is recomputed."""

    def __init__(self, workdir: Path):
        self.dir = workdir / ".done"
        self.dir.mkdir(parents=True, exist_ok=True)

    def done(self, name: str, key: str | None = None) -> bool:
        """True only when the marker's key equals `key` AND every output recorded with the marker still has the
        recorded checksum (an output changed or removed since then makes the step run again)."""
        f = self.dir / name
        if not f.exists():
            return False
        text = f.read_text().strip()
        if not text:
            log.warning("step %s: finished by an older version without recorded parameters; recomputing it", name)
            return False
        try:
            rec = json.loads(text)
        except ValueError:
            rec = {"key": text, "outputs": {}}
        if key is not None and rec.get("key") != key:
            log.info("step %s: parameters or inputs changed since it last ran; recomputing", name)
            return False
        if not rec.get("outputs"):
            # every step records its outputs; a marker without them cannot show that they are still intact
            log.info("step %s: its marker records no outputs; recomputing", name)
            return False
        for out, digest in (rec.get("outputs") or {}).items():
            if not Path(out).exists() or _sha256(out) != digest:
                log.info("step %s: its output %s changed since it ran; recomputing", name, out)
                return False
        return True

    def invalidate(self, name: str) -> None:
        """Remove a step's marker BEFORE its outputs are changed: if the step then fails, no marker of an earlier
        run can make a later run reuse half-replaced outputs."""
        (self.dir / name).unlink(missing_ok=True)

    def mark(self, name: str, key: str | None = None, outputs=()) -> None:
        """Publish the marker last, after every output of the step is in place: the key and the sha256 of each
        output (verified by `done` before the step is reused)."""
        from .io_utils import staged
        rec = {"key": key or "", "outputs": {str(Path(o).resolve()): _sha256(o) for o in outputs if o is not None}}
        if not rec["outputs"]:
            raise ValueError(f"step {name}: a step marker must record the step's output files")
        with staged(self.dir / name) as tmp:
            tmp.write_text(json.dumps(rec))

    def output_digest(self, name: str, output) -> str | None:
        """The checksum recorded for one output of a finished step (feeds the keys of the steps that use it)."""
        f = self.dir / name
        try:
            return json.loads(f.read_text()).get("outputs", {}).get(str(Path(output).resolve()))
        except (OSError, ValueError):
            return None


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()
