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
