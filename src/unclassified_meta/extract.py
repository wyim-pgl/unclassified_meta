"""Extract the reads a classifier left unclassified (Kraken2, Centrifuge, Centrifuger per-read outputs)."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from .io_utils import normalise_read_id, open_text


def parse_kraken2(path: str | Path) -> tuple[set[str], int, int]:
    """Kraken2 per-read output: 'C|U <tab> readID <tab> taxID ...'. Returns (unclassified IDs, total, classified)."""
    unclassified, total, classified = set(), 0, 0
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 3 or f[0] not in ("C", "U"):
                raise ValueError(f"{path}: not Kraken2 per-read output (line: {line[:80]!r}); "
                                 "pass the --output file, not the --report file")
            total += 1
            if f[0] == "U":
                unclassified.add(normalise_read_id(f[1]))
            else:
                classified += 1
    return unclassified, total, classified


def parse_centrifuge(path: str | Path) -> tuple[set[str], int, int]:
    """Centrifuge / Centrifuger per-read output: 'readID seqID taxID score ...' (header line optional).

    A read may have several lines (one per equally good assignment); it is classified when any line
    carries a taxID other than 0 and a seqID other than 'unclassified'.
    """
    status: dict[str, bool] = defaultdict(bool)
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "readID":
                continue
            if len(f) < 3:
                raise ValueError(f"{path}: not Centrifuge per-read output (line: {line[:80]!r})")
            rid = normalise_read_id(f[0])
            hit = f[1] != "unclassified" and f[2].strip() not in ("0", "")
            status[rid] = status[rid] or hit
    unclassified = {r for r, hit in status.items() if not hit}
    return unclassified, len(status), len(status) - len(unclassified)


def parse_classifier(path: str | Path, classifier: str) -> tuple[set[str], int, int]:
    if classifier == "kraken2":
        return parse_kraken2(path)
    if classifier in ("centrifuge", "centrifuger"):
        return parse_centrifuge(path)
    raise ValueError(f"unknown classifier {classifier!r}")


def filter_fastq(src: str | Path, dst: str | Path, keep: set[str]) -> int:
    """Copy FASTQ records whose normalised ID is in `keep`. Returns the number written."""
    n = 0
    with open_text(src) as fi, open_text(dst, "wt") as fo:
        while True:
            head = fi.readline()
            if not head:
                break
            seq, plus, qual = fi.readline(), fi.readline(), fi.readline()
            if not qual:
                raise ValueError(f"{src}: truncated FASTQ record at {head.strip()!r}")
            if normalise_read_id(head) in keep:
                fo.write(head + seq + plus + qual)
                n += 1
    return n
