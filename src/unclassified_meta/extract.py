"""Extract the reads a classifier left unclassified (Kraken2, Centrifuge, Centrifuger per-read outputs).

The optional score / hit-length thresholds for Centrifuge(r) output (idea contributed by Eli Holschbach): a
classification below the thresholds used for `centrifuger-quant` must not disappear from both the quantified and
the unclassified reads.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import NamedTuple

from .io_utils import normalise_read_id, open_text, refuse_overlap, staged


class Parsed(NamedTuple):
    """Result of reading one per-read classifier output."""
    unclassified: set       # read IDs to extract: no hit at all, or only hits below the thresholds
    total: int              # reads (or read pairs) in the classifier output
    classified: int         # reads with at least one assignment that passes the thresholds
    low_score: int          # reads assigned by the classifier, but only below --min-score / --min-hitlen
    all_ids: set = set()    # every read ID of the classifier output (the FASTQ must hold exactly these)


def parse_kraken2(path: str | Path) -> Parsed:
    """Kraken2 per-read output: 'C|U <tab> readID <tab> taxID ...'."""
    unclassified, total, classified = set(), 0, 0
    seen: set[str] = set()
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 3 or f[0] not in ("C", "U"):
                raise ValueError(f"{path}: not Kraken2 per-read output (line: {line[:80]!r}); "
                                 "pass the --output file, not the --report file")
            rid = normalise_read_id(f[1])
            if rid in seen:
                raise ValueError(f"{path}: read {rid!r} occurs twice; Kraken2 writes one line per read")
            seen.add(rid)
            total += 1
            if f[0] == "U":
                unclassified.add(rid)
            else:
                classified += 1
    return Parsed(unclassified, total, classified, 0, seen)


def centrifuge_reads(path: str | Path, min_score: int = 0, min_hitlen: int = 0):
    """The one place that decides what a Centrifuge / Centrifuger read is. Yields (read ID, status, taxIDs).

    Per-read output (header line optional): readID seqID taxID score 2ndBestScore hitLength queryLength numMatches.
    A read may have several lines (one per reported assignment, `-k` > 1); the tools write them one after the
    other, and consecutive lines with the same read ID are taken as one read. A read ID that reappears after
    other reads (a concatenated or re-sorted file) is an error, whatever its status: otherwise one read could be
    counted as classified from one block and extracted as unclassified from another.

    A line is an accepted assignment when seqID is not 'unclassified', taxID is not 0, score >= min_score and
    hitLength >= min_hitlen: the rule of `centrifuger-quant --min-score/--min-length`.
    status 2  classified: at least one accepted line; taxIDs = the taxIDs of the accepted lines, in file order
    status 1  low score:  assignments, but none accepted
    status 0  unclassified: no assignment
    `extract` (which reads to write) and `taxa-table --per-read` (which reads to count) both use this function,
    so a read is always on exactly one side.
    """
    need = 6 if (min_score or min_hitlen) else 3
    current, status, taxids = None, 0, []
    seen: set[str] = set()
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "readID":
                continue
            if len(f) < need:
                raise ValueError(f"{path}: not Centrifuge per-read output with {need} or more columns "
                                 f"(line: {line[:80]!r})")
            rid = normalise_read_id(f[0])
            if rid != current:
                if rid in seen:
                    raise ValueError(f"{path}: read {rid!r} occurs in separate blocks; the lines of a read must be "
                                     "consecutive and each read must occur once (was the file concatenated or sorted?)")
                seen.add(rid)
                if current is not None:
                    yield current, status, taxids
                current, status, taxids = rid, 0, []
            taxid = f[2].strip()
            if f[1] != "unclassified" and taxid not in ("0", ""):
                if need == 6 and (int(float(f[3])) < min_score or int(float(f[5])) < min_hitlen):
                    status = max(status, 1)
                else:
                    status = 2
                    if taxid not in taxids:
                        taxids.append(taxid)
    if current is not None:
        yield current, status, taxids


def parse_centrifuge(path: str | Path, min_score: int = 0, min_hitlen: int = 0) -> Parsed:
    """Reads to extract from Centrifuge / Centrifuger per-read output: those without an accepted assignment
    (see `centrifuge_reads`). Low-score and unclassified reads are both extracted, so no read is lost between
    quantification and this pipeline."""
    unclassified, total, n_unclassified, low = set(), 0, 0, 0
    all_ids: set[str] = set()
    for rid, status, _ in centrifuge_reads(path, min_score, min_hitlen):     # refuses reopened read IDs
        total += 1
        all_ids.add(rid)
        if status < 2:
            unclassified.add(rid)
            n_unclassified += 1
            low += status == 1
    return Parsed(unclassified, total, total - n_unclassified, low, all_ids)


def parse_classifier(path: str | Path, classifier: str, min_score: int = 0, min_hitlen: int = 0) -> Parsed:
    if classifier == "kraken2":
        if min_score or min_hitlen:
            raise ValueError("--min-score / --min-hitlen apply to Centrifuge and Centrifuger per-read output only; "
                             "Kraken2 output has no score or hit-length column (use Kraken2 --confidence when classifying)")
        return parse_kraken2(path)
    if classifier in ("centrifuge", "centrifuger"):
        return parse_centrifuge(path, min_score, min_hitlen)
    raise ValueError(f"unknown classifier {classifier!r}")


def check_output_is_not_an_input(dst: str | Path, inputs) -> None:
    """Refuse an output that is, or links to, one of the inputs (see io_utils.refuse_overlap). Commands also run
    the same check over all outputs and all inputs of all samples before they write anything."""
    refuse_overlap([dst], inputs, "extract")


def _records(path: str | Path):
    """Yield (normalised read ID, the four FASTQ lines) of every record."""
    with open_text(path) as fi:
        while True:
            head = fi.readline()
            if not head:
                return
            seq, plus, qual = fi.readline(), fi.readline(), fi.readline()
            if not qual or not head.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"{path}: truncated or malformed FASTQ record at {head.strip()!r}")
            yield normalise_read_id(head), head + seq + plus + qual


def _claim(rid: str, allowed: set | None, where: str) -> None:
    """Each FASTQ read must be a read of the classifier output, and only once (`allowed` is consumed)."""
    if allowed is None:
        return
    if rid not in allowed:
        raise ValueError(f"{where}: read {rid!r} is not in the classifier output, or occurs twice in the FASTQ; "
                         "the FASTQ must hold exactly the reads the classifier saw")
    allowed.discard(rid)


def _check_left(allowed: set | None, where: str) -> None:
    if allowed:
        example = sorted(allowed)[:3]
        raise ValueError(f"{where}: {len(allowed)} reads of the classifier output are not in the FASTQ "
                         f"(e.g. {', '.join(example)}); are the reads and the classifier output from the same run?")


def filter_fastq(src: str | Path, dst: str | Path, keep: set[str] | None, allowed: set | None = None) -> int:
    """Copy FASTQ records whose normalised ID is in `keep` (all when keep is None) to `dst`, staged and renamed.
    With `allowed` (the classifier's read IDs, consumed) every record must be one of them, exactly once, and
    none may be missing."""
    dst = Path(dst)
    check_output_is_not_an_input(dst, [src])
    n = 0
    with staged(dst) as tmp, open_text(tmp, "wt") as fo:
        for rid, rec in _records(src):
            _claim(rid, allowed, str(src))
            if keep is None or rid in keep:
                fo.write(rec)
                n += 1
        _check_left(allowed, str(src))
    return n


def filter_fastq_pair(r1: str | Path, r2: str | Path, o1: str | Path, o2: str | Path, keep: set[str] | None,
                      allowed: set | None = None) -> int:
    """Paired version of `filter_fastq`: both files are read in step, and the normalised read IDs (without /1, /2
    and comments) must agree record by record, so the mates written stay pairs. Used with keep=None (all reads)
    for samples without classifier output, so their mates are checked too."""
    o1, o2 = Path(o1), Path(o2)
    for o in (o1, o2):
        check_output_is_not_an_input(o, [r1, r2])
    n_rec = n = 0
    with staged(o1) as t1, staged(o2) as t2, open_text(t1, "wt") as f1, open_text(t2, "wt") as f2:
        it1, it2 = _records(r1), _records(r2)
        while True:
            a, b = next(it1, None), next(it2, None)
            if a is None and b is None:
                break
            if a is None or b is None:
                raise ValueError(f"{r1} and {r2} hold different numbers of records: the mate files are not in sync")
            n_rec += 1
            if a[0] != b[0]:
                raise ValueError(f"mate files out of sync at record {n_rec}: {a[0]!r} in {Path(r1).name} but {b[0]!r} "
                                 f"in {Path(r2).name}; the two files must list the same reads in the same order")
            _claim(a[0], allowed, str(r1))
            if keep is None or a[0] in keep:
                f1.write(a[1])
                f2.write(b[1])
                n += 1
        _check_left(allowed, str(r1))
    return n


def copy_reads(src: str | Path, dst: str | Path) -> None:
    """All reads of a single-end input (classifier 'none'): a real gzip copy, never a link, staged and renamed."""
    import gzip
    import shutil
    dst = Path(dst)
    check_output_is_not_an_input(dst, [src])
    with staged(dst) as tmp:
        if str(src).endswith(".gz"):
            shutil.copyfile(src, tmp)
        else:
            with open(src, "rb") as fi, gzip.open(tmp, "wb", compresslevel=1) as fo:
                shutil.copyfileobj(fi, fo)
