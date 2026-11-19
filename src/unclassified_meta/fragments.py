"""One fragment-counting policy for Bowtie 2 output, shared by the read-back step and `panel-quant`.

A *fragment* is one read pair (paired-end) or one read (single-end). Bowtie 2 writes all records of a fragment
next to each other, so the SAM stream is grouped by read name without sorting.

Two modes:

primary  (read-back in `run`)   the fragment goes to the sequence of mate 1's primary alignment, or to mate 2's
                                when mate 1 did not align. R1-only and R2-only fragments are both counted.
feature  (`panel-quant`)        with `bowtie2 -k N`, the best-scoring alignments of each mate (all alignments
                                whose score AS equals that mate's best) are looked up in a sequence -> feature
                                map. The fragment counts once for a feature only when all of them fall in
                                that single feature; everything else goes to a ledger outcome (ambiguous,
                                ambiguous_decoy, decoy_*, unmapped). Lower-scoring alignments are ignored: a
                                read that matches its own contig exactly and a divergent reference poorly
                                belongs to the contig.
"""
from __future__ import annotations

import logging
import re
import resource
import shlex
import subprocess
import time
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator

from .runner import StepError, require

log = logging.getLogger("unclassified_meta")

DECOY_PREFIX = "decoy_"          # feature-map values starting with this are decoy classes, not features
UNMAPPED, AMBIGUOUS, AMBIGUOUS_DECOY, DECOY_MIXED = "unmapped", "ambiguous", "ambiguous_decoy", "decoy_mixed"


def sam_fragments(lines: Iterable[str], scores: bool = False) -> Iterator[list[tuple]]:
    """Group consecutive SAM records by read name -> [(flag, reference name[, alignment score]), ...] per
    fragment. The alignment score (AS:i tag; None when absent) is parsed only when `scores` is set."""
    name, recs = None, []
    for line in lines:
        if not line or line[0] == "@":
            continue
        f = line.split("\t") if scores else line.split("\t", 3)
        if len(f) < 3:
            continue
        if f[0] != name:
            if name is not None:
                yield recs
            name, recs = f[0], []
        if scores:
            a = None
            for tag in f[11:]:
                if tag.startswith("AS:i:"):
                    a = int(tag[5:])
                    break
            recs.append((int(f[1]), f[2], a))
        else:
            recs.append((int(f[1]), f[2]))
    if name is not None:
        yield recs


def primary_target(recs: list[tuple[int, str]]) -> str | None:
    """Sequence of the fragment under the 'primary' policy (mate 1 first, then mate 2), or None when unaligned."""
    first = second = None
    for flag, ref, *_ in recs:
        if flag & 0x904 or ref == "*":      # unmapped, secondary or supplementary
            continue
        if flag & 0x80:                     # mate 2
            second = second or ref
        else:                               # mate 1 or a single-end read
            first = first or ref
    return first or second


def best_refs(recs: list[tuple]) -> set[str]:
    """Reference names of the best-scoring alignments of each mate (ties at the mate's best AS are all kept)."""
    best: dict[int, int] = {}
    aligned = []
    for flag, ref, *rest in recs:
        if flag & 0x4 or ref == "*":
            continue
        mate = 2 if flag & 0x80 else 1
        score = rest[0] if rest and rest[0] is not None else 0
        aligned.append((mate, ref, score))
        if mate not in best or score > best[mate]:
            best[mate] = score
    return {ref for mate, ref, score in aligned if score == best[mate]}


def feature_outcome(recs: list[tuple], seq2feature: dict[str, str]) -> str:
    """Outcome of one fragment under the 'feature' policy: a feature id or a ledger class."""
    targets, decoys = set(), set()
    for ref in best_refs(recs):
        feat = seq2feature.get(ref)
        if feat is None:
            raise KeyError(f"reference sequence {ref!r} is not in the sequence -> feature map")
        (decoys if feat.startswith(DECOY_PREFIX) else targets).add(feat)
    if not targets and not decoys:
        return UNMAPPED
    if targets and decoys:
        return AMBIGUOUS_DECOY
    if targets:
        return next(iter(targets)) if len(targets) == 1 else AMBIGUOUS
    return next(iter(decoys)) if len(decoys) == 1 else DECOY_MIXED


def count_primary(lines: Iterable[str]) -> tuple[Counter, int]:
    """(fragments per sequence, number of fragments seen) under the 'primary' policy."""
    counts: Counter = Counter()
    n = 0
    for recs in sam_fragments(lines):
        n += 1
        t = primary_target(recs)
        if t is not None:
            counts[t] += 1
    return counts, n


def count_features(lines: Iterable[str], seq2feature: dict[str, str]) -> tuple[Counter, Counter, int]:
    """(fragments per feature, ledger of the other outcomes, number of fragments seen) under the 'feature' policy."""
    feats: Counter = Counter()
    ledger: Counter = Counter()
    n = 0
    for recs in sam_fragments(lines, scores=True):
        n += 1
        o = feature_outcome(recs, seq2feature)
        if o in (UNMAPPED, AMBIGUOUS, AMBIGUOUS_DECOY, DECOY_MIXED) or o.startswith(DECOY_PREFIX):
            ledger[o] += 1
        else:
            feats[o] += 1
    return feats, ledger, n


def bowtie2_input_reads(logfile: Path) -> int | None:
    """Number of reads/pairs Bowtie 2 says it processed ('N reads; of these:')."""
    m = re.search(r"^(\d+) reads; of these:", Path(logfile).read_text(errors="replace"), re.M)
    return int(m.group(1)) if m else None


def run_bowtie2(prefix: Path, r1: Path, r2: Path | None, logfile: Path, threads: int, options: str,
                consume) -> dict:
    """Run Bowtie 2, hand its SAM stream (text lines) to `consume`, and return {'result', 'wall_s', 'max_rss_gb',
    'input_reads'}. Fails when Bowtie 2 fails or when the fragments seen differ from Bowtie 2's own count."""
    require("bowtie2")
    reads = ["-1", str(r1), "-2", str(r2)] if r2 else ["-U", str(r1)]
    cmd = ["bowtie2", "-p", str(threads), *shlex.split(options), "--no-hd", "-x", str(prefix), *reads]
    log.info("RUN %s 2> %s | <count fragments>", " ".join(shlex.quote(c) for c in cmd), logfile)
    t0 = time.time()
    with open(logfile, "w") as err:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True, bufsize=1 << 20)
        try:
            result = consume(proc.stdout)
        finally:
            proc.stdout.close()
            rc = proc.wait()
    if rc != 0:
        tail = "\n".join(Path(logfile).read_text(errors="replace").strip().splitlines()[-15:])
        raise StepError(f"bowtie2 failed (exit {rc}); see {logfile}\n{tail}")
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024   # GB; maximum over children so far
    return {"result": result, "wall_s": round(time.time() - t0, 1), "max_rss_gb": round(rss, 3),
            "input_reads": bowtie2_input_reads(logfile)}
