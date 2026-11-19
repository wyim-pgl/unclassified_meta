"""Align each sample's unclassified reads back to the contigs and count fragments per contig."""
from __future__ import annotations

import shlex
from collections import Counter
from pathlib import Path

from .fragments import count_primary, run_bowtie2
from .runner import StepError, require, sh


def build_index(contigs: Path, prefix: Path, threads: int) -> None:
    require("bowtie2-build")
    if not Path(str(prefix) + ".1.bt2").exists() and not Path(str(prefix) + ".1.bt2l").exists():
        q = shlex.quote
        sh(f"bowtie2-build --threads {threads} {q(str(contigs))} {q(str(prefix))} > {q(str(prefix) + '.build.log')} 2>&1")


def count_sample(prefix: Path, r1: Path, r2: Path | None, out_counts: Path, threads: int,
                 preset: str = "--very-fast") -> dict:
    """Write 'count<TAB>contig' lines, one fragment (pair or single read) counted once: the contig of mate 1's
    primary alignment, or of mate 2's when mate 1 did not align (see fragments.py, 'primary' policy)."""
    res = run_bowtie2(prefix, r1, r2, Path(str(out_counts) + ".bowtie2.log"), threads, preset, count_primary)
    counts, seen = res["result"]
    if res["input_reads"] is not None and res["input_reads"] != seen:
        raise StepError(f"{out_counts.name}: {seen} fragments in the SAM stream but Bowtie 2 read {res['input_reads']}")
    with open(out_counts, "w") as fo:
        for contig, n in sorted(counts.items()):
            fo.write(f"{n}\t{contig}\n")
    return res


def read_counts(path: Path) -> Counter:
    c = Counter()
    for line in open(path):
        if line.strip():
            n, contig = line.rstrip("\n").split("\t")
            c[contig] += int(n)
    return c


def category_counts(contig_counts: Counter, contig_category: dict[str, str], n_reads: int) -> dict[str, int]:
    """Sum read counts per category; reads that did not align are 'not_assembled'."""
    out: Counter = Counter()
    for contig, n in contig_counts.items():
        out[contig_category.get(contig, "unclassified")] += n
    aligned = sum(out.values())
    if aligned > n_reads:
        raise ValueError(f"{aligned} aligned reads exceed the {n_reads} input reads")
    out["not_assembled"] = n_reads - aligned
    return dict(out)
