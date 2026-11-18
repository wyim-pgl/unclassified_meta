"""Align each sample's unclassified reads back to the contigs and count reads per contig."""
from __future__ import annotations

import shlex
from collections import Counter
from pathlib import Path

from .runner import require, sh


def build_index(contigs: Path, prefix: Path, threads: int) -> None:
    require("bowtie2-build")
    if not Path(str(prefix) + ".1.bt2").exists() and not Path(str(prefix) + ".1.bt2l").exists():
        q = shlex.quote
        sh(f"bowtie2-build --threads {threads} {q(str(contigs))} {q(str(prefix))} > {q(str(prefix) + '.build.log')} 2>&1")


def count_sample(prefix: Path, r1: Path, r2: Path | None, out_counts: Path, threads: int,
                 preset: str = "--very-fast") -> None:
    """Write 'count<TAB>contig' lines: primary alignments of the first mate (paired) or of each read (single)."""
    require("bowtie2", "samtools")
    q = lambda x: shlex.quote(str(x))
    reads = f"-1 {q(r1)} -2 {q(r2)}" if r2 else f"-U {q(r1)}"
    flt = "-F 2308 -f 64" if r2 else "-F 2308"   # 2308 = unmapped + secondary + supplementary
    sh(f"bowtie2 -p {threads} {preset} --no-unal -x {q(prefix)} {reads} 2> {q(str(out_counts) + '.bowtie2.log')} "
       f"| samtools view {flt} - | cut -f3 | sort | uniq -c | awk '{{print $1\"\\t\"$2}}' > {q(out_counts)}")


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
