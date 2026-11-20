"""strain-evidence: overlap-aware viral sequence clusters, local read phasing and evidence-graded strain labels.

Why: a 97/98/99 % identity cluster of contigs is an *operational sequence cluster*, not a strain. Two contigs of
one virus that cover different parts of the genome are not "different": they are a coverage difference. This
step therefore compares contigs only inside the reference interval both cover, and grades every cluster by the
evidence it actually has (docs/DESIGN_robustness.md, section 1):

  operational_sequence_cluster      single linkage at the identity threshold (always)
  candidate_strain_lineage          differs from another cluster inside an overlap >= --min-overlap, and no
                                    member is flagged as a possible chimera / artifact
  phased_viral_haplotype            a candidate whose distinguishing variants (>= 2 of them) are carried together
                                    by one supported local read haplotype (reads link them)
  reference_supported_viral_strain  no flagged member, and the cluster consensus matches its reference at
                                    >= --reference-identity over >= --reference-cover of the reference length
                                    (this level does not need reads; the separate evidence columns say what
                                    else holds)

Method
1. Contigs of the chosen categories with a nucleotide reference hit (BLASTN against the virus genomes or the
   extra virus FASTA: method blastn / blastn_extra) are grouped by that reference. Protein-only, Kraken2 and
   unmatched contigs cannot be placed on a genome and are reported as not analysable (contig_status.tsv).
2. All analysable contigs are aligned to the selected references with `minimap2 -c --cs -x asm20
   --secondary=no` (supplementary alignments kept). The cs strings are parsed here into reference intervals and
   per-position differences (SNVs, insertions, deletions) on reference coordinates.
3. Per contig, flags for possible chimeras / artifacts: split_alignment (several collinear pieces),
   non_collinear (pieces in another order or strand), multi_reference (>= 100 query bp aligned best to another
   selected reference), identity_jump (windowed identity on the reference, 200 bp windows, spread >= 5 points),
   low_aligned_fraction (< 80 % of the contig aligned to its reference), not_aligned. Flagged contigs are never
   used to link a cluster (a chimera would bridge two lineages); they stay singletons and keep their flags.
4. Every pair of contigs on one reference: overlap of their reference intervals, SNVs/indels between them inside
   that overlap only, identity over the overlap = 1 - (SNVs + indel bp) / overlap bp. relation:
   overlapping_identical, overlapping_divergent, or disjoint (overlap < --min-overlap: a coverage difference).
5. Clusters per --identity threshold: single linkage on unflagged pairs with overlap >= --min-overlap and overlap
   identity >= threshold. A singleton with no overlap >= --min-overlap to any contig is `unresolved_disjoint`
   (it may belong to any cluster). Two members of one cluster that overlap below the threshold are
   `internal_divergence` (single linkage chained them).
6. Reads (optional): every sample is aligned with Bowtie 2 to the selected references (pairs together,
   `--no-unal`); the SAM is saved (gzip) and read twice in Python (CIGAR walked against the reference sequence,
   which is equivalent to the MD tag here). Pass 1: pileup (bases with quality >= --min-baseq; where the mates of
   a fragment overlap only mate 1 counts, so one molecule is one observation). A site is polymorphic when a
   non-major allele (SNV or indel) has depth >= --min-depth, >= --min-alt reads, AF >= --min-af, reads on both
   strands, and a binomial upper-tail p-value against --error-rate below 0.05 / m (Bonferroni over the m tested
   alleles of that sample and reference). Sites where the major allele differs from the reference are reported
   as fixed. Pass 2: the allele of every fragment at every polymorphic site (mates merged; a mate conflict drops
   the site for that fragment). Pairs of sites: 2x2 tables of major / non-major alleles with r2 and D'. Local
   blocks: from left to right, consecutive sites are added while >= --min-link fragments carry an allele at every
   site of the block. Haplotypes of a block: complete allele patterns of those fragments, counts, frequencies and
   a percentile bootstrap interval (multinomial over the fragments, seeded). A pattern is `possible_error` when a
   more frequent pattern at Hamming distance 1 could produce its count by sequencing error (binomial,
   p >= 0.001 at --error-rate). Supported haplotype = count >= --min-hap-reads and not possible_error.

Limits (written into the tables as well)
- Blocks are LOCAL: their length is bounded by the fragment length (phasing_limits.tsv gives read length,
  insert-size median / 95th percentile, depth and the longest block). Separate blocks are never joined, so no
  genome-wide haplotype is produced.
- Reads are aligned to the selected references only; reads of strains far from them (end-to-end Bowtie 2) do not
  align, and reads of an unselected relative may align here.
- multi_reference sees only the references selected by the run's contigs, not the whole reference database.
- --max-fragments keeps the FIRST N fragments of each sample (Bowtie 2 -u), not a random subsample.
- The pileup is pure Python/numpy: fine for a few million fragments per sample on small viral genomes; use
  --max-fragments for deeper samples.
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import os
import re
import shutil
import time
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import __version__
from . import manifest
from .io_utils import make_staging_dir, open_text, refuse_overlap, remove_staging_dir, write_fasta
from .runner import file_id, require, sh

log = logging.getLogger("unclassified_meta")

NT_METHODS = ("blastn", "blastn_extra")      # contig methods that place a contig on a reference genome
WINDOW = 200                                 # bp, windowed identity of a contig along its reference
IDENTITY_JUMP = 5.0                          # percentage points between the best and worst window
MIN_ALIGNED_FRACTION = 0.8
MULTI_REF_MIN_BP = 100
MAX_INSERT = 1000                            # bowtie2 -X
SIG_ALPHA = 0.05                             # family-wise level of the variant test (Bonferroni)
ERROR_ALPHA = 0.001                          # possible_error: p >= this under the error model
BOOTSTRAPS = 1000
MAX_LINKAGE_FIGURES = 40
MAX_TRACK_ROWS = 200
CLASSES = ("operational_sequence_cluster", "candidate_strain_lineage", "phased_viral_haplotype",
           "reference_supported_viral_strain")
OUTPUT_TABLES = ("viral_cluster_overlap.tsv", "viral_strain_clusters.tsv", "contig_status.tsv", "variant_sites.tsv",
                 "variant_linkage.tsv", "viral_haplotype_evidence.tsv", "phasing_limits.tsv")


# ---------------------------------------------------------------- contig alignments (pure, unit-tested)

@dataclass(frozen=True)
class Diff:
    """One difference of a query from the reference. pos is 0-based on the reference: the first deleted base
    for a deletion, the base after the inserted bases for an insertion (the same convention as a CIGAR walk)."""
    pos: int
    kind: str           # snv | ins | del
    ref: str
    alt: str

    @property
    def site(self) -> tuple[int, str]:
        return (self.pos, "snv" if self.kind == "snv" else "indel")

    @property
    def allele(self) -> str:
        """Allele code shared with the read alleles: a base for an SNV, '+SEQ' for an insertion, '-L' for a
        deletion; '.' is the reference (no indel) allele of an indel site."""
        return self.alt if self.kind == "snv" else ("+" + self.alt if self.kind == "ins" else f"-{len(self.ref)}")


CS_TOKEN = re.compile(r":\d+|\*[a-z][a-z]|[+\-][a-z]+|=[A-Za-z]+|~[a-z]{2}\d+[a-z]{2}")


def parse_cs(cs: str, tstart: int) -> tuple[list[Diff], int]:
    """Differences of a minimap2 short cs string (`--cs`) starting at reference position `tstart` (0-based);
    returns (differences, reference end). cs is on the reference's forward strand for both query strands."""
    toks = CS_TOKEN.findall(cs)
    if "".join(toks) != cs:
        raise ValueError(f"unreadable cs string (not minimap2 --cs short form): {cs[:60]}")
    pos, out = tstart, []
    for t in toks:
        c = t[0]
        if c == ":":
            pos += int(t[1:])
        elif c == "=":
            pos += len(t) - 1
        elif c == "*":
            out.append(Diff(pos, "snv", t[1].upper(), t[2].upper()))
            pos += 1
        elif c == "+":
            out.append(Diff(pos, "ins", "", t[1:].upper()))
        elif c == "-":
            out.append(Diff(pos, "del", t[1:].upper(), ""))
            pos += len(t) - 1
        else:
            raise ValueError(f"spliced alignment in a cs string ({t}); contigs are aligned with -x asm20")
    return out, pos


@dataclass
class Segment:
    contig: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    ref: str
    rlen: int
    rstart: int
    rend: int
    mapq: int
    diffs: list = field(default_factory=list)


def parse_paf(lines) -> list[Segment]:
    """Alignments of a `minimap2 -c --cs` PAF (secondary alignments, tp:A:S, are skipped)."""
    out = []
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        tags = {t[:2]: t[5:] for t in f[12:]}
        if tags.get("tp") == "S":
            continue
        if "cs" not in tags:
            raise ValueError("PAF line without a cs tag: align with minimap2 -c --cs")
        rstart = int(f[7])
        diffs, rend = parse_cs(tags["cs"], rstart)
        if rend != int(f[8]):
            raise ValueError(f"cs string of {f[0]} on {f[5]} ends at {rend}, PAF says {f[8]}")
        out.append(Segment(f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), rstart, rend, int(f[11]),
                           diffs))
    return out


def merge_intervals(iv) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def intersect_intervals(a, b) -> list[tuple[int, int]]:
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def interval_length(iv) -> int:
    return sum(e - s for s, e in iv)


def window_identities(rstart: int, rend: int, diffs: list[Diff], window: int = WINDOW) -> list[float]:
    """Identity (%) of each full `window` of an aligned reference interval (difference columns: an SNV 1, an
    indel its length, counted in the window where it starts)."""
    n = (rend - rstart) // window
    cols = [0] * n
    for d in diffs:
        k = (d.pos - rstart) // window
        if 0 <= k < n:
            cols[k] += 1 if d.kind == "snv" else len(d.ref or d.alt)
    return [100.0 * (1 - c / window) for c in cols]


def _collinear(segs: list[Segment]) -> bool:
    """Pieces of one contig in the same order and orientation on the reference (a big indel, not a chimera)."""
    if len({s.strand for s in segs}) > 1:
        return False
    by_q = sorted(segs, key=lambda s: s.qstart)
    starts = [s.rstart for s in by_q]
    if by_q[0].strand == "-":
        starts = starts[::-1]
        by_q = by_q[::-1]
    return starts == sorted(starts) and all(a.rend <= b.rstart + 50 for a, b in zip(by_q, by_q[1:]))


@dataclass
class Track:
    """A contig placed on its reference: covered reference intervals and its alleles at the sites where it
    differs from the reference."""
    contig: str
    ref: str
    length: int
    intervals: list = field(default_factory=list)
    alleles: dict = field(default_factory=dict)        # site -> allele (differences from the reference only)
    ref_alleles: dict = field(default_factory=dict)    # site -> reference allele (base, or '.' for an indel site)
    deleted: list = field(default_factory=list)        # reference intervals deleted in the contig
    flags: list = field(default_factory=list)
    aligned_fraction: float = 0.0
    other_refs: list = field(default_factory=list)
    n_segments: int = 0

    @property
    def linkable(self) -> bool:
        return not self.flags


def build_track(contig: str, length: int, ref: str, segments: list[Segment]) -> Track:
    """Place one contig on its reference (its BLASTN best hit) from all of its alignments and set its flags."""
    on = [s for s in segments if s.ref == ref]
    other = [s for s in segments if s.ref != ref]
    t = Track(contig, ref, length, n_segments=len(on))
    if not on:
        t.flags.append("not_aligned")
        return t
    if len(on) > 1:
        t.flags.append("split_alignment" if _collinear(on) else "non_collinear")
    q_on = merge_intervals([(s.qstart, s.qend) for s in on])
    t.aligned_fraction = interval_length(q_on) / length if length else 0.0
    if t.aligned_fraction < MIN_ALIGNED_FRACTION:
        t.flags.append("low_aligned_fraction")
    if other:
        q_other = merge_intervals([(s.qstart, s.qend) for s in other])
        extra = interval_length(q_other) - interval_length(intersect_intervals(q_other, q_on))
        if extra >= MULTI_REF_MIN_BP:
            t.flags.append("multi_reference")
            t.other_refs = sorted({s.ref for s in other})
    wins = [w for s in on for w in window_identities(s.rstart, s.rend, s.diffs)]
    if len(wins) >= 2 and max(wins) - min(wins) >= IDENTITY_JUMP:
        t.flags.append("identity_jump")
    t.intervals = merge_intervals([(s.rstart, s.rend) for s in on])
    for s in sorted(on, key=lambda s: -(s.rend - s.rstart)):      # overlapping pieces: the longest one wins
        for d in s.diffs:
            if d.site not in t.alleles:
                t.alleles[d.site] = d.allele
                t.ref_alleles[d.site] = d.ref if d.kind == "snv" else "."
                if d.kind == "del":
                    t.deleted.append((d.pos, d.pos + len(d.ref)))
    return t


def _allele_len(a: str) -> int:
    if a.startswith("-"):
        return int(a[1:])
    if a.startswith("+"):
        return len(a) - 1
    return 0


def _inside(pos: int, kind: str, length: int, intervals) -> bool:
    """A site lies inside the intervals: an SNV its base; an indel its anchor base and every base it touches."""
    for s, e in intervals:
        if kind == "snv" and s <= pos < e:
            return True
        if kind == "indel" and s < pos and pos + max(length, 1) <= e:
            return True
    return False


def _in_any(pos: int, intervals) -> bool:
    return any(s <= pos < e for s, e in intervals)


def compare_tracks(a: Track, b: Track) -> dict:
    """Differences of two contigs on one reference INSIDE their overlap only. Returns overlap bp and intervals,
    n_snv, n_indel, indel_bp, identity (%, None without overlap) and diffs [(site, allele_a, allele_b)]."""
    ov = intersect_intervals(a.intervals, b.intervals)
    obp = interval_length(ov)
    diffs = []
    n_snv = n_indel = indel_bp = 0
    if obp:
        for site in sorted(set(a.alleles) | set(b.alleles)):
            pos, kind = site
            ref = a.ref_alleles.get(site) or b.ref_alleles.get(site)
            x, y = a.alleles.get(site, ref), b.alleles.get(site, ref)
            if x == y:
                continue
            if kind == "snv":
                if not _inside(pos, kind, 0, ov) or _in_any(pos, a.deleted) or _in_any(pos, b.deleted):
                    continue
                n_snv += 1
            else:
                ln = max(_allele_len(x), _allele_len(y))
                if not _inside(pos, kind, ln, ov):
                    continue
                n_indel += 1
                indel_bp += ln
            diffs.append((site, x, y))
    ident = 100.0 * (1 - (n_snv + indel_bp) / obp) if obp else None
    return {"overlap_bp": obp, "overlap": ov, "n_snv": n_snv, "n_indel": n_indel, "indel_bp": indel_bp,
            "identity": ident, "diffs": diffs}


def relation(overlap_bp: int, n_diff: int, min_overlap: int) -> str:
    if overlap_bp < min_overlap:
        return "disjoint"
    return "overlapping_identical" if n_diff == 0 else "overlapping_divergent"


def cluster_members(names: list[str], pairs: list[dict], threshold: float, min_overlap: int,
                    linkable: set[str]) -> list[list[str]]:
    """Single-linkage clusters of one reference's contigs (in the order of `names`) on pairs of linkable contigs
    with overlap >= min_overlap and overlap identity >= threshold."""
    from .panel import single_linkage
    edges = [(p["contig_a"], p["contig_b"]) for p in pairs
             if p["contig_a"] in linkable and p["contig_b"] in linkable and p["overlap_bp"] >= min_overlap
             and p["identity"] is not None and p["identity"] >= threshold]
    num = single_linkage(names, edges)
    groups: dict[int, list[str]] = defaultdict(list)
    for n in names:
        groups[num[n]].append(n)
    return [groups[k] for k in sorted(groups)]


def consensus_stats(members: list[Track], ref_len: int) -> tuple[float | None, float]:
    """Identity (%) of the cluster consensus to the reference over the union of its members' intervals, and the
    covered fraction of the reference. At each differing site the allele of the majority of the covering members
    is the consensus (a tie keeps the reference)."""
    union = merge_intervals([iv for m in members for iv in m.intervals])
    ulen = interval_length(union)
    if not ulen or not ref_len:
        return None, 0.0
    cols = 0
    for site in {s for m in members for s in m.alleles}:
        pos, kind = site
        votes: Counter = Counter()
        for m in members:
            if site in m.alleles:
                votes[m.alleles[site]] += 1
            elif _in_any(pos, m.intervals) and not _in_any(pos, m.deleted):
                votes["<ref>"] += 1
        if not votes:
            continue
        top, n = votes.most_common(1)[0]
        if top != "<ref>" and n > votes.get("<ref>", 0):
            cols += 1 if kind == "snv" else _allele_len(top)
    return 100.0 * (1 - cols / ulen), ulen / ref_len


def distinguishing_sites(cluster: set[str], tracks: dict[str, Track], pairs: list[dict],
                         min_overlap: int) -> tuple[dict, list[str]]:
    """Sites where an unflagged member of `cluster` differs from an unflagged contig of ANOTHER cluster inside an
    overlap >= min_overlap: {site: the cluster's allele}, and the ids of those divergent pairs."""
    alleles: dict = defaultdict(Counter)
    used = []
    for p in pairs:
        a, b = p["contig_a"], p["contig_b"]
        if (a in cluster) == (b in cluster) or not p["diffs"] or p["overlap_bp"] < min_overlap:
            continue
        if not (tracks[a].linkable and tracks[b].linkable):
            continue
        used.append(p["pair_id"])
        for site, x, y in p["diffs"]:
            alleles[site][x if a in cluster else y] += 1
    return {s: c.most_common(1)[0][0] for s, c in alleles.items()}, used


def haplotype_carries(dist: dict, block_sites: list, pattern: tuple) -> bool:
    """A local haplotype carries a cluster's variants: >= 2 of the cluster's distinguishing sites lie in the block
    and the haplotype has the cluster's allele at every one of them."""
    hits = [(dist[s], a) for s, a in zip(block_sites, pattern) if s in dist]
    return len(hits) >= 2 and all(x == y for x, y in hits)


def classify(flagged: bool, candidate: bool, phased: bool, reference_supported: bool) -> str:
    """Evidence grade of one cluster (rule printed in viral_strain_clusters.tsv)."""
    if flagged:
        return CLASSES[0]
    if reference_supported:
        return CLASSES[3]
    if candidate and phased:
        return CLASSES[2]
    return CLASSES[1] if candidate else CLASSES[0]


# ---------------------------------------------------------------- statistics (pure, unit-tested)

def binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p), in log space (no scipy)."""
    if k <= 0:
        return 1.0
    if k > n or p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    lp, lq, ln = math.log(p), math.log1p(-p), math.lgamma(n + 1)

    def pmf(i: int) -> float:
        return math.exp(ln - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lp + (n - i) * lq)
    mode = int((n + 1) * p)
    if k <= mode:
        return max(0.0, min(1.0, 1.0 - sum(pmf(i) for i in range(k))))
    total = 0.0
    for i in range(k, n + 1):
        t = pmf(i)
        total += t
        if t < total * 1e-16:
            break
    return min(1.0, total)


def ld_stats(n11: int, n12: int, n21: int, n22: int) -> tuple[float | None, float | None]:
    """r2 and D' of a 2x2 table (rows: site 1 major / other; columns: site 2 major / other); None when a site is
    monomorphic among these fragments."""
    n = n11 + n12 + n21 + n22
    if not n:
        return None, None
    pa, pb = (n11 + n12) / n, (n11 + n21) / n
    den = pa * (1 - pa) * pb * (1 - pb)
    if den <= 0:
        return None, None
    d = n11 / n - pa * pb
    dmax = min(pa * (1 - pb), (1 - pa) * pb) if d > 0 else min(pa * pb, (1 - pa) * (1 - pb))
    return d * d / den, (d / dmax if dmax > 0 else None)


def possible_errors(patterns: dict, error_rate: float, alpha: float = ERROR_ALPHA) -> dict:
    """{pattern: the more frequent pattern at Hamming distance 1 that could produce its count by sequencing error
    (binomial P(X >= count | count_parent + count, error_rate) >= alpha), or None}."""
    ranked = sorted(patterns.items(), key=lambda kv: -kv[1])
    out = {}
    for p, c in patterns.items():
        out[p] = None
        for q, cq in ranked:
            if cq <= c:
                break
            if sum(x != y for x, y in zip(p, q)) == 1 and binom_sf(c, cq + c, error_rate) >= alpha:
                out[p] = q
                break
    return out


def find_blocks(n_sites: int, frags: list[dict], min_link: int) -> list[list[int]]:
    """Local blocks over sites 0..n_sites-1 (in position order): from left to right, a block grows by the next
    site while >= min_link fragments carry an allele at every site of it; blocks of one site are not kept and
    blocks never share a site."""
    carriers: list[set] = [set() for _ in range(n_sites)]
    for k, f in enumerate(frags):
        for i in f:
            carriers[i].add(k)
    blocks, i = [], 0
    while i < n_sites:
        cur, j = carriers[i], i
        while j + 1 < n_sites:
            nxt = cur & carriers[j + 1]
            if len(nxt) < min_link:
                break
            cur, j = nxt, j + 1
        if j > i:
            blocks.append(list(range(i, j + 1)))
            i = j + 1
        else:
            i += 1
    return blocks


def block_haplotypes(block: list[int], frags: list[dict], error_rate: float, min_hap_reads: int,
                     rng: np.random.Generator, n_boot: int = BOOTSTRAPS) -> list[dict]:
    """Complete allele patterns of the fragments that carry every site of the block, most frequent first:
    count, frequency, percentile bootstrap interval (multinomial over those fragments), possible_error parent
    and `supported`."""
    pats = Counter(tuple(f[i] for i in block) for f in frags if all(i in f for i in block))
    if not pats:
        return []
    order = sorted(pats, key=lambda p: (-pats[p], p))
    counts = np.array([pats[p] for p in order], dtype=float)
    n = int(counts.sum())
    draws = rng.multinomial(n, counts / n, size=n_boot) / n
    lo, hi = np.percentile(draws, [2.5, 97.5], axis=0)
    err = possible_errors(pats, error_rate)
    return [{"pattern": p, "count": pats[p], "frequency": pats[p] / n, "ci_low": float(lo[k]), "ci_high": float(hi[k]),
             "error_parent": err[p], "supported": pats[p] >= min_hap_reads and err[p] is None, "n": n}
            for k, p in enumerate(order)]


# ---------------------------------------------------------------- reads: CIGAR walk and variant calling

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
BASE_INDEX = np.full(256, 4, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    BASE_INDEX[ord(_b)] = _i
    BASE_INDEX[ord(_b.lower())] = _i


def walk_cigar(pos0: int, cigar: str) -> list[tuple[str, int, int, int]]:
    """(op, reference position, query position, length) of the M (=, X), I and D operations of a CIGAR."""
    out, r, q = [], pos0, 0
    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op in "M=X":
            out.append(("M", r, q, n))
            r += n
            q += n
        elif op == "I":
            out.append(("I", r, q, n))
            q += n
        elif op == "D":
            out.append(("D", r, q, n))
            r += n
        elif op == "N":
            r += n
        elif op == "S":
            q += n
    return out


def _quals(qual: str, n: int) -> np.ndarray:
    if qual == "*":
        return np.full(n, 40, dtype=np.int16)
    return np.frombuffer(qual.encode(), dtype=np.uint8).astype(np.int16) - 33


def _indel_events(blocks, seq: str, q: np.ndarray) -> list[tuple[int, str, int]]:
    """(reference position, allele, quality) of the indels of a read; the quality of an insertion is its lowest
    base quality, of a deletion the lower of its two flanking bases."""
    out = []
    for op, r, qi, n in blocks:
        if op == "I":
            out.append((r, "+" + seq[qi:qi + n].upper(), int(q[qi:qi + n].min()) if n else 0))
        elif op == "D":
            fl = [int(q[x]) for x in (qi - 1, qi) if 0 <= x < len(q)]
            out.append((r, f"-{n}", min(fl) if fl else 0))
    return out


def read_site_alleles(pos0: int, cigar: str, seq: str, qual: str, site_pos: list[int], site_kind: list[str],
                      site_alleles: list, min_baseq: int) -> dict[int, str]:
    """{site index: allele} of one read at the given sites (sorted positions). An SNV site needs an aligned base
    with quality >= min_baseq; an indel site gets the read's indel there, or '.' when the read aligns the bases
    on both sides of the site without an indel. Alleles not among the site's alleles are left out."""
    blocks = walk_cigar(pos0, cigar)
    if not blocks:
        return {}
    q = _quals(qual, len(seq))
    end = max(r + (n if op in "MD" else 0) for op, r, _qi, n in blocks)
    lo, hi = bisect_left(site_pos, pos0), bisect_right(site_pos, end)
    if lo >= hi:
        return {}
    mblocks = [(r, r + n, qi) for op, r, qi, n in blocks if op == "M"]
    events = {r: (a, qv) for r, a, qv in _indel_events(blocks, seq, q)}

    def qpos(p):
        for s, e, qi in mblocks:
            if s <= p < e:
                return qi + p - s
        return None
    out = {}
    for i in range(lo, hi):
        p = site_pos[i]
        if site_kind[i] == "snv":
            qi = qpos(p)
            if qi is None or q[qi] < min_baseq:
                continue
            a = seq[qi].upper()
        elif p in events:
            a, qv = events[p]
            if qv < min_baseq:
                continue
        else:
            q1, q0 = qpos(p), qpos(p - 1)
            if q1 is None or q0 is None or q1 != q0 + 1:
                continue
            a = "."
        if a in site_alleles[i]:
            out[i] = a
    return out


def sam_fragments(path: Path):
    """Fragments of a saved Bowtie 2 SAM (records of one read name are adjacent): lists of parsed records
    (flag, reference, 0-based position, cigar, tlen, seq, qual) without unmapped, secondary, supplementary or
    QC-failed records."""
    name, recs = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            if f[0] != name:
                if recs:
                    yield recs
                name, recs = f[0], []
            flag = int(f[1])
            if flag & 0xB04:
                continue
            recs.append((flag, f[2], int(f[3]) - 1, f[5], int(f[8]), f[9], f[10]))
    if recs:
        yield recs


def _mate_order(recs):
    """Mate 1 first; mate 2's positions inside mate 1's reference span are not counted again (one molecule)."""
    return sorted(recs, key=lambda r: 0 if r[0] & 64 or not r[0] & 128 else 1)


def _span(rec) -> tuple[str, int, int]:
    blocks = walk_cigar(rec[2], rec[3])
    end = max((r + n for op, r, _q, n in blocks if op in "MD"), default=rec[2])
    return rec[1], rec[2], end


@dataclass
class Pileup:
    counts: dict                 # ref -> int32 array (L, 4 bases, 2 strands), quality-filtered bases
    indels: dict                 # ref -> {(pos, allele): [forward, reverse]}
    read_lengths: Counter = field(default_factory=Counter)
    inserts: dict = field(default_factory=lambda: defaultdict(Counter))       # ref -> Counter of |TLEN|
    fragments: Counter = field(default_factory=Counter)                      # ref -> fragments aligned


def pileup(sam: Path, ref_seqs: dict[str, str], min_baseq: int) -> Pileup:
    """Pass 1 over the saved SAM: base and indel counts per strand (one observation per fragment and position)."""
    pu = Pileup({r: np.zeros((len(s), 4, 2), dtype=np.int32) for r, s in ref_seqs.items()},
                {r: defaultdict(lambda: [0, 0]) for r in ref_seqs})
    for recs in sam_fragments(sam):
        recs = _mate_order(recs)
        mask = None
        seen_refs = set()
        for k, (flag, ref, pos0, cigar, tlen, seq, qual) in enumerate(recs):
            if ref not in pu.counts or seq == "*":
                continue
            seen_refs.add(ref)
            pu.read_lengths[len(seq)] += 1
            if flag & 2 and flag & 64 and tlen:
                pu.inserts[ref][abs(tlen)] += 1
            strand = 1 if flag & 16 else 0
            blocks = walk_cigar(pos0, cigar)
            q = _quals(qual, len(seq))
            sb = BASE_INDEX[np.frombuffer(seq.encode(), dtype=np.uint8)]
            ms = me = -1
            if mask is not None and mask[0] == ref:
                ms, me = mask[1], mask[2]
            for op, r, qi, n in blocks:
                if op != "M":
                    continue
                pos = np.arange(r, r + n)
                b = sb[qi:qi + n]
                keep = (b < 4) & (q[qi:qi + n] >= min_baseq) & ~((pos >= ms) & (pos < me))
                np.add.at(pu.counts[ref], (pos[keep], b[keep], strand), 1)
            for r, a, qv in _indel_events(blocks, seq, q):
                if qv >= min_baseq and not (ms <= r < me):
                    pu.indels[ref][(r, a)][strand] += 1
            if k == 0:
                mask = _span(recs[0])
        for ref in seen_refs:
            pu.fragments[ref] += 1
    return pu


def call_sites(counts: np.ndarray, indels: dict, ref_seq: str, min_depth: int, min_alt: int, min_af: float,
               error_rate: float, alpha: float = SIG_ALPHA) -> tuple[list[dict], int]:
    """Variant sites of one sample and reference from pileup counts. Every non-major allele at a position with
    depth >= min_depth is one test (m tests: 3 per SNV position, one per indel allele); it passes with
    >= min_alt reads, AF >= min_af, reads on both strands and binomial P(X >= alt | depth, error_rate) < alpha/m.
    Returns (sites, m). status: polymorphic (a passing non-major allele) or fixed (major allele != reference)."""
    tot = counts.sum(axis=2)
    depth = tot.sum(axis=1)
    major = tot.argmax(axis=1)
    tested = depth >= min_depth
    by_pos: dict[int, dict] = defaultdict(dict)
    for (p, a), fr in indels.items():
        by_pos[p][a] = fr
    indel_tables = {}
    for p, alle in by_pos.items():
        if p >= len(depth):
            continue
        n_ins = sum(sum(fr) for a, fr in alle.items() if a.startswith("+"))
        n_del = sum(sum(fr) for a, fr in alle.items() if a.startswith("-"))
        d = int(depth[p]) + n_del
        if d < min_depth:
            continue
        tab = {a: list(fr) for a, fr in alle.items()}
        ref_fr = counts[p].sum(axis=0)
        tab["."] = [max(0, int(ref_fr[0]) - sum(fr[0] for a, fr in alle.items() if a.startswith("+"))),
                    max(0, int(ref_fr[1]) - sum(fr[1] for a, fr in alle.items() if a.startswith("+")))]
        if sum(tab["."]) == 0 and d - n_ins - n_del > 0:
            tab["."] = [d - n_ins - n_del, 0]
        indel_tables[p] = (d, tab)
    m = 3 * int(tested.sum()) + sum(len(t) - 1 for _d, t in indel_tables.values())
    thr = alpha / max(m, 1)
    sites = []

    def assess(pos, kind, ref_allele, d, tab):
        cnt = {a: sum(fr) for a, fr in tab.items()}
        maj = max(cnt, key=lambda a: (cnt[a], a == ref_allele))
        passing, best_p = [], None
        for a, fr in tab.items():
            c = cnt[a]
            if a == maj or c < min_alt or c < min_af * d or min(fr) < 1:
                continue
            pv = binom_sf(c, d, error_rate)
            best_p = pv if best_p is None else min(best_p, pv)
            if pv < thr:
                passing.append(a)
        if not passing and maj == ref_allele:
            return
        alle = [maj] + sorted(passing, key=lambda a: -cnt[a])
        sites.append({"pos": pos, "kind": kind, "ref_allele": ref_allele, "major": maj, "alleles": alle,
                      "depth": d, "allele_counts": ",".join(f"{a}:{cnt[a]}" for a in sorted(cnt, key=lambda a: -cnt[a])
                                                            if cnt[a]),
                      "minor_af": (cnt[passing[0]] / d) if passing else None,
                      "p_value": best_p, "status": "polymorphic" if passing else "fixed"})
    refidx = BASE_INDEX[np.frombuffer(ref_seq.encode(), dtype=np.uint8)]
    cand = set(np.nonzero(tested & (refidx < 4) & (major != refidx))[0].tolist())      # fixed differences (not at N)
    for b in range(4):
        a = tot[:, b]
        ok = tested & (major != b) & (a >= min_alt) & (a >= min_af * depth) & (counts[:, b, 0] >= 1) & \
            (counts[:, b, 1] >= 1)
        cand.update(np.nonzero(ok)[0].tolist())
    for p in sorted(int(x) for x in cand):
        tab = {"ACGT"[b]: [int(counts[p, b, 0]), int(counts[p, b, 1])] for b in range(4)}
        assess(p, "snv", ref_seq[p].upper(), int(depth[p]), tab)
    for p in sorted(indel_tables):
        d, tab = indel_tables[p]
        assess(p, "indel", ".", d, tab)
    sites.sort(key=lambda s: (s["pos"], s["kind"]))
    return sites, m


def fragment_alleles(sam: Path, sites_by_ref: dict[str, list[dict]], min_baseq: int) -> dict[str, list[dict]]:
    """Pass 2: per reference, the {site index: allele} of every fragment with alleles at >= 2 polymorphic sites
    (mates merged; a site where the mates disagree is dropped for that fragment)."""
    idx = {r: ([s["pos"] for s in ss], [s["kind"] for s in ss], [set(s["alleles"]) for s in ss])
           for r, ss in sites_by_ref.items() if ss}
    out: dict[str, list[dict]] = {r: [] for r in idx}
    for recs in sam_fragments(sam):
        per_ref: dict[str, dict] = defaultdict(dict)
        bad: dict[str, set] = defaultdict(set)
        for flag, ref, pos0, cigar, _tlen, seq, qual in recs:
            if ref not in idx or seq == "*":
                continue
            got = read_site_alleles(pos0, cigar, seq, qual, *idx[ref], min_baseq)
            cur = per_ref[ref]
            for i, a in got.items():
                if i in cur and cur[i] != a:
                    bad[ref].add(i)
                cur[i] = a
        for ref, cur in per_ref.items():
            for i in bad[ref]:
                cur.pop(i, None)
            if len(cur) >= 2:
                out[ref].append(cur)
    return out


def linkage_rows(sites: list[dict], frags: list[dict], min_link: int) -> list[dict]:
    """2x2 linkage of every pair of sites carried together by >= min_link fragments (major vs other allele)."""
    tabs: dict[tuple, list] = defaultdict(lambda: [0, 0, 0, 0])
    for f in frags:
        for i, j in combinations(sorted(f), 2):
            x, y = f[i] == sites[i]["major"], f[j] == sites[j]["major"]
            tabs[(i, j)][(0 if x else 2) + (0 if y else 1)] += 1
    rows = []
    for (i, j), t in sorted(tabs.items()):
        if sum(t) < min_link:
            continue
        r2, dp = ld_stats(*t)
        rows.append({"site_a": _site_label(sites[i]), "site_b": _site_label(sites[j]), "n_fragments": sum(t),
                     "n_major_major": t[0], "n_major_other": t[1], "n_other_major": t[2], "n_other_other": t[3],
                     "r2": _round(r2, 4), "d_prime": _round(dp, 4)})
    return rows


# ---------------------------------------------------------------- helpers

def _round(x, k=3):
    return "" if x is None else round(float(x), k)


def _site_label(s: dict) -> str:
    return f"{s['pos'] + 1}:{s['kind']}"


def _site_key_label(site, allele=None) -> str:
    return f"{site[0] + 1}:{site[1]}" + (f":{allele}" if allele is not None else "")


def _slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9.]+", "_", str(text)).strip("_")


def resolution_tag(threshold: float) -> str:
    return f"c{threshold:g}"


def parse_identities(text: str) -> list[float]:
    try:
        vals = sorted({float(x) for x in str(text).split(",") if x.strip()})
    except ValueError:
        raise ValueError(f"--identity {text!r}: give comma-separated percentages, e.g. 97,98,99") from None
    if not vals or any(not 50 < v <= 100 for v in vals):
        raise ValueError(f"--identity {text!r}: percentages above 50 and at most 100")
    return vals


def extract_sequences(fasta: str | Path, wanted: set[str]) -> dict[str, str]:
    """The records of `fasta` whose first header word is in `wanted` (streamed; other records are not kept)."""
    out, cur, chunks = {}, None, []
    with open_text(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(chunks).upper()
                sid = line[1:].split()[0] if line[1:].strip() else ""
                cur, chunks = (sid if sid in wanted and sid not in out else None), []
            elif cur is not None:
                chunks.append(line.strip())
    if cur is not None:
        out[cur] = "".join(chunks).upper()
    return out


def _write_tsv(df: pd.DataFrame, path: Path, notes=()) -> None:
    with open(path, "w") as fo:
        for n in notes:
            fo.write(f"# {n}\n")
        df.to_csv(fo, sep="\t", index=False)


def read_commented_tsv(path: str | Path) -> pd.DataFrame:
    """A TSV whose leading '#' lines document its rules."""
    import io
    lines = [l for l in open_text(path)]
    k = 0
    while k < len(lines) and lines[k].startswith("#"):
        k += 1
    return pd.read_csv(io.StringIO("".join(lines[k:])), sep="\t", dtype=str, keep_default_na=False)


def read_strain_clusters(path: str | Path, resolution: float) -> tuple[dict[str, str], dict]:
    """contig -> cluster id at `resolution` from a sealed viral_strain_clusters.tsv, and the manifest record of
    its strain-evidence folder (panel-build --strain-clusters)."""
    rec = manifest.verify_file(path, what="strain-evidence outputs")
    df = read_commented_tsv(path)
    need = {"resolution", "cluster_id", "members"}
    if not need <= set(df.columns):
        miss = ", ".join(sorted(need - set(df.columns)))
        raise ValueError(f"{path}: not a viral_strain_clusters.tsv (missing {miss})")
    res = pd.to_numeric(df["resolution"], errors="coerce")
    sel = df[(res - float(resolution)).abs() < 1e-9]
    if sel.empty:
        have = ", ".join(f"{x:g}" for x in sorted(res.dropna().unique()))
        raise ValueError(f"{path} has no clusters at --strain-resolution {resolution:g} (resolutions: {have or '-'})")
    out: dict[str, str] = {}
    for r in sel.itertuples(index=False):
        for c in (m.strip() for m in str(r.members).split(",")):
            if not c:
                continue
            if c in out and out[c] != r.cluster_id:
                raise ValueError(f"{path}: contig {c} is in two clusters at {resolution:g} ({out[c]}, {r.cluster_id})")
            out[c] = r.cluster_id
    return out, rec


def _tool_line(cmd: list[str]) -> str:
    import subprocess
    try:
        p = subprocess.run(cmd, capture_output=True, text=True)
        lines = (p.stdout or p.stderr).strip().splitlines()
        return lines[0].split("/")[-1] if lines else "unknown"           # bowtie2 prints its binary path
    except OSError:
        return "not found"


# ---------------------------------------------------------------- figures

def _track_figure(ref: str, ref_len: int, tracks: list[Track], cluster_of: dict[str, str], res: float,
                  path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = sorted(tracks, key=lambda t: (cluster_of.get(t.contig, ""), t.intervals[0][0] if t.intervals else 0))
    if len(rows) > MAX_TRACK_ROWS:
        log.info("overlap tracks of %s: the first %d of %d contigs are drawn", ref, MAX_TRACK_ROWS, len(rows))
        rows = rows[:MAX_TRACK_ROWS]
    cids = sorted({cluster_of.get(t.contig, "") for t in rows})
    cmap = plt.get_cmap("tab10")
    col = {c: cmap(k % 10) for k, c in enumerate(cids)}
    fig, ax = plt.subplots(figsize=(10, 1.2 + 0.32 * len(rows)))
    for y, t in enumerate(rows):
        c = "#9E9E9E" if t.flags else col[cluster_of.get(t.contig, "")]
        for s, e in t.intervals:
            ax.plot([s, e], [y, y], lw=6, color=c, solid_capstyle="butt")
        xs = [site[0] for site in t.alleles]
        if xs:
            ax.vlines(xs, y - 0.3, y + 0.3, color="black", lw=0.6)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{t.contig} [{cluster_of.get(t.contig, '-')}]" + (" *" if t.flags else "") for t in rows],
                       fontsize=7)
    ax.set_xlim(0, ref_len)
    ax.set_ylim(-1, len(rows))
    ax.invert_yaxis()
    ax.set_xlabel(f"{ref} position (bp)")
    ax.set_title(f"Contigs on {ref}, coloured by {res:g} % cluster (ticks: differences from the reference; "
                 "grey *: flagged)", fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _linkage_figure(block: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    haps = block["haplotypes"]
    fig, ax = plt.subplots(figsize=(7, 1.0 + 0.35 * len(haps)))
    labels = ["|".join(h["pattern"]) for h in haps]
    colors = ["#3C5488" if h["supported"] else ("#B3B3B3" if h["error_parent"] else "#91D1C2") for h in haps]
    ax.barh(range(len(haps)), [h["count"] for h in haps], color=colors)
    ax.set_yticks(range(len(haps)))
    ax.set_yticklabels(labels, fontsize=7, family="monospace")
    ax.invert_yaxis()
    ax.set_xlabel("fragments carrying the pattern")
    ax.set_title(f"{block['block_id']}  sites {', '.join(_site_key_label(s) for s in block['sites'])}\n"
                 "dark: supported; grey: possible sequencing error; light: below --min-hap-reads", fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- the command

@dataclass(frozen=True)
class Params:
    identities: tuple = (97.0, 98.0, 99.0)
    min_overlap: int = 200
    categories: tuple = ("virus",)
    min_depth: int = 10
    min_alt: int = 5
    min_af: float = 0.02
    min_baseq: int = 20
    error_rate: float = 0.005
    min_link: int = 5
    min_hap_reads: int = 3
    reference_identity: float = 99.0
    reference_cover: float = 0.8
    max_fragments: int | None = None
    seed: int = 1

    def check(self) -> None:
        if self.min_overlap < 1:
            raise ValueError("--min-overlap must be at least 1")
        if not 0 < self.error_rate < 1 or not 0 <= self.min_af <= 1:
            raise ValueError("--error-rate must lie in (0, 1) and --min-af in [0, 1]")
        if min(self.min_depth, self.min_alt, self.min_link, self.min_hap_reads) < 1:
            raise ValueError("--min-depth, --min-alt, --min-link and --min-hap-reads must be at least 1")
        if not 0 < self.reference_cover <= 1 or not 0 < self.reference_identity <= 100:
            raise ValueError("--reference-cover is a fraction in (0, 1], --reference-identity a percentage")
        if self.max_fragments is not None and self.max_fragments < 1:
            raise ValueError("--max-fragments must be at least 1")
        if self.seed < 0:
            raise ValueError("--seed must be 0 or a positive whole number")


def _run_references(run_dir: Path, reference: str | None) -> tuple[list[Path], dict]:
    """The FASTA files reference genomes are pulled from: --reference (or the run's virus genomes), then the
    run's extra virus FASTA; the run's own files must still have the content the run recorded."""
    from .panel import check_references_against_run
    cfg = {}
    if (run_dir / "config.used.yaml").exists():
        cfg = (yaml.safe_load((run_dir / "config.used.yaml").read_text()) or {}).get("references", {}) or {}
    recorded = json.loads((run_dir / "references.used.json").read_text())
    virus, extra = reference or cfg.get("virus"), cfg.get("extra_virus_fasta")
    check_references_against_run(recorded.get("references") or {}, None if reference else virus, None, extra, [])
    if reference:
        log.info("reference genomes from --reference %s (not checked against the run's record)", reference)
    files = [Path(f) for f in (virus, extra) if f]
    for f in files:
        if not f.is_file():
            raise FileNotFoundError(f"reference FASTA {f} not found (pass --reference FASTA)")
    return files, {"virus": str(virus) if virus else None, "extra_virus_fasta": str(extra) if extra else None}


def strain_evidence(run_dir: str | Path, outdir: str | Path, params: Params, reference: str | None = None,
                    samples: list | None = None, threads: int = 4) -> Path:
    """Run the whole step (module docstring) and seal `outdir` with MANIFEST.json last."""
    params.check()
    require("minimap2")
    if samples:
        require("bowtie2", "bowtie2-build")
        from .config import check_sample_names
        check_sample_names([s.name for s in samples], "strain-evidence samples")
    run_dir, outdir = Path(run_dir), Path(outdir)
    ctab, cfa = run_dir / "tables" / "contig_category.tsv", run_dir / "contigs" / "contigs.fa"
    rec_run = manifest.verify(run_dir, required=("contigs/contigs.fa", "tables/contig_category.tsv",
                                                 "references.used.json"), what="run outputs")
    ref_files, ref_used = _run_references(run_dir, reference)
    read_inputs = [x for s in (samples or []) for x in (s.r1, s.r2) if x]
    refuse_overlap([outdir / n for n in OUTPUT_TABLES + (manifest.NAME,)], [ctab, cfa, *ref_files, *read_inputs],
                   "strain-evidence")

    contigs = pd.read_csv(ctab, sep="\t", dtype={"best_hit": str, "best_hit_title": str, "method": str, "layer": str,
                                                 "contig": str})
    sel = contigs[contigs["category"].isin(params.categories)]
    if sel.empty:
        raise ValueError(f"no contigs in categories {sorted(params.categories)}; categories present: "
                         f"{sorted(contigs['category'].astype(str).unique())}")
    status, best = [], {}
    for r in sel.itertuples(index=False):
        hit = r.best_hit if isinstance(r.best_hit, str) and r.best_hit not in ("", "nan") else ""
        method = r.method if isinstance(r.method, str) else "-"
        if hit and method in NT_METHODS:
            best[r.contig] = hit
        else:
            why = "no reference hit" if not hit else f"{method} hit (not placed on a genome)"
            status.append({"contig": r.contig, "length": int(r.length), "category": r.category, "method": method,
                           "reference": hit, "status": f"not_analysable: {why}", "n_segments": 0, "ref_start": "",
                           "ref_end": "", "aligned_fraction": "", "flags": "", "other_references": ""})
    if not best:
        raise ValueError(f"no contig of categories {sorted(params.categories)} has a BLASTN reference hit; nothing "
                         "to place on a genome")
    wanted = set(best.values())
    refs: dict[str, str] = {}
    for f in ref_files:
        refs.update({k: v for k, v in extract_sequences(f, wanted - set(refs)).items() if v})
    missing = sorted(wanted - set(refs))
    if missing:
        raise ValueError(f"{len(missing)} reference genome(s) named by the contigs were not found in "
                         f"{', '.join(map(str, ref_files))}: {', '.join(missing[:5])}; pass the FASTA the run "
                         "classified with (--reference)")

    previous = manifest.previous_artifacts(outdir)
    manifest.begin(outdir, "strain-evidence")
    stage = make_staging_dir(outdir, "strain-evidence")
    work = stage / "work"
    work.mkdir()
    (stage / "figures").mkdir()
    try:
        written = _compute(stage, work, contigs, cfa, best, refs, status, params, samples or [], threads)
        (outdir / "figures").mkdir(exist_ok=True)
        for rel in written:
            os.replace(stage / rel, outdir / rel)
    finally:
        shutil.rmtree(work, ignore_errors=True)
        remove_staging_dir(stage)
    manifest.remove_stale(outdir, previous, written)
    record = {
        "step": "strain-evidence", "unclassified_meta_version": __version__,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "parameters": {**{k: (list(v) if isinstance(v, tuple) else v) for k, v in params.__dict__.items()},
                       "window_bp": WINDOW, "identity_jump": IDENTITY_JUMP,
                       "min_aligned_fraction": MIN_ALIGNED_FRACTION, "multi_reference_min_bp": MULTI_REF_MIN_BP,
                       "bowtie2_max_insert": MAX_INSERT, "variant_alpha_bonferroni": SIG_ALPHA,
                       "possible_error_alpha": ERROR_ALPHA, "bootstraps": BOOTSTRAPS,
                       "minimap2": "-c --cs -x asm20 --secondary=no",
                       "bowtie2": f"--no-unal -X {MAX_INSERT}" + (f" -u {params.max_fragments}"
                                                                   if params.max_fragments else "")},
        "seeds": {"bootstrap": "numpy.random.default_rng([seed, sample_index, reference_index, block_index])",
                  "seed": params.seed},
        "run": {"path": str(run_dir.resolve()), "manifest_sealed": rec_run.get("sealed")},
        "references": ref_used,
        "inputs": {"contigs.fa": {"path": str(cfa.resolve()), "sha256": manifest.sha256(cfa)},
                   "contig_category.tsv": {"path": str(ctab.resolve()), "sha256": manifest.sha256(ctab)},
                   **{f"reference:{f.name}": {"path": str(f.resolve()), "sha256": manifest.sha256(f)}
                      for f in ref_files},
                   **{f"reads:{s.name}": [file_id(s.r1), file_id(s.r2)] for s in (samples or [])}},
        "tools": {"minimap2": _tool_line(["minimap2", "--version"]),
                  **({"bowtie2": _tool_line(["bowtie2", "--version"])} if samples else {})},
    }
    manifest.seal(outdir, written, record)
    log.info("strain-evidence: %d contigs on %d reference(s) -> %s", len(best), len(refs), outdir)
    return outdir


def _compute(stage: Path, work: Path, contigs: pd.DataFrame, cfa: Path, best: dict, refs: dict, status: list,
             params: Params, samples: list, threads: int) -> list[str]:
    """Everything between the checked inputs and the staged outputs; returns the staged files (relative)."""
    lengths = dict(zip(contigs["contig"], contigs["length"].astype(int)))
    q = work / "contigs.fa"
    seqs = {}
    from .io_utils import read_fasta
    for n, s in read_fasta(cfa):
        if n in best:
            seqs[n] = s
    lost = sorted(set(best) - set(seqs))
    if lost:
        raise ValueError(f"{len(lost)} contigs of contig_category.tsv are not in contigs.fa: {', '.join(lost[:5])}")
    write_fasta(((n, seqs[n]) for n in sorted(seqs)), q)
    rfa = work / "references.fa"
    ref_order = sorted(refs)
    write_fasta(((r, refs[r]) for r in ref_order), rfa)
    paf = work / "contigs.paf"
    sh(["minimap2", "-c", "--cs", "-x", "asm20", "--secondary=no", "-t", str(threads), str(rfa), str(q)], stdout=paf)
    by_contig: dict[str, list[Segment]] = defaultdict(list)
    with open(paf) as fh:
        for s in parse_paf(fh):
            by_contig[s.contig].append(s)

    tracks: dict[str, Track] = {}
    for c in sorted(best):
        t = build_track(c, lengths.get(c, len(seqs[c])), best[c], by_contig.get(c, []))
        tracks[c] = t
        status.append({"contig": c, "length": t.length, "category": "", "method": "", "reference": best[c],
                       "status": "analysed", "n_segments": t.n_segments,
                       "ref_start": t.intervals[0][0] + 1 if t.intervals else "",
                       "ref_end": t.intervals[-1][1] if t.intervals else "",
                       "aligned_fraction": round(t.aligned_fraction, 4), "flags": ",".join(t.flags),
                       "other_references": ",".join(t.other_refs)})
    meta = contigs.set_index("contig")
    for row in status:
        if row["status"] == "analysed":
            row["category"], row["method"] = meta.at[row["contig"], "category"], meta.at[row["contig"], "method"]

    # pairs per reference, then clusters per threshold
    by_ref: dict[str, list[str]] = defaultdict(list)
    for c in sorted(tracks, key=lambda c: (-tracks[c].length, c)):
        by_ref[best[c]].append(c)
    pairs_by_ref: dict[str, list[dict]] = defaultdict(list)
    pair_rows, k = [], 0
    for ref in ref_order:
        names = by_ref.get(ref, [])
        for a, b in combinations(sorted(names), 2):
            ta, tb = tracks[a], tracks[b]
            cmp = compare_tracks(ta, tb)
            if not cmp["overlap_bp"]:
                rel = "disjoint"
            else:
                rel = relation(cmp["overlap_bp"], cmp["n_snv"] + cmp["n_indel"], params.min_overlap)
            k += 1
            p = {"pair_id": f"P{k:05d}", "reference": ref, "contig_a": a, "contig_b": b, **cmp, "relation": rel}
            pairs_by_ref[ref].append(p)
            pair_rows.append({
                "pair_id": p["pair_id"], "reference": ref, "contig_a": a, "contig_b": b,
                "a_start": ta.intervals[0][0] + 1 if ta.intervals else "",
                "a_end": ta.intervals[-1][1] if ta.intervals else "",
                "b_start": tb.intervals[0][0] + 1 if tb.intervals else "",
                "b_end": tb.intervals[-1][1] if tb.intervals else "",
                "overlap_bp": cmp["overlap_bp"],
                "overlap_start": cmp["overlap"][0][0] + 1 if cmp["overlap"] else "",
                "overlap_end": cmp["overlap"][-1][1] if cmp["overlap"] else "",
                "n_snv": cmp["n_snv"], "n_indel": cmp["n_indel"], "indel_bp": cmp["indel_bp"],
                "identity": _round(cmp["identity"], 3), "relation": rel,
                "flags_a": ",".join(ta.flags), "flags_b": ",".join(tb.flags),
                "differences": ",".join(f"{_site_key_label(s)}:{x}/{y}" for s, x, y in cmp["diffs"])})

    clusters: list[dict] = []
    for res in params.identities:
        tag, n = resolution_tag(res), 0
        for ref in ref_order:
            names = by_ref.get(ref, [])
            if not names:
                continue
            pairs = pairs_by_ref[ref]
            linkable = {c for c in names if tracks[c].linkable}
            groups = cluster_members(names, pairs, res, params.min_overlap, linkable)
            groups.sort(key=lambda g: (-sum(tracks[c].length for c in g), g[0]))
            for g in groups:
                n += 1
                gs = set(g)
                cflags = []
                near = [p for p in pairs if (p["contig_a"] in gs) != (p["contig_b"] in gs)
                        and p["overlap_bp"] >= params.min_overlap]
                if len(g) == 1 and not near:
                    cflags.append("unresolved_disjoint")
                if any(p["contig_a"] in gs and p["contig_b"] in gs and p["overlap_bp"] >= params.min_overlap
                       and p["identity"] is not None and p["identity"] < res for p in pairs):
                    cflags.append("internal_divergence")
                flagged = any(tracks[c].flags for c in g)
                dist, dpairs = distinguishing_sites(gs, tracks, pairs, params.min_overlap)
                ident, cover = consensus_stats([tracks[c] for c in g], len(refs[ref]))
                clusters.append({"resolution": res, "cluster_id": f"{tag}_{n:03d}", "reference": ref, "members": g,
                                 "flagged": flagged, "cluster_flags": cflags, "dist": dist, "divergent_pairs": dpairs,
                                 "consensus_identity": ident, "reference_cover": cover, "blocks": []})

    # reads: local phasing per sample and reference
    blocks: list[dict] = []
    site_rows, link_rows, limit_rows = [], [], []
    if samples:
        blocks = _phase(work, stage, samples, refs, ref_order, params, threads, site_rows, link_rows, limit_rows)

    # haplotypes that carry a cluster's distinguishing variants
    for cl in clusters:
        if not cl["dist"] or cl["flagged"]:
            continue
        for b in blocks:
            if b["reference"] != cl["reference"]:
                continue
            for k, h in enumerate(b["haplotypes"], 1):
                if h["supported"] and haplotype_carries(cl["dist"], b["sites"], h["pattern"]):
                    cl["blocks"].append(f"{b['block_id']}/H{k}")
                    h.setdefault("clusters", []).append(cl["cluster_id"])

    cl_rows = []
    for cl in clusters:
        candidate = bool(cl["divergent_pairs"]) and not cl["flagged"]
        ident = cl["consensus_identity"]
        ref_ok = ident is not None and ident >= params.reference_identity and \
            cl["reference_cover"] >= params.reference_cover
        cls = classify(cl["flagged"], candidate, bool(cl["blocks"]), ref_ok)
        ev = [f"members={len(cl['members'])}"]
        if cl["divergent_pairs"]:
            ev.append("divergent_pairs=" + ",".join(cl["divergent_pairs"]))
        if cl["blocks"]:
            ev.append("haplotype_blocks=" + ",".join(cl["blocks"]))
        if ident is not None:
            ev.append(f"consensus={ident:.2f}%/{100 * cl['reference_cover']:.1f}% of reference"
                      + (" (meets reference rule)" if ref_ok else ""))
        cl_rows.append({
            "resolution": f"{cl['resolution']:g}", "cluster_id": cl["cluster_id"], "reference": cl["reference"],
            "members": ",".join(cl["members"]), "n_members": len(cl["members"]), "class": cls,
            "evidence": "; ".join(ev), "divergent_pairs": ",".join(cl["divergent_pairs"]),
            "chimera_flags": ";".join(f"{c}:{','.join(tracks[c].flags)}" for c in cl["members"] if tracks[c].flags),
            "cluster_flags": ",".join(cl["cluster_flags"]),
            "consensus_identity": _round(ident, 3), "reference_cover": _round(cl["reference_cover"], 4),
            "distinguishing_sites": ",".join(_site_key_label(s, a) for s, a in sorted(cl["dist"].items())),
            "haplotype_blocks": ",".join(cl["blocks"])})

    written = []

    def put(df, name, notes):
        _write_tsv(df, stage / name, notes)
        written.append(name)

    pair_cols = ["pair_id", "reference", "contig_a", "contig_b", "a_start", "a_end", "b_start", "b_end", "overlap_bp",
                 "overlap_start", "overlap_end", "n_snv", "n_indel", "indel_bp", "identity", "relation", "flags_a",
                 "flags_b", "differences"]
    put(pd.DataFrame(pair_rows, columns=pair_cols), "viral_cluster_overlap.tsv", [
        "pairs of contigs on one reference; positions are 1-based reference coordinates (start..end inclusive)",
        "SNVs / indels are counted between the two contigs INSIDE their overlap only; identity = 100 * (1 - (SNVs + "
        "indel bp) / overlap bp)",
        f"relation: disjoint = overlap < --min-overlap {params.min_overlap} bp (a coverage difference, not "
        "divergence); overlapping_identical / overlapping_divergent otherwise",
        "differences: position:type:allele_a/allele_b (base; +SEQ insertion; -L deletion; . reference)"])
    cl_cols = ["resolution", "cluster_id", "reference", "members", "n_members", "class", "evidence", "divergent_pairs",
               "chimera_flags", "cluster_flags", "consensus_identity", "reference_cover", "distinguishing_sites",
               "haplotype_blocks"]
    put(pd.DataFrame(cl_rows, columns=cl_cols), "viral_strain_clusters.tsv", [
        "operational clusters: single linkage on unflagged contig pairs with overlap >= --min-overlap and overlap "
        "identity >= the resolution; flagged contigs (chimera_flags) stay singletons",
        "class rule: flagged member -> operational_sequence_cluster; else consensus identity >= "
        f"{params.reference_identity:g} % over >= {params.reference_cover:g} of the reference -> "
        "reference_supported_viral_strain; else divergent from another cluster in an overlap (divergent_pairs) and "
        ">= 2 distinguishing sites carried by one supported local read haplotype (haplotype_blocks) -> "
        "phased_viral_haplotype; else divergent_pairs only -> candidate_strain_lineage; else "
        "operational_sequence_cluster",
        "cluster_flags: unresolved_disjoint = a singleton without overlap to any contig (it may belong to any "
        "cluster); internal_divergence = two members overlap below the resolution (chained by single linkage)",
        "distinguishing_sites: position:type:allele of this cluster where it differs from another cluster",
        "a sequence cluster is not a haplotype: see docs/DESIGN_robustness.md section 1"])
    st_cols = ["contig", "length", "category", "method", "reference", "status", "n_segments", "ref_start", "ref_end",
               "aligned_fraction", "flags", "other_references"]
    put(pd.DataFrame(sorted(status, key=lambda r: r["contig"]), columns=st_cols), "contig_status.tsv", [
        "contigs of the chosen categories: analysed (placed on their BLASTN reference) or not_analysable",
        f"flags: split_alignment, non_collinear, multi_reference (>= {MULTI_REF_MIN_BP} bp best on another selected "
        f"reference), identity_jump (>= {IDENTITY_JUMP:g} points between {WINDOW} bp windows), "
        f"low_aligned_fraction (< {MIN_ALIGNED_FRACTION:g}), not_aligned"])
    no_reads = [] if samples else ["no reads were given (--samples / --input-dir): no phasing was done"]
    site_cols = ["sample", "reference", "position", "type", "ref_allele", "major_allele", "alleles", "depth",
                 "allele_counts", "minor_af", "p_value", "status"]
    put(pd.DataFrame(site_rows, columns=site_cols), "variant_sites.tsv", no_reads + [
        f"pileup per fragment (mate overlap counted once), base quality >= {params.min_baseq}; a non-major allele "
        f"passes with depth >= {params.min_depth}, >= {params.min_alt} reads, AF >= {params.min_af:g}, reads on both "
        f"strands and binomial P(X >= alt | depth, {params.error_rate:g}) < {SIG_ALPHA:g} / m (Bonferroni over the m "
        "tested alleles of the sample and reference); fixed = major allele differs from the reference"])
    link_cols = ["sample", "reference", "site_a", "site_b", "n_fragments", "n_major_major", "n_major_other",
                 "n_other_major", "n_other_other", "r2", "d_prime"]
    put(pd.DataFrame(link_rows, columns=link_cols), "variant_linkage.tsv", no_reads + [
        f"pairs of polymorphic sites carried together by >= --min-link {params.min_link} fragments; 2x2 of major "
        "vs other allele"])
    hap_rows = []
    for b in blocks:
        for k, h in enumerate(b["haplotypes"], 1):
            hap_rows.append({
                "sample": b["sample"], "reference": b["reference"], "block_id": b["block_id"],
                "haplotype_id": f"{b['block_id']}/H{k}", "block_start": b["sites"][0][0] + 1,
                "block_end": b["sites"][-1][0] + 1, "sites": ",".join(_site_key_label(s) for s in b["sites"]),
                "n_sites": len(b["sites"]), "n_fragments": h["n"], "pattern": "|".join(h["pattern"]),
                "count": h["count"], "frequency": round(h["frequency"], 4), "ci_low": round(h["ci_low"], 4),
                "ci_high": round(h["ci_high"], 4), "possible_error": h["error_parent"] is not None,
                "error_parent": "|".join(h["error_parent"]) if h["error_parent"] else "",
                "supported": h["supported"], "carries_clusters": ",".join(h.get("clusters", []))})
    hap_cols = ["sample", "reference", "block_id", "haplotype_id", "block_start", "block_end", "sites", "n_sites",
                "n_fragments", "pattern", "count", "frequency", "ci_low", "ci_high", "possible_error", "error_parent",
                "supported", "carries_clusters"]
    put(pd.DataFrame(hap_rows, columns=hap_cols), "viral_haplotype_evidence.tsv", no_reads + [
        "LOCAL haplotypes: blocks are consecutive polymorphic sites carried together by >= --min-link "
        f"{params.min_link} fragments; separate blocks are never joined (no genome-wide haplotype)",
        f"frequency interval: 2.5-97.5 percentiles of {BOOTSTRAPS} multinomial bootstrap draws over the block's "
        f"fragments; possible_error: explainable as error of a more frequent pattern at Hamming distance 1 "
        f"(binomial p >= {ERROR_ALPHA:g}); supported: count >= {params.min_hap_reads} and not possible_error"])
    lim_cols = ["sample", "reference", "reference_length", "n_fragments", "read_length_median", "insert_median",
                "insert_p95", "mean_depth", "n_polymorphic_sites", "n_linked_pairs", "n_blocks",
                "max_block_span", "note"]
    put(pd.DataFrame(limit_rows, columns=lim_cols), "phasing_limits.tsv", no_reads + [
        "what phasing can reach: blocks cannot be longer than a fragment (insert_p95); max_block_span is the "
        "longest block found (first to last site, bp)"])

    # figures
    plot_res = 98.0 if 98.0 in params.identities else params.identities[len(params.identities) // 2]
    cof = {c: cl["cluster_id"] for cl in clusters if cl["resolution"] == plot_res for c in cl["members"]}
    for ref in ref_order:
        if by_ref.get(ref):
            rel = f"figures/overlap_tracks_{_slug(ref)}.png"
            _track_figure(ref, len(refs[ref]), [tracks[c] for c in by_ref[ref]], cof, plot_res, stage / rel)
            written.append(rel)
    drawn = 0
    for b in blocks:
        if sum(h["supported"] for h in b["haplotypes"]) < 2:
            continue
        if drawn >= MAX_LINKAGE_FIGURES:
            log.info("linkage figures: the first %d blocks with >= 2 supported haplotypes are drawn",
                     MAX_LINKAGE_FIGURES)
            break
        rel = f"figures/linkage_{_slug(b['reference'])}_{_slug(b['sample'])}_B{b['number']:03d}.png"
        _linkage_figure(b, stage / rel)
        written.append(rel)
        drawn += 1
    return written


def _phase(work: Path, stage: Path, samples: list, refs: dict, ref_order: list, params: Params, threads: int,
           site_rows: list, link_rows: list, limit_rows: list) -> list[dict]:
    """Bowtie 2 alignment of every sample to the selected references and local phasing (module docstring)."""
    from .fragments import run_bowtie2
    idx = work / "refs_index"
    idx.mkdir()
    write_fasta(((r, refs[r]) for r in ref_order), idx / "refs.fa")
    sh(["bowtie2-build", "--threads", str(threads), str(idx / "refs.fa"), str(idx / "refs")],
       stdout=idx / "build.log")
    opts = f"--no-unal -X {MAX_INSERT}" + (f" -u {params.max_fragments}" if params.max_fragments else "")
    blocks = []
    for si, s in enumerate(samples):
        sam = work / f"{_slug(s.name)}.sam.gz"

        def consume(lines, sam=sam):
            n = 0
            with gzip.open(sam, "wt", compresslevel=1) as fo:
                for line in lines:
                    fo.write(line)
                    n += 1
            return n
        run_bowtie2(idx / "refs", s.r1, s.r2, work / f"{_slug(s.name)}.bowtie2.log", threads, opts, consume)
        pu = pileup(sam, refs, params.min_baseq)
        sites_by_ref = {}
        for ref in ref_order:
            sites, _m = call_sites(pu.counts[ref], pu.indels[ref], refs[ref], params.min_depth, params.min_alt,
                                   params.min_af, params.error_rate)
            for x in sites:
                site_rows.append({"sample": s.name, "reference": ref, "position": x["pos"] + 1, "type": x["kind"],
                                  "ref_allele": x["ref_allele"], "major_allele": x["major"],
                                  "alleles": ",".join(x["alleles"]), "depth": x["depth"],
                                  "allele_counts": x["allele_counts"], "minor_af": _round(x["minor_af"], 4),
                                  "p_value": "" if x["p_value"] is None else f"{x['p_value']:.3g}",
                                  "status": x["status"]})
            sites_by_ref[ref] = [x for x in sites if x["status"] == "polymorphic"]
        frags = fragment_alleles(sam, sites_by_ref, params.min_baseq)
        for ri, ref in enumerate(ref_order):
            ss = sites_by_ref[ref]
            ff = frags.get(ref, [])
            lrows = linkage_rows(ss, ff, params.min_link) if ss else []
            link_rows.extend({"sample": s.name, "reference": ref, **r} for r in lrows)
            found = find_blocks(len(ss), ff, params.min_link) if ss else []
            spans = []
            for bi, blk in enumerate(found, 1):
                rng = np.random.default_rng([params.seed, si, ri, bi])
                haps = block_haplotypes(blk, ff, params.error_rate, params.min_hap_reads, rng)
                sites = [(ss[i]["pos"], ss[i]["kind"]) for i in blk]
                spans.append(sites[-1][0] - sites[0][0] + 1)
                blocks.append({"block_id": f"{s.name}:{ref}:B{bi:03d}", "number": bi, "sample": s.name,
                               "reference": ref, "sites": sites, "haplotypes": haps})
            ins = pu.inserts.get(ref) or Counter()
            ins_vals = np.repeat(np.array(list(ins.keys()), dtype=float), list(ins.values())) if ins else None
            rl = pu.read_lengths
            rl_vals = np.repeat(np.array(list(rl.keys()), dtype=float), list(rl.values())) if rl else None
            depth = pu.counts[ref].sum(axis=(1, 2))
            limit_rows.append({
                "sample": s.name, "reference": ref, "reference_length": len(refs[ref]),
                "n_fragments": pu.fragments.get(ref, 0),
                "read_length_median": _round(np.median(rl_vals), 1) if rl_vals is not None else "",
                "insert_median": _round(np.median(ins_vals), 1) if ins_vals is not None else "",
                "insert_p95": _round(np.percentile(ins_vals, 95), 1) if ins_vals is not None else "",
                "mean_depth": _round(depth.mean(), 2), "n_polymorphic_sites": len(ss), "n_linked_pairs": len(lrows),
                "n_blocks": len(found), "max_block_span": max(spans) if spans else 0,
                "note": "blocks are local (bounded by the fragment length), not genome-wide haplotypes"})
        sam.unlink(missing_ok=True)
    return blocks


# ---------------------------------------------------------------- command line (wired as `strain-evidence`)

def add_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument("--run", required=True, help="output folder of a finished `unclassified-meta run`")
    p.add_argument("-o", "--outdir", default="strain_evidence", help="output folder (default: strain_evidence)")
    p.add_argument("--reference", help="viral genome FASTA (default: the run's virus genomes, then its extra virus "
                                       "FASTA)")
    p.add_argument("--identity", default="97,98,99", help="cluster resolutions, %% identity (default: 97,98,99)")
    p.add_argument("--min-overlap", type=int, default=200, help="bp of shared reference interval for a comparison "
                                                                "(default: 200)")
    p.add_argument("--categories", default="virus", help="contig categories to use (default: virus)")
    g = p.add_argument_group("reads for local phasing (optional; choose one)")
    g.add_argument("--samples", help="sample sheet TSV: sample, r1, [r2]")
    g.add_argument("--input-dir", help="folder with FASTQ files; samples are detected")
    g.add_argument("--max-fragments", type=int, default=None,
                   help="use the first N fragments of each sample (Bowtie 2 -u; default: all)")
    v = p.add_argument_group("variant sites and haplotypes")
    v.add_argument("--min-depth", type=int, default=10)
    v.add_argument("--min-alt", type=int, default=5)
    v.add_argument("--min-af", type=float, default=0.02)
    v.add_argument("--min-baseq", type=int, default=20)
    v.add_argument("--error-rate", type=float, default=0.005, help="per-base error rate of the binomial tests")
    v.add_argument("--min-link", type=int, default=5, help="fragments that must carry every site of a block")
    v.add_argument("--min-hap-reads", type=int, default=3, help="fragments of a supported haplotype")
    v.add_argument("--seed", type=int, default=1, help="seed of the bootstrap intervals")
    c = p.add_argument_group("reference-supported strains")
    c.add_argument("--reference-identity", type=float, default=99.0)
    c.add_argument("--reference-cover", type=float, default=0.8)
    p.add_argument("-t", "--threads", type=int, default=4)


def params_from_args(args) -> Params:
    return Params(identities=tuple(parse_identities(args.identity)), min_overlap=args.min_overlap,
                  categories=tuple(c.strip() for c in args.categories.split(",") if c.strip()),
                  min_depth=args.min_depth, min_alt=args.min_alt, min_af=args.min_af, min_baseq=args.min_baseq,
                  error_rate=args.error_rate, min_link=args.min_link, min_hap_reads=args.min_hap_reads,
                  reference_identity=args.reference_identity, reference_cover=args.reference_cover,
                  max_fragments=args.max_fragments, seed=args.seed)


def main(args) -> int:
    """`strain-evidence`: checks every input against the output folder, locks it, runs and seals."""
    from contextlib import ExitStack

    from .cli import _guard, _inputs, _run_references as run_refs
    from .config import detect_samples, read_samplesheet
    from .runner import setup_logging
    if args.samples and args.input_dir:
        raise SystemExit("give at most one of --samples SHEET.tsv or --input-dir FOLDER")
    params = params_from_args(args)
    params.check()
    samples = read_samplesheet(args.samples) if args.samples else (detect_samples(args.input_dir)
                                                                   if args.input_dir else [])
    with ExitStack() as stack:
        out = _guard(stack, "strain-evidence", args.outdir, _inputs(args, samples, extra=[run_refs(args.run)]))
        setup_logging(out / "strain-evidence.log")
        strain_evidence(args.run, out, params, args.reference, samples, args.threads)
    return 0
