"""Tiny synthetic data set for the end-to-end test (no downloads).

Writes into OUTDIR:
  refs/  the download-refs layout with fixture files (rRNA, RefSeq-like viral genomes and proteins, plant proteins)
  host.fa                    random 'host' genome
  reads/S1_R1.fastq.gz ...   paired samples S1, S2 and single-end S3
  reads/S1.kraken2           Kraken2 per-read output (bacteria + part of host 'classified')
  reads/S2.centrifuge        Centrifuge per-read output
Truth: unclassified reads are mostly virus (virusA in the reference set, virusB only via protein similarity),
some host, a few bacterial reads and random noise.
"""
from __future__ import annotations

import gzip
import random
import sys
from pathlib import Path

BASES = "ACGT"
CODON = {}
for a in BASES:
    for b in BASES:
        for c in BASES:
            CODON[a + b + c] = None
AA = "ACDEFGHIKLMNPQRSTVWY"
TABLE = dict(zip(sorted(CODON),
                 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"))


def rand_seq(rng, n, gc=0.45):
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices(BASES, weights=w, k=n))


def orf_genome(rng, n_codons):
    """Random coding sequence without internal stop codons (so it translates to one long protein)."""
    stops = {"TAA", "TAG", "TGA"}
    out = ["ATG"]
    while len(out) < n_codons:
        c = rand_seq(rng, 3, 0.5)
        if c not in stops:
            out.append(c)
    return "".join(out)


def translate(nt):
    return "".join(TABLE[nt[i:i + 3]] for i in range(0, len(nt) - 2, 3))


def mutate(rng, seq, rate):
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < rate:
            s[i] = rng.choice([b for b in BASES if b != s[i]])
    return "".join(s)


def synonymous_mutate(rng, nt, rate):
    """Change third codon positions only where the amino acid stays the same (keeps the protein, lowers nt identity)."""
    codons = [nt[i:i + 3] for i in range(0, len(nt), 3)]
    for k, c in enumerate(codons):
        if k and rng.random() < rate:
            alts = [c[:2] + b for b in BASES if b != c[2] and TABLE[c[:2] + b] == TABLE[c]]
            if alts:
                codons[k] = rng.choice(alts)
    return "".join(codons)


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def fragments(rng, genome, n, rl=100, ins=(250, 350), err=0.005):
    for _ in range(n):
        L = rng.randint(*ins)
        if L > len(genome):
            L = len(genome)
        st = rng.randint(0, len(genome) - L)
        f = genome[st:st + L]
        if rng.random() < 0.5:
            f = revcomp(f)
        yield mutate(rng, f[:rl], err), mutate(rng, revcomp(f)[:rl], err)


def wfa(path, recs):
    with open(path, "w") as fo:
        for n, s in recs:
            fo.write(f">{n}\n")
            for i in range(0, len(s), 70):
                fo.write(s[i:i + 70] + "\n")


def seal_refs(out: Path) -> None:
    """Seal each reference folder like download-refs does (MANIFEST.json written after the files)."""
    from unclassified_meta import manifest
    for sub in ("rrna", "viral", "protein"):
        d = out / "refs" / sub
        manifest.seal(d, sorted(f.name for f in d.iterdir() if f.is_file() and f.name != manifest.NAME),
                      {"step": "test fixture"})


def main(outdir: str, seed: int = 7) -> None:
    rng = random.Random(seed)
    out = Path(outdir)
    for d in ("refs/rrna", "refs/viral", "refs/protein", "reads"):
        (out / d).mkdir(parents=True, exist_ok=True)

    host = rand_seq(rng, 40000, 0.38)
    bact = rand_seq(rng, 15000, 0.60)
    virus_a = orf_genome(rng, 2000)                       # 6 kb, in the viral reference (as a 90 % variant)
    virus_b_ref = orf_genome(rng, 1500)                   # 4.5 kb; only its protein is in the viral protein set
    virus_b = synonymous_mutate(rng, virus_b_ref, 0.95)   # the sample carries a nucleotide-divergent variant
    plant_gene = orf_genome(rng, 400)                     # a plant gene absent from the host assembly
    rrna = rand_seq(rng, 1800, 0.55)
    noise = rand_seq(rng, 3000, 0.5)

    wfa(out / "host.fa", [("chr1", host)])
    wfa(out / "refs/rrna/plant_rrna.fasta", [("RRNA1 Eukaryota;Archaeplastida;Chloroplastida;Testplant 18S", rrna)])
    wfa(out / "refs/viral/refseq_viral.fna",
        [("NC_900001.1 Test virus A, complete genome", mutate(rng, virus_a, 0.10)),
         ("NC_900002.1 Unrelated virus Z, complete genome", rand_seq(rng, 5000))])
    wfa(out / "refs/viral/refseq_viral_protein.faa",
        [("YP_900002.1 polyprotein [Test virus B]", translate(virus_b_ref).rstrip("*")),
         ("YP_900003.1 capsid [Unrelated virus Z]", "".join(rng.choices(AA, k=300)))])
    wfa(out / "refs/protein/sprot_viridiplantae.faa",
        [("sp|P99999|TEST_PLANT Test plant protein OS=Testplant OX=33090", translate(plant_gene).rstrip("*")),
         ("sp|P99998|OTHER_PLANT Other protein OS=Testplant OX=33090", "".join(rng.choices(AA, k=250)))])
    seal_refs(out)

    # sample composition: (source, genome, pairs)
    designs = {
        "S1": [("host", host, 1500), ("bact", bact, 800), ("virA", virus_a, 3000), ("virB", virus_b, 1500),
               ("plant", plant_gene, 300), ("rrna", rrna, 200), ("noise", noise, 30)],
        "S2": [("host", host, 1000), ("bact", bact, 500), ("virA", virus_a, 4000), ("virB", virus_b, 800),
               ("plant", plant_gene, 200), ("rrna", rrna, 150), ("noise", noise, 30)],
        "S3": [("host", host, 1200), ("bact", bact, 400), ("virA", virus_a, 2500), ("virB", virus_b, 1200),
               ("plant", plant_gene, 200), ("rrna", rrna, 100), ("noise", noise, 30)],
    }
    for sample, parts in designs.items():
        reads = []
        for src, genome, n in parts:
            for i, (a, b) in enumerate(fragments(rng, genome, n)):
                reads.append((f"{sample}_{src}_{i}", a, b, src))
        rng.shuffle(reads)
        q = "I" * 100
        if sample == "S3":   # single-end
            with gzip.open(out / "reads" / f"{sample}.fastq.gz", "wt") as f1:
                for rid, a, _, _ in reads:
                    f1.write(f"@{rid}\n{a}\n+\n{q}\n")
        else:
            with gzip.open(out / "reads" / f"{sample}_R1.fastq.gz", "wt") as f1, \
                    gzip.open(out / "reads" / f"{sample}_R2.fastq.gz", "wt") as f2:
                for rid, a, b, _ in reads:
                    f1.write(f"@{rid}/1\n{a}\n+\n{q}\n")
                    f2.write(f"@{rid}/2\n{b}\n+\n{q}\n")
        # classifier: bacteria always classified; a third of host reads classified; the rest unclassified
        def classified(rid, src):
            return src == "bact" or (src == "host" and int(rid.rsplit("_", 1)[1]) % 3 == 0)
        if sample == "S1":
            with open(out / "reads" / "S1.kraken2", "w") as fo:
                for rid, a, b, src in reads:
                    c = classified(rid, src)
                    fo.write(f"{'C' if c else 'U'}\t{rid}\t{2 if c else 0}\t100|100\t0:66 |:| 0:66\n")
        elif sample == "S2":
            with open(out / "reads" / "S2.centrifuge", "w") as fo:
                fo.write("readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n")
                for rid, a, b, src in reads:
                    if classified(rid, src):
                        # bacterial reads: strong hits; the "host" reads the classifier assigned: weak, short hits
                        sc, hl = (900, 80) if src == "bact" else (49, 22)
                        fo.write(f"{rid}\tNZ_CP000001.1\t562\t{sc}\t0\t{hl}\t200\t1\n")
                    else:
                        fo.write(f"{rid}\tunclassified\t0\t0\t0\t0\t200\t1\n")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 7)
