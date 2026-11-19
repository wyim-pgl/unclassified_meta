"""Synthetic two-group data set for the panel / differential-abundance test (no downloads).

12 paired-end samples (A1-A6, B1-B6), 8 viruses, a host and a bacterium. Virus V4 is six times more abundant in
group B; nothing else differs between the groups. Reference evidence per virus:
  V1-V3  genome in the viral reference at ~98 % identity   -> "X" bins
  V4-V6  genome in the viral reference at ~88 % identity   -> "X-like" bins
  V7     only two proteins in the viral protein reference   -> one "X-like (protein only)" bin
  V8     nothing in any reference                           -> unmatched sequence cluster

The "bacterial" reads are the ones the classifier assigned (Kraken2 per-read output, `C`). For every sample a
Kraken2-style report is written to reports/<sample>.kreport whose totals agree with that per-read output: the
classified reads are spread over 60 species in 12 genera (the same proportions in both groups, so none is
differential), a few stay at genus or domain level, and a small "Viruses" clade stands for viral reads the
classifier did assign. The report is simulated at the table level (the reads themselves are random sequence).
"""
from __future__ import annotations

import gzip
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from make_synthetic import (AA, fragments, seal_refs, mutate, orf_genome, rand_seq, synonymous_mutate, translate,  # noqa: E402
                            wfa)

TRUE_DA = "V4"
FOLD = 6.0
N_GENERA, SPECIES_PER_GENUS = 12, 5


def write_kreport(path: Path, n_total: int, n_classified: int, weights: list[float], rng: random.Random) -> None:
    """Kraken2-style report (pct, clade reads, own reads, rank, taxID, indented name) for one sample."""
    n_virus = round(0.03 * n_classified)
    n_domain = round(0.04 * n_classified)                 # assigned to Bacteria itself (above species)
    n_genus = round(0.06 * n_classified)                  # assigned to a genus (above species)
    n_species = n_classified - n_virus - n_domain - n_genus
    species = [0] * len(weights)
    for i in rng.choices(range(len(weights)), weights=weights, k=n_species):
        species[i] += 1
    genus_own = [0] * N_GENERA
    for g in rng.choices(range(N_GENERA), k=n_genus):
        genus_own[g] += 1
    rows = [(n_total - n_classified, n_total - n_classified, "U", 0, 0, "unclassified"),
            (n_classified, 0, "R", 1, 0, "root"),
            (n_classified - n_virus, n_domain, "D", 2, 1, "Bacteria"),
            (n_classified - n_virus - n_domain, 0, "P", 1224, 2, "Pseudomonadota")]
    for g in range(N_GENERA):
        sp = species[g * SPECIES_PER_GENUS:(g + 1) * SPECIES_PER_GENUS]
        if genus_own[g] + sum(sp) == 0:
            continue
        rows.append((genus_own[g] + sum(sp), genus_own[g], "G", 9000 + g, 3, f"Testbacter{g + 1:02d}"))
        for k, n in enumerate(sp):
            if n:
                rows.append((n, n, "S", 80000 + g * SPECIES_PER_GENUS + k, 4,
                             f"Testbacter{g + 1:02d} species{k + 1}"))
    if n_virus:
        rows += [(n_virus, 0, "D", 10239, 1, "Viruses"), (n_virus, n_virus, "S", 70001, 2, "Test phage 1")]
    with open(path, "w") as fh:
        for clade, own, code, taxid, depth, name in rows:
            fh.write(f"{100 * clade / n_total:6.2f}\t{clade}\t{own}\t{code}\t{taxid}\t{'  ' * depth}{name}\n")


def main(outdir: str, seed: int = 23) -> None:
    rng = random.Random(seed)
    out = Path(outdir)
    rng_rep = random.Random(seed + 1000)                 # reports use their own stream
    sp_weights = [rng_rep.lognormvariate(0, 0.6) for _ in range(N_GENERA * SPECIES_PER_GENUS)]
    for d in ("refs/rrna", "refs/viral", "refs/protein", "reads", "reports"):
        (out / d).mkdir(parents=True, exist_ok=True)
    host = rand_seq(rng, 40000, 0.38)
    bact = rand_seq(rng, 15000, 0.60)
    rrna = rand_seq(rng, 1800, 0.55)
    virus = {f"V{i}": orf_genome(rng, rng.randint(1100, 1500)) for i in range(1, 9)}
    v7_ref = virus["V7"]
    virus["V7"] = synonymous_mutate(rng, v7_ref, 0.95)       # nucleotide-divergent, protein kept
    prot7 = translate(v7_ref).rstrip("*")
    half = len(prot7) // 2

    wfa(out / "host.fa", [("chr1", host)])
    wfa(out / "refs/rrna/plant_rrna.fasta", [("RRNA1 Eukaryota;Archaeplastida;Chloroplastida;Testplant 18S", rrna)])
    genomes = [(f"NC_90000{i}.1 Test virus {i}, complete genome", mutate(rng, virus[f"V{i}"], 0.02)) for i in (1, 2, 3)]
    genomes += [(f"NC_90000{i}.1 Test virus {i}, complete genome", mutate(rng, virus[f"V{i}"], 0.12)) for i in (4, 5, 6)]
    genomes.append(("NC_900099.1 Unrelated virus Z, complete genome", rand_seq(rng, 5000)))
    wfa(out / "refs/viral/refseq_viral.fna", genomes)
    wfa(out / "refs/viral/refseq_viral_protein.faa",
        [("YP_900071.1 replicase [Test virus 7]", prot7[:half]),
         ("YP_900072.1 coat protein [Test virus 7]", prot7[half:]),
         ("YP_900099.1 capsid [Unrelated virus Z]", "".join(rng.choices(AA, k=300)))])
    wfa(out / "refs/protein/sprot_viridiplantae.faa",
        [("sp|P99998|OTHER_PLANT Other protein OS=Testplant OX=33090", "".join(rng.choices(AA, k=250)))])
    (out / "refs" / "REFS_VERSION.tsv").write_text(
        "name\tversion\tsource\nrefseq_viral_genomic\tRefSeq release 999\ttest fixture\n"
        "refseq_viral_protein\tRefSeq release 999\ttest fixture\n")
    seal_refs(out)

    base = {v: rng.randint(180, 380) for v in virus}
    meta = ["sample\tgroup"]
    truth = ["sample\tgroup\t" + "\t".join(virus)]
    q = "I" * 100
    for g in "AB":
        for k in range(1, 7):
            sample = f"{g}{k}"
            depth = rng.uniform(0.7, 1.4)
            parts = [("host", host, int(600 * depth)), ("bact", bact, int(3000 * depth)),
                     ("rrna", rrna, int(60 * depth))]
            n_virus = {}
            for v, seq in virus.items():
                n = base[v] * depth * rng.lognormvariate(0, 0.2) * (FOLD if (v == TRUE_DA and g == "B") else 1.0)
                n_virus[v] = int(round(n))
                parts.append((v, seq, n_virus[v]))
            reads = []
            for src, genome, n in parts:
                for i, (a, b) in enumerate(fragments(rng, genome, n)):
                    reads.append((f"{sample}_{src}_{i}", a, b, src))
            rng.shuffle(reads)
            with gzip.open(out / "reads" / f"{sample}_R1.fastq.gz", "wt") as f1, \
                    gzip.open(out / "reads" / f"{sample}_R2.fastq.gz", "wt") as f2, \
                    open(out / "reads" / f"{sample}.kraken2", "w") as fk:
                for rid, a, b, src in reads:
                    f1.write(f"@{rid}/1\n{a}\n+\n{q}\n")
                    f2.write(f"@{rid}/2\n{b}\n+\n{q}\n")
                    c = src == "bact"
                    fk.write(f"{'C' if c else 'U'}\t{rid}\t{2 if c else 0}\t100|100\t0:66 |:| 0:66\n")
            write_kreport(out / "reports" / f"{sample}.kreport", len(reads),
                          sum(1 for r in reads if r[3] == "bact"), sp_weights, rng_rep)
            meta.append(f"{sample}\t{g}")
            truth.append(f"{sample}\t{g}\t" + "\t".join(str(n_virus[v]) for v in virus))
    (out / "metadata.tsv").write_text("\n".join(meta) + "\n")
    (out / "truth_pairs.tsv").write_text("\n".join(truth) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 23)
