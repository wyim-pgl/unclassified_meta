"""Synthetic repeated-measures count table for the `da --rand-formula` test (no reads, no downloads).

Design: 2 groups (H, M) x 6 units per group x 2 replicate samples per unit = 24 samples, 18 features.
Every feature has a unit effect shared by the replicates of a unit (the replicates are not independent);
two features are truly 8 times more abundant in group M. Sample depth varies.

Writes counts.tsv, features.tsv, metadata.tsv and truth.tsv in the layout `unclassified-meta da` reads.
"""
from __future__ import annotations

import math
import random
import sys
from pathlib import Path

TRUE_DA = ("taxon_003", "taxon_017")
FOLD = 8.0


def poisson(rng: random.Random, lam: float) -> int:
    if lam > 50:
        return max(0, int(round(rng.gauss(lam, math.sqrt(lam)))))
    limit, k, p = math.exp(-lam), 0, 1.0
    while True:
        p *= rng.random()
        if p <= limit:
            return k
        k += 1


def main(outdir: str, seed: int = 5, unit_sd: float = 0.3, n_features: int = 18, units: int = 6, reps: int = 2) -> None:
    rng = random.Random(seed)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    feats = [f"taxon_{i:03d}" for i in range(1, n_features + 1)]
    base = {f: rng.lognormvariate(5.5, 1.0) for f in feats}
    samples, meta = [], ["sample\tgroup\tunit\trep"]
    counts = {f: [] for f in feats}
    for g in "HM":
        for u in range(1, units + 1):
            unit = f"{g}_plant{u}"
            unit_eff = {f: rng.gauss(0, unit_sd) for f in feats}
            for r in range(1, reps + 1):
                s = f"{unit}_rep{r}"
                depth = rng.uniform(0.6, 1.6)
                samples.append(s)
                meta.append(f"{s}\t{g}\t{unit}\t{r}")
                for f in feats:
                    lam = base[f] * depth * math.exp(unit_eff[f] + rng.gauss(0, 0.15))
                    if f in TRUE_DA and g == "M":
                        lam *= FOLD
                    counts[f].append(poisson(rng, lam))
    (out / "metadata.tsv").write_text("\n".join(meta) + "\n")
    with open(out / "counts.tsv", "w") as fh:
        fh.write("feature_id\t" + "\t".join(samples) + "\n")
        for f in feats:
            fh.write(f + "\t" + "\t".join(map(str, counts[f])) + "\n")
    with open(out / "features.tsv", "w") as fh:
        fh.write("feature_id\tlabel\tnamespace\taccession\tidentity_class\n")
        for f in feats:
            fh.write(f"{f}\t{f.replace('_', ' ')}\tsynthetic\t{f}\ttaxon\n")
    (out / "truth.tsv").write_text("feature_id\tfold_M_vs_H\n" + "".join(f"{f}\t{FOLD}\n" for f in TRUE_DA))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 5)
