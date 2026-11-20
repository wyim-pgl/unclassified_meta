"""Write feature subsets of a count table after the prv_cut filter: zero-free features, and features whose every
condition x plant cell has an observation (full-rank feature-specific design)."""
import sys

import pandas as pd

counts_path, meta_path, outdir = sys.argv[1:4]
counts = pd.read_csv(counts_path, sep="\t", dtype={"feature_id": str}).set_index("feature_id")
meta = pd.read_csv(meta_path, sep="\t", dtype=str).set_index("sample")
counts = counts[meta.index]
counts = counts.loc[(counts != 0).mean(axis=1) >= 0.10]
cells = (meta["condition"] + "_" + meta["plant"]).values
empty = counts.ne(0).T.groupby(cells).sum().T.eq(0).sum(axis=1)
counts.loc[(counts != 0).all(axis=1)].to_csv(f"{outdir}/sub_zerofree.tsv", sep="\t")
counts.loc[empty == 0].to_csv(f"{outdir}/sub_fullrank.tsv", sep="\t")
counts.loc[empty <= 2].to_csv(f"{outdir}/sub_le2empty.tsv", sep="\t")
print((counts != 0).all(axis=1).sum(), (empty == 0).sum(), (empty <= 2).sum())
