"""Break the R vs Python differences down by the number of zero samples of each feature.

usage: strata.py R.tsv PY.tsv COUNTS.tsv META.tsv TERM
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[3] / "src"))
from unclassified_meta.da import read_results  # noqa: E402

r = read_results(sys.argv[1]).set_index("feature_id")
py = read_results(sys.argv[2]).set_index("feature_id")
counts = pd.read_csv(sys.argv[3], sep="\t", dtype={"feature_id": str}).set_index("feature_id")
meta = pd.read_csv(sys.argv[4], sep="\t", dtype=str).set_index("sample")
t = sys.argv[5]
counts = counts.loc[r.index, meta.index]
nz = (counts == 0).sum(axis=1)
# a plant x condition cell without any observation makes the feature-specific design rank deficient
cells = meta["condition"] + "_" + meta["plant"]
empty_cells = counts.ne(0).T.groupby(cells.values).sum().T.eq(0).sum(axis=1)
df = pd.DataFrame({"zeros": nz, "empty_cells": empty_cells,
                   "lfc_r": r[f"lfc_{t}"], "lfc_py": py[f"lfc_{t}"], "se_r": r[f"se_{t}"], "se_py": py[f"se_{t}"],
                   "p_r": r[f"p_{t}"], "p_py": py[f"p_{t}"]})
df["rel_lfc"] = (df.lfc_r - df.lfc_py).abs() / df.lfc_r.abs().clip(lower=1e-12)
df["rel_se"] = (df.se_r - df.se_py).abs() / df.se_r.abs()
df["dp"] = (df.p_r - df.p_py).abs()
bins = pd.cut(df.zeros, [-1, 0, 1, 2, 4, 8, 12, 18], labels=["0", "1", "2", "3-4", "5-8", "9-12", "13-18"])
g = df.groupby(bins, observed=True)
print(f"term {t}")
print(pd.DataFrame({"n": g.size(), "lfc_na_both": g.apply(lambda x: (x.lfc_r.isna() & x.lfc_py.isna()).sum()),
                    "p_py_nan": g.apply(lambda x: x.p_py.isna().sum()),
                    "max_rel_lfc": g.rel_lfc.max(), "med_rel_lfc": g.rel_lfc.median(),
                    "max_rel_se": g.rel_se.max(), "max_dp": g.dp.max(), "med_dp": g.dp.median()}).to_string())
print("by empty condition x plant cells:")
g2 = df.groupby("empty_cells")
print(pd.DataFrame({"n": g2.size(), "p_py_nan": g2.apply(lambda x: x.p_py.isna().sum()),
                    "lfc_na_r": g2.apply(lambda x: x.lfc_r.isna().sum()),
                    "max_rel_lfc": g2.rel_lfc.max(), "max_dp": g2.dp.max()}).to_string())
