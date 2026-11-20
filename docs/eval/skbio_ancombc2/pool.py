"""Pool per-imputation fits with the worktree's da.combine_imputations (commit 537f358 rule included) and compare
the pooled tables of two backends.

usage: pool.py OUTDIR TERM --a R1.tsv R2.tsv ... --b P1.tsv P2.tsv ... [--p-adjust holm]
"""
import argparse
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[3] / "src"))
from unclassified_meta import da  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("outdir"); ap.add_argument("term")
ap.add_argument("--a", nargs="+", required=True); ap.add_argument("--b", nargs="+", required=True)
ap.add_argument("--p-adjust", default="holm"); ap.add_argument("--alpha", type=float, default=0.05)
a = ap.parse_args()
t = a.term
pa, _ = da.combine_imputations([da.read_results(f) for f in a.a], a.alpha, a.p_adjust)
pb, _ = da.combine_imputations([da.read_results(f) for f in a.b], a.alpha, a.p_adjust)
pa.to_csv(f"{a.outdir}/pooled_a.tsv", sep="\t", index=False)
pb.to_csv(f"{a.outdir}/pooled_b.tsv", sep="\t", index=False)
pa, pb = pa.set_index("feature_id"), pb.set_index("feature_id")
common = pa.index.intersection(pb.index)
print(f"features a {len(pa)} b {len(pb)} common {len(common)}")
pa, pb = pa.loc[common], pb.loc[common]
na, nb = pa[f"n_fits_{t}"], pb[f"n_fits_{t}"]
M = len(a.a)
print(f"n_fits=={M}: a {int((na == M).sum())} b {int((nb == M).sum())}; n_fits equal {int((na == nb).sum())}/{len(common)}")
print(pd.crosstab(na.rename("n_fits_a"), nb.rename("n_fits_b")).to_string())
pooled = (na == M) & (nb == M)
for s in ("lfc", "se", "df", "fmi", "p", "q"):
    x, y = pa.loc[pooled, f"{s}_{t}"].astype(float), pb.loc[pooled, f"{s}_{t}"].astype(float)
    fin = np.isfinite(x) & np.isfinite(y)
    d = (x - y).abs()[fin]
    print(f"{s:4s} pooled in both n={int(pooled.sum())} finite {int(fin.sum())}  max|diff| {d.max():.3g} "
          f"median {d.median():.3g}  pearson {np.corrcoef(x[fin], y[fin])[0, 1]:.6f}  inf a {int(np.isinf(x).sum())} b {int(np.isinf(y).sum())}")
for c in ("diff", "diff_robust"):
    x, y = pa[f"{c}_{t}"].astype(bool), pb[f"{c}_{t}"].astype(bool)
    print(f"{c}: a {int(x.sum())} b {int(y.sum())} both {int((x & y).sum())} a_only {list(pa.index[x & ~y].map(lambda i: pa.loc[i, 'label']))} "
          f"b_only {list(pb.index[~x & y].map(lambda i: pa.loc[i, 'label']))}")
sig = pa.index[(pa[f"q_{t}"] < 0.1) | (pb[f"q_{t}"] < 0.1)]
print(pd.DataFrame({"label": pa.loc[sig, "label"].str[:30], "lfc_a": pa.loc[sig, f"lfc_{t}"], "lfc_b": pb.loc[sig, f"lfc_{t}"],
                    "fmi_a": pa.loc[sig, f"fmi_{t}"], "fmi_b": pb.loc[sig, f"fmi_{t}"], "df_a": pa.loc[sig, f"df_{t}"],
                    "df_b": pb.loc[sig, f"df_{t}"], "q_a": pa.loc[sig, f"q_{t}"], "q_b": pb.loc[sig, f"q_{t}"],
                    "n_a": pa.loc[sig, f"n_fits_{t}"], "n_b": pb.loc[sig, f"n_fits_{t}"]}).round(4).to_string(index=False))
