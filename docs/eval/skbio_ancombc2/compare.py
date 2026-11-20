"""Compare an R ancombc2_results.tsv with a scikit-bio result in the same layout.

usage: compare.py R.tsv PY.tsv [--json OUT.json] [--detail OUT.tsv]
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[3] / "src"))


def read(path):
    from unclassified_meta.da import read_results
    return read_results(path).drop_duplicates("feature_id").set_index("feature_id")


def tf(s):
    return s.astype(str).str.upper().eq("TRUE")


def untestable(df, t):
    """R's 'could not test': p = 1 with a finite nonzero W (dof NA)."""
    p, w = df[f"p_{t}"].astype(float), df[f"W_{t}"].astype(float)
    return (p >= 1) & np.isfinite(w) & (w != 0)


def compare(r, py):
    common = r.index.intersection(py.index)
    out = {"features_r": len(r), "features_py": len(py), "common": len(common),
           "only_r": len(r.index.difference(py.index)), "only_py": len(py.index.difference(r.index)), "terms": {}}
    r, py = r.loc[common], py.loc[common]
    terms = [c[4:] for c in r.columns if c.startswith("lfc_")]
    for t in terms:
        d = {}
        for s in ("lfc", "se", "W", "p", "q"):
            a, b = r[f"{s}_{t}"].astype(float).to_numpy(), py[f"{s}_{t}"].astype(float).to_numpy()
            both = np.isfinite(a) & np.isfinite(b)
            na_r, na_py = ~np.isfinite(a), ~np.isfinite(b)
            diff = np.abs(a[both] - b[both])
            rel = diff / np.maximum(np.abs(a[both]), 1e-300)
            d[s] = {"n_both_finite": int(both.sum()), "na_r_only": int((na_r & ~na_py).sum()),
                    "na_py_only": int((na_py & ~na_r).sum()), "na_both": int((na_r & na_py).sum()),
                    "max_abs_diff": float(diff.max()) if diff.size else None,
                    "median_abs_diff": float(np.median(diff)) if diff.size else None,
                    "max_rel_diff": float(rel.max()) if rel.size else None,
                    "n_rel_diff_gt_1e-6": int((rel > 1e-6).sum()),
                    "pearson": float(np.corrcoef(a[both], b[both])[0, 1]) if both.sum() > 2 else None}
        for call in ("diff", "passed_ss", "diff_robust"):
            col = f"{call}_{t}"
            if col in r.columns and col in py.columns:
                a, b = tf(r[col]), tf(py[col])
                d[call] = {"r_true": int(a.sum()), "py_true": int(b.sum()), "both": int((a & b).sum()),
                           "r_only": int((a & ~b).sum()), "py_only": int((~a & b).sum())}
        if f"W_{t}" in r.columns:
            ur, up = untestable(r, t), untestable(py, t)
            pr = py[f"p_{t}"].astype(float)
            d["untestable_r"] = {"n": int(ur.sum()),
                                 "py_p_is_nan": int((ur & pr.isna()).sum()),
                                 "py_p_is_1": int((ur & (pr >= 1)).sum()),
                                 "py_p_lt_1": int((ur & (pr < 1)).sum()),
                                 "py_p_lt_0.05": int((ur & (pr < 0.05)).sum())}
            d["untestable_py"] = int(up.sum())
            d["py_p_nan_total"] = int(pr.isna().sum())
        out["terms"][t] = d
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("r"); ap.add_argument("py"); ap.add_argument("--json")
    a = ap.parse_args()
    res = compare(read(a.r), read(a.py))
    txt = json.dumps(res, indent=1)
    if a.json:
        open(a.json, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
