import sys, glob, os
sys.path.insert(0, "/mnt/f/scratch/develop/unclassified_meta-p257fix/src")
import numpy as np, pandas as pd
from unclassified_meta import da
root = "/home/wyim/.claude/jobs/d4415b8c/tmp/sim/runs"
n_runs = n_diff_call = n_diff_p = 0; maxdq = 0.0; fails = 0
for d in sorted(glob.glob(f"{root}/*/da_em")):
    fits = sorted(glob.glob(f"{d}/imputations/ancombc2_results_*.tsv"))
    if not fits or not os.path.exists(f"{d}/mi_results.tsv"): continue
    old = da.read_results(f"{d}/mi_results.tsv")
    new, terms = da.combine_imputations([da.read_results(f) for f in fits], 0.05, "BH")
    t = terms[0]; o = old.set_index("feature_id"); n = new.set_index("feature_id")
    n_runs += 1
    n_diff_call += int((o[f"diff_{t}"].astype(str).str.upper().eq("TRUE") != n.loc[o.index, f"diff_{t}"].astype(bool)).sum())
    dp = (o[f"p_{t}"] - n.loc[o.index, f"p_{t}"]).abs()
    n_diff_p += int((dp > 1e-9).sum()); maxdq = max(maxdq, float((o[f"q_{t}"] - n.loc[o.index, f"q_{t}"]).abs().max()))
    for f in fits:
        r = da.read_results(f); fails += int(((r[f"p_{t}"] >= 1) & (r[f"W_{t}"].abs() > 0)).sum()) if f"W_{t}" in r else 0
print(f"runs {n_runs}; features whose call changed {n_diff_call}; p changed {n_diff_p}; max |dq| {maxdq:.3g}; per-fit p=1 with W!=0: {fails}")
