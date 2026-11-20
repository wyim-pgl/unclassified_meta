"""One summary row per fit: R vs scikit-bio on the tested term, with Holm (as fitted) and BH (re-adjusted from the
same p-values, p = 1 for untestable features as R does).

usage: summary.py TERM NAME=R.tsv,PY.tsv [NAME=R.tsv,PY.tsv ...]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from unclassified_meta.da import p_adjust, read_results  # noqa: E402


def tf(s):
    return s.astype(str).str.upper().eq("TRUE")


def calls(a, b):
    return f"{int(a.sum())}/{int(b.sum())}/{int((a & b).sum())}"


t = sys.argv[1]
rows = []
for arg in sys.argv[2:]:
    name, paths = arg.split("=", 1)
    rp, pp = paths.split(",")
    r = read_results(rp).set_index("feature_id")
    p = read_results(pp).set_index("feature_id").reindex(r.index)
    lr, lp = r[f"lfc_{t}"].astype(float), p[f"lfc_{t}"].astype(float)
    both = lr.notna() & lp.notna()
    d = (lr - lp)[both]
    se_rel = ((r[f"se_{t}"] - p[f"se_{t}"]).abs() / r[f"se_{t}"])[both]
    pr, pq = r[f"p_{t}"].astype(float), p[f"p_{t}"].astype(float).fillna(1.0)
    unt_r = (pr >= 1) & r[f"W_{t}"].astype(float).pipe(lambda w: np.isfinite(w) & (w != 0))
    unt_p = (pq >= 1) & p[f"W_{t}"].astype(float).pipe(lambda w: np.isfinite(w) & (w != 0))
    bh_r, bh_p = p_adjust(pr.to_numpy(), "BH") < 0.05, p_adjust(pq.to_numpy(), "BH") < 0.05
    rows.append({"fit": name, "features": len(r), "lfc_est_both": int(both.sum()),
                 "lfc_pearson": round(float(np.corrcoef(lr[both], lp[both])[0, 1]), 4),
                 "median|dlfc|": round(float(d.abs().median()), 4), "max|dlfc|": round(float(d.abs().max()), 3),
                 "max_rel_dse": round(float(se_rel.max()), 4),
                 "untestable R/Py/both": calls(unt_r, unt_p),
                 "holm diff R/Py/both": calls(tf(r[f"diff_{t}"]), tf(p[f"diff_{t}"])),
                 "holm robust R/Py/both": calls(tf(r[f"diff_robust_{t}"]), tf(p[f"diff_robust_{t}"])),
                 "BH q<0.05 R/Py/both": calls(pd.Series(bh_r), pd.Series(bh_p))})
print(pd.DataFrame(rows).to_string(index=False))
