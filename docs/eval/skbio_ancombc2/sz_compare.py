"""Compare scikit-bio struc_zero with R's .get_struc_zero output.  usage: sz_compare.py counts meta r_prefix"""
import sys

import pandas as pd
from skbio.stats.composition import struc_zero

counts_path, meta_path, rp = sys.argv[1:4]
counts = pd.read_csv(counts_path, sep="\t", dtype={"feature_id": str}).set_index("feature_id")
meta = pd.read_csv(meta_path, sep="\t", dtype=str).set_index("sample")
table = counts[meta.index].T
for nl, suffix in ((False, ""), (True, "_neglb")):
    py = struc_zero(table, meta, "condition", neg_lb=nl)
    r = pd.read_csv(f"{rp}{suffix}.tsv", sep="\t", dtype={"feature_id": str}).set_index("feature_id")
    r = r.astype(str).apply(lambda s: s.str.upper().eq("TRUE")).loc[py.index]
    for g in ("healthy", "macho"):
        a, b = r[g], py[g].astype(bool)
        print(f"neg_lb={nl} {g:8s} features {len(a)}  R {int(a.sum())}  Py {int(b.sum())}  "
              f"agree {int((a == b).sum())}  disagree {int((a != b).sum())}")
