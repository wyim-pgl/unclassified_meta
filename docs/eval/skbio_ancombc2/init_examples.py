"""Print features whose initial coefficients differ most between R and scikit-bio."""
import sys

import numpy as np
import pandas as pd

d, counts_path, meta_path = sys.argv[1:4]
rd = lambda n: pd.read_csv(f"{d}/{n}", sep="\t", dtype={"feature_id": str}).set_index("feature_id")
rb, pb, est = rd("r_beta.tsv"), rd("py_beta.tsv"), rd("py_estimable.tsv")
rt, pt = np.loadtxt(f"{d}/r_theta.txt"), np.loadtxt(f"{d}/py_theta.txt")
print("initial theta max|R-Py|", np.abs(rt - pt).max())
c = pd.read_csv(counts_path, sep="\t", dtype={"feature_id": str}).set_index("feature_id")
m = pd.read_csv(meta_path, sep="\t", dtype=str).set_index("sample")
c = c.loc[rb.index, m.index]
db = (rb["conditionmacho"] - pb["conditionmacho"]).abs()
picks = db.sort_values(ascending=False).index[:2].tolist() + [db[(db > 1e-6) & (db < 0.05)].index[0]]
for f in picks:
    obs = c.loc[f].ne(0).groupby((m.condition + "_" + m.plant).values).sum()
    print("\n", f[:70])
    print("  observed samples per cell:", dict(obs))
    print("  R  :", rb.loc[f].round(4).to_dict())
    print("  Py :", pb.loc[f].round(4).to_dict())
    print("  Py estimable:", est.loc[f].to_dict())
