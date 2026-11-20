"""Compare R and scikit-bio initial (pre-bias-correction) coefficients per feature, grouped by design deficiency."""
import sys

import numpy as np
import pandas as pd

d, counts_path, meta_path = sys.argv[1:4]
rb = pd.read_csv(f"{d}/r_beta.tsv", sep="\t", dtype={"feature_id": str}).set_index("feature_id")
pb = pd.read_csv(f"{d}/py_beta.tsv", sep="\t", dtype={"feature_id": str}).set_index("feature_id")
rv = pd.read_csv(f"{d}/r_var.tsv", sep="\t", dtype={"feature_id": str}).set_index("feature_id")
pv = pd.read_csv(f"{d}/py_var.tsv", sep="\t", dtype={"feature_id": str}).set_index("feature_id")
est = pd.read_csv(f"{d}/py_estimable.tsv", sep="\t", dtype={"feature_id": str}).set_index("feature_id")
counts = pd.read_csv(counts_path, sep="\t", dtype={"feature_id": str}).set_index("feature_id")
meta = pd.read_csv(meta_path, sep="\t", dtype=str).set_index("sample")
counts = counts.loc[rb.index, meta.index]
obs = counts.ne(0)
# what is missing for each feature
plant_absent = obs.T.groupby(meta["plant"].values).sum().T.eq(0)         # a whole plant unobserved
cond_absent = obs.T.groupby(meta["condition"].values).sum().T.eq(0)      # a whole condition unobserved
cells = (meta["condition"] + "_" + meta["plant"]).values
empty = obs.T.groupby(cells).sum().T.eq(0).sum(axis=1)
cls = np.where(cond_absent.any(axis=1), "condition level absent",
       np.where(plant_absent["P257"], "reference plant P257 absent",
       np.where(plant_absent.any(axis=1), "other plant level absent",
       np.where(empty > 0, "cells empty, levels present", "full rank"))))
t = "conditionmacho"
df = pd.DataFrame({"cls": cls, "b_r": rb[t], "b_py": pb[t], "v_r": rv[t], "v_py": pv[t], "est_py": est[t]})
df["used_r"] = df.b_r.notna() & df.v_r.notna()
df["used_py"] = df.est_py & np.isfinite(df.v_py)
df["db"] = (df.b_r - df.b_py).abs()
g = df.groupby("cls")
print(f"initial coefficient {t}: R vs scikit-bio, by what the feature's design lacks")
print(pd.DataFrame({"n": g.size(), "in_R_EM": g.used_r.sum(), "in_Py_EM": g.used_py.sum(),
                    "max|db|_both": g.apply(lambda x: x.db[x.used_r & x.used_py].max()),
                    "n|db|>1e-6": g.apply(lambda x: int((x.db[x.used_r & x.used_py] > 1e-6).sum()))}).to_string())
ex = df[(df.cls == "reference plant P257 absent") & df.used_r].head(3)
print(ex.round(4).to_string())
