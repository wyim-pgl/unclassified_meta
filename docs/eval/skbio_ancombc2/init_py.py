"""scikit-bio initial ANCOM-BC2 estimates (the _ancombc_core steps before sampling-fraction correction).
usage: init_py.py counts.tsv meta.tsv out_prefix"""
import sys

import numpy as np
import pandas as pd
from skbio.stats.composition import _ancombc as A

counts_path, meta_path, out = sys.argv[1:4]
sys.path.insert(0, sys.path[0])
from skbio_fit import load, r_name  # noqa: E402

counts, meta = load(counts_path, meta_path, "condition", "healthy", "condition + plant")
counts = counts.loc[(counts != 0).mean(axis=1) >= 0.10]
matrix = counts.T.to_numpy(dtype=float)
data, missing = A._transform_data(matrix, 0, center=True)
dmat = A._build_dmatrix("condition + plant", meta, dtype=data.dtype)
covars = [r_name(c) for c in dmat.design_info.column_names]
var_hat, beta, theta, _, estimable, _ = A._estimate_params(data, dmat, None, missing, True, True, True)
pd.DataFrame(beta.T, index=counts.index, columns=covars).rename_axis("feature_id").to_csv(f"{out}_beta.tsv", sep="\t")
pd.DataFrame(var_hat, index=counts.index, columns=covars).rename_axis("feature_id").to_csv(f"{out}_var.tsv", sep="\t")
if estimable is not None:
    pd.DataFrame(estimable, index=counts.index, columns=covars).rename_axis("feature_id").to_csv(f"{out}_estimable.tsv", sep="\t")
rows = []
for i, c in enumerate(covars):
    if estimable is None:
        b, v = beta[i], var_hat[:, i]
    else:
        ok = estimable[:, i] & np.isfinite(var_hat[:, i])
        b, v = beta[i, ok], var_hat[ok, i]
    rows.append([c, *A._estimate_bias_em(b, v, tol=1e-5, max_iter=100)])
pd.DataFrame(rows, columns=["term", "delta_em", "delta_wls", "var_delta"]).to_csv(f"{out}_bias.tsv", sep="\t", index=False)
np.savetxt(f"{out}_theta.txt", theta)
print(pd.DataFrame(rows, columns=["term", "delta_em", "delta_wls", "var_delta"]).to_string(index=False))
