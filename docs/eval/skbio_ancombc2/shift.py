"""Is the R - Python LFC difference a constant shift per term (a different bias estimate delta_em)?"""
import sys

import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[3] / "src"))
from unclassified_meta.da import read_results  # noqa: E402

r = read_results(sys.argv[1]).set_index("feature_id")
py = read_results(sys.argv[2]).set_index("feature_id")
counts = pd.read_csv(sys.argv[3], sep="\t", dtype={"feature_id": str}).set_index("feature_id").loc[r.index]
full = counts.index[(counts != 0).all(axis=1)]
for t in [c[4:] for c in r.columns if c.startswith("lfc_")]:
    d = (r.loc[full, f"lfc_{t}"] - py.loc[full, f"lfc_{t}"])
    ratio = (r.loc[full, f"se_{t}"] / py.loc[full, f"se_{t}"])
    print(f"{t:16s} zero-free n={len(full)}  lfc R-Py: mean {d.mean():+.6f} sd {d.std():.2e} min {d.min():+.6f} "
          f"max {d.max():+.6f} | se R/Py: mean {ratio.mean():.5f} sd {ratio.std():.2e}")
