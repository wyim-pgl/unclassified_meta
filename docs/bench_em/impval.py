"""Imputation distribution: baseline vs optimised, with optimised vs optimised (other seed) as the null reference.
impval.py <dir with npz> <sample stem>...  Each npz: M independent imputations of one sample.

Per target (mean >= 1 read in either set): z = mean difference / its standard error, and the ratio of SDs.
Under equal distributions z ~ N(0, 1) (M = 200 independent chains per set), so mean z^2 ~ 1 and |z| > 3 in
~0.27 % of targets; the opt-vs-opt row shows the same statistics for two runs of one sampler."""
import sys

import numpy as np


def load(path):
    d = np.load(path)
    draws, targets = d["draws"], d["targets"]          # read each array from the archive once
    return {t: draws[:, i] for i, t in enumerate(targets)}, int(d["reads"]), float(d["seconds"]), draws.shape[0]


def compare(a, b):
    keys = sorted(set(a) | set(b))
    m = len(next(iter(a.values())))
    A = np.array([a.get(k, np.zeros(m)) for k in keys], dtype=float)
    B = np.array([b.get(k, np.zeros(m)) for k in keys], dtype=float)
    ma, mb, va, vb = A.mean(1), B.mean(1), A.var(1, ddof=1), B.var(1, ddof=1)
    keep = (np.maximum(ma, mb) >= 1) & ((va + vb) > 0)
    z = (ma - mb)[keep] / np.sqrt((va + vb)[keep] / m)
    both = keep & (va > 0) & (vb > 0)
    r = np.sqrt(va[both] / vb[both])
    lr = np.log(r)
    return {"targets": int(keep.sum()), "mean_z2": float(np.mean(z ** 2)), "frac_|z|>3": float(np.mean(np.abs(z) > 3)),
            "max_|z|": float(np.abs(z).max()), "sd_ratio_median": float(np.median(r)),
            "sd_ratio_q01_q99": (float(np.quantile(r, 0.01)), float(np.quantile(r, 0.99))),
            "mean_log_sd_ratio": float(lr.mean()), "se_mean_log_sd_ratio": float(lr.std() / np.sqrt(len(lr)))}


d = sys.argv[1]
for stem in sys.argv[2:]:
    base, nb, tb, m = load(f"{d}/baseline_src_{stem}_M200_s101.npz")
    o1, n1, t1, _ = load(f"{d}/opt3_src_{stem}_M200_s202.npz")
    o2, n2, t2, _ = load(f"{d}/opt3_src_{stem}_M200_s303.npz")
    mass = all((np.sum([v for v in x.values()], axis=0) == n).all() for x, n in ((base, nb), (o1, n1), (o2, n2)))
    print(f"## {stem}: M {m}, reads {nb}, every draw sums to the reads: {mass}; time base {tb:.0f}s opt {t1:.0f}s "
          f"/ {t2:.0f}s (x{tb / t1:.1f})")
    for name, (x, y) in {"base vs opt": (base, o1), "opt vs opt (null)": (o2, o1)}.items():
        print(f"  {name}: " + ", ".join(f"{k} {v:.3g}" if isinstance(v, float) else
                                        f"{k} {v}" if not isinstance(v, tuple) else
                                        f"{k} {v[0]:.3g}-{v[1]:.3g}" for k, v in compare(x, y).items()))
