"""The targets with the largest |z| between baseline and optimised imputations, with both optimised seeds and a
permutation (rank) check: impval_top.py <dir> <stem>"""
import sys

import numpy as np

d, stem = sys.argv[1], sys.argv[2]


def load(tag):
    z = np.load(f"{d}/{tag}_{stem}_M200_{'s101' if tag == 'baseline_src' else 's202'}.npz") if tag != "opt2" else \
        np.load(f"{d}/opt3_src_{stem}_M200_s303.npz")
    dr, tg = z["draws"], z["targets"]
    return {t: dr[:, i] for i, t in enumerate(tg)}


base, o1, o2 = load("baseline_src"), load("opt3_src"), load("opt2")
keys = sorted(set(base) | set(o1) | set(o2))
zero = np.zeros(200)
rows = []
for k in keys:
    a, b, c = base.get(k, zero), o1.get(k, zero), o2.get(k, zero)
    if max(a.mean(), b.mean()) < 1 or a.var() + b.var() == 0:
        continue
    z = (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 200)
    rows.append((abs(z), z, k, a, b, c))
rows.sort(key=lambda r: -r[0])
for az, z, k, a, b, c in rows[:6]:
    pooled = np.concatenate([a, b, c])
    # Mann-Whitney-type rank comparison of base against both optimised runs together
    ranks = pooled.argsort().argsort()
    u = ranks[:200].mean() / len(pooled)
    print(f"{k}: z {z:+.2f} | base mean {a.mean():.1f} sd {a.std():.1f} | opt {b.mean():.1f} sd {b.std():.1f} | "
          f"opt2 {c.mean():.1f} sd {c.std():.1f} | base rank share {u:.3f} (0.5 if equal) | "
          f"zeros base/opt/opt2 {(a == 0).sum()}/{(b == 0).sum()}/{(c == 0).sum()}")
