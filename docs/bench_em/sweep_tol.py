"""Dual-certificate tolerance sweep against the baseline: sweep_tol.py <classes.pkl> <baseline result.pkl>"""
import pickle
import sys
import time

import numpy as np

from unclassified_meta import em

classes = pickle.load(open(sys.argv[1], "rb"))["classes"]
base = pickle.load(open(sys.argv[2], "rb"))
c = em._Classes(classes)


def ll(counts):
    th = np.array([counts.get(t, 0.0) for t in c.targets]) / c.total
    return c.loglik(th)


lb = ll(base["counts"])
for tol in [0.1, 0.01, 1e-3, 1e-4]:
    t = time.time()
    fit = em.expected_counts(classes, dual_tol=tol)
    dt = time.time() - t
    keys = sorted(set(fit.counts) | set(base["counts"]))
    d = np.array([fit.counts.get(k, 0.0) - base["counts"].get(k, 0.0) for k in keys])
    kept = np.array([not k.startswith("#") for k in keys])
    print(f"{base['sample']} tol {tol:g}: {dt:.1f}s it {fit.iterations} conv {fit.converged} (base {base['t_em']:.0f}s "
          f"it {base['iterations']}) ll new-base {ll(fit.counts) - lb:+.5f} max|diff| all {np.abs(d).max():.3g} "
          f"taxa {np.abs(d[kept]).max():.3g} sum|diff| {np.abs(d).sum():.3g}", flush=True)
