"""Run em.expected_counts and em.imputations on one dumped sample with stage timers; pickle the results.

usage: run_em.py <classes.pkl> <out.pkl> <M> <seed> <sample_index>
The code under test is whatever `unclassified_meta` PYTHONPATH points at (baseline or optimised)."""
import collections
import functools
import pickle
import resource
import sys
import time

import numpy as np

from unclassified_meta import em

src, out, M, seed, si = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
T = collections.defaultdict(float)
N = collections.Counter()
ITERS = []


def timed(name, owner, attr):
    f = getattr(owner, attr)

    @functools.wraps(f)
    def g(*a, **k):
        t = time.perf_counter()
        try:
            r = f(*a, **k)
        finally:
            T[name] += time.perf_counter() - t
            N[name] += 1
        if name == "squarem":
            ITERS.append(r[1])
        return r
    setattr(owner, attr, g)


for name in ["_class_components", "dominated", "_squarem", "_chain", "_draw_assignments", "_tree_theta",
             "_laminar_part", "_gap", "_sir_starts", "_chains"]:
    if hasattr(em, name):
        timed(name.strip("_"), em, name)
timed("estep", em._Classes, "estep")
timed("loglik", em._Classes, "loglik")
timed("classes_init", em._Classes, "__init__")
timed("tree_flow", em._Tree, "flow")
timed("tree_init", em._Tree, "__init__")

d = pickle.load(open(src, "rb"))
classes = d["classes"]
t0 = time.perf_counter()
fit = em.expected_counts(classes)
t_em = time.perf_counter() - t0
snap = dict(T)
t0 = time.perf_counter()
imps = em.imputations(classes, em.Fit({}, fit.theta, 0, True), M, np.random.default_rng([seed, 1, si])) if M else []
t_imp = time.perf_counter() - t0
imp_stage = {k: T[k] - snap.get(k, 0.0) for k in T}
res = {"sample": d["sample"], "counts": fit.counts, "theta": fit.theta, "iterations": fit.iterations,
       "converged": fit.converged, "imputations": imps, "t_em": t_em, "t_imp": t_imp, "parse_s": d["parse_s"],
       "em_stage": snap, "imp_stage": imp_stage, "calls": dict(N), "squarem_iters": ITERS,
       "maxrss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
pickle.dump(res, open(out, "wb"))
print(d["sample"], f"em {t_em:.1f}s it {fit.iterations} conv {fit.converged} imp{M} {t_imp:.1f}s "
      f"rss {res['maxrss_kb'] / 1024:.0f}MB", flush=True)
