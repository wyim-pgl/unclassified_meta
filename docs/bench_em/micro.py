"""Micro-benchmark of one Gibbs sweep's pieces on a sample's largest crossing group: micro.py <classes.pkl>"""
import pickle
import sys
import timeit

import numpy as np

from unclassified_meta import em

c = pickle.load(open(sys.argv[1], "rb"))["classes"]
cross = [x for x in em._class_components(c) if any(len(k) > 1 for k in x) and not em._Tree(x).ok]
big = max(cross, key=len)
cl = em._Classes(big)
rng = np.random.default_rng(1)
M = 10
th = rng.standard_gamma(cl.single + 0.5, size=(M, cl.T))
print("K", cl.K, "T", cl.T, "width", cl.width, "reads", cl.n.sum(), "n==1 classes", int((cl.n == 1).sum()))
valid = cl.pad >= 0
rest = np.broadcast_to(np.rint(cl.n).astype(np.int64), (M, cl.K)).copy()
p = rng.random((M, cl.K))
w = np.where(valid, th[:, np.clip(cl.pad, 0, None)], 0.0)
idx = (np.arange(M)[:, None] * cl.T + np.clip(cl.pad[:, 0], 0, None)).ravel()
x = rng.binomial(rest, p)
tests = {
    "draw_assignments M=10": lambda: em._draw_assignments(cl, th, rng),
    "gather w": lambda: np.where(valid, th[:, np.clip(cl.pad, 0, None)], 0.0),
    "binomial (M,K)": lambda: rng.binomial(rest, p),
    "bincount (M*K)": lambda: np.bincount(idx, weights=x.ravel(), minlength=M * cl.T),
    "standard_gamma (M,T)": lambda: rng.standard_gamma(th + 0.5),
    "uniform (M,K)": lambda: rng.random((M, cl.K)),
}
for name, f in tests.items():
    n = 20
    print(f"{name}: {timeit.timeit(f, number=n) / n * 1e3:.2f} ms", flush=True)
