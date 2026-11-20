"""Many imputations of one sample with the code on PYTHONPATH: imp_draws.py <classes.pkl> <out.npz> <M> <seed>
(em.imputations does not use the EM fit, so only the classes are needed.)"""
import pickle
import sys
import time

import numpy as np

from unclassified_meta import em

classes = pickle.load(open(sys.argv[1], "rb"))["classes"]
m, seed = int(sys.argv[3]), int(sys.argv[4])
t = time.time()
draws = em.imputations(classes, em.Fit({}, {}, 0, True), m, np.random.default_rng(seed))
dt = time.time() - t
targets = sorted({k for d in draws for k in d})
mat = np.array([[d.get(k, 0) for k in targets] for d in draws], dtype=np.int64)
np.savez_compressed(sys.argv[2], targets=np.array(targets), draws=mat, seconds=dt,
                    reads=sum(classes.values()))
print(sys.argv[1].rsplit("/", 1)[1], f"M {m} {dt:.1f}s", flush=True)
