"""Expected counts, optimised vs baseline, for every sample: cmp_counts.py <base dir> <new dir> <classes dir>
Differences in reads, relative to the count and to the baseline imputation SD of each target (M draws)."""
import glob
import pickle
import sys

import numpy as np

from unclassified_meta import em

print("sample\ttargets\treads\tsum_abs_diff\tmax_abs_diff_taxa\tmax_abs_diff_tokens\tmax_rel_diff_taxa_ge10reads\t"
      "max_diff_over_impSD\ttargets_diff_gt_0.1SD\tll_new_minus_base\tmass_new_minus_reads")
for f in sorted(glob.glob(sys.argv[1] + "/*.pkl"), key=lambda p: int(p.rsplit("/", 1)[1].split("_")[2])):
    name = f.rsplit("/", 1)[1]
    a, b = pickle.load(open(f, "rb")), pickle.load(open(f"{sys.argv[2]}/{name}", "rb"))
    classes = pickle.load(open(f"{sys.argv[3]}/{name}", "rb"))["classes"]
    keys = sorted(set(a["counts"]) | set(b["counts"]))
    x = np.array([a["counts"].get(k, 0.0) for k in keys])
    y = np.array([b["counts"].get(k, 0.0) for k in keys])
    d = np.abs(y - x)
    taxa = np.array([not k.startswith("#") for k in keys])
    imp = np.array([[dd.get(k, 0) for k in keys] for dd in a["imputations"]], dtype=float)
    sd = imp.std(axis=0, ddof=1)
    ratio = np.where(sd > 0, d / np.where(sd > 0, sd, 1), np.where(d > 1e-6, np.inf, 0))
    big = taxa & (x >= 10)
    c = em._Classes(classes)
    ll = lambda cnt: c.loglik(np.array([cnt.get(t, 0.0) for t in c.targets]) / c.total)
    print(a["sample"], len(keys), int(c.total), f"{d.sum():.3g}", f"{d[taxa].max():.3g}", f"{d[~taxa].max():.3g}",
          f"{(d[big] / x[big]).max():.2e}", f"{ratio.max():.3g}", int((ratio > 0.1).sum()),
          f"{ll(b['counts']) - ll(a['counts']):+.2e}", f"{y.sum() - c.total:+.2e}", sep="\t")
