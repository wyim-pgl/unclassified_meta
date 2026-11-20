"""Per-sample structure of the EM input: reads, targets, ambiguous reads, classes, crossing components."""
import glob
import pickle
import sys

from unclassified_meta import em

prefix = sys.argv[1]           # e.g. data/RNA_S_
print("sample\treads\ttargets\tambiguous_reads\tclasses\tambiguous_classes\tmean_candidates_ambiguous\t"
      "components\tcrossing\tlargest_crossing_targets\tlargest_crossing_classes\tparse_s")
for f in sorted(glob.glob(prefix + "*.pkl"), key=lambda p: int(p.split("_")[2])):
    d = pickle.load(open(f, "rb"))
    c = d["classes"]
    amb = {k: n for k, n in c.items() if len(k) > 1}
    na = sum(amb.values())
    comps = em._class_components(c)
    cross = [x for x in comps if any(len(k) > 1 for k in x) and not em._Tree(x).ok]
    big = max(cross, key=len) if cross else {}
    print(d["sample"], sum(c.values()), len({t for k in c for t in k}), na, len(c), len(amb),
          round(sum(len(k) * n for k, n in amb.items()) / na, 3) if na else 0, len(comps), len(cross),
          len({t for k in big for t in k}), len(big), round(d["parse_s"], 1), sep="\t")
