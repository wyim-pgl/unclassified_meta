"""Run the real taxa-table --assign em parse for every sample and pickle the equivalence classes it hands to em.
em.expected_counts is replaced by a recorder (even split) so only the parse runs."""
import csv
import pickle
import sys
import time
from pathlib import Path

from unclassified_meta import em, taxa

T, R, OUT = sys.argv[1], sys.argv[2], Path(sys.argv[3])
W = "/data/gpfs/assoc/pgl/develop/unclassified_meta_test/OSD18"
TX = "/data/gpfs/assoc/pgl/develop/unclassified_meta_test/osd_taxonomy"
db = "cfr_nr" if T == "RNA" else "cfr_core_nt"
thr = dict(min_score=64) if T == "RNA" else dict(min_hitlen=25)
reports = {r["sample"]: Path(r["report"]) for r in csv.DictReader(open(f"{W}/reports_{T}.tsv"), delimiter="\t")}
OUT.mkdir(parents=True, exist_ok=True)
order, done, t_last = list(reports), [0], [time.time()]


def rec(classes, *a, **k):
    s = order[done[0]]
    done[0] += 1
    now = time.time()
    with open(OUT / f"{T}_{R}_{s}.pkl", "wb") as fh:
        pickle.dump({"sample": s, "classes": dict(classes), "parse_s": now - t_last[0]}, fh)
    print(s, "parse", round(now - t_last[0], 1), "s classes", len(classes), flush=True)
    counts = {}
    for key, n in classes.items():
        for t in key:
            counts[t] = counts.get(t, 0) + n / len(key)
    tot = sum(classes.values())
    t_last[0] = time.time()
    return em.Fit(counts, {t: v / tot for t, v in counts.items()}, 0, True)


em.expected_counts = rec
taxa.build_table_from_per_read(reports, OUT / f"tbl_{T}_{R}", f"{TX}/{db}.nodes.dmp", f"{TX}/{db}.names.tsv",
                               rank=R, drop=("Viridiplantae", "Metazoa"), namespace=f"{db}_dump", assign="em", **thr)
