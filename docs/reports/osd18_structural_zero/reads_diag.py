"""Unique vs shared reads for target species (rank S), from the per-read Centrifuger files used by taxa-table."""
import sys, csv, collections
sys.path.insert(0, "/data/gpfs/assoc/pgl/develop/unclassified_meta_test/em_code_fix/src")
from unclassified_meta.extract import centrifuge_reads
from unclassified_meta.taxa import Taxonomy, scan_names, RANK_NAMES
TX = "/data/gpfs/assoc/pgl/develop/unclassified_meta_test/OSD18"
nodes, names, reports, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
TARGETS = {"3432864": "Tobamovirus cacti", "1297894": "Citrus yellow vein-associated virus",
           "3393506": "Tobamovirus muricaudae", "352851": "Sclerotinia nivalis", "1187928": "Alternaria incomplexa"}
tax = Taxonomy(nodes)
_, hits = scan_names(names, want_names={"viridiplantae", "metazoa"})
drop = set().union(*hits.values())
target_ranks = RANK_NAMES["S"]
cache = {}
def tok(t):
    h = cache.get(t)
    if h is None:
        p = tax.path(t)
        if not p: h = "unknown"
        elif drop.intersection(p): h = "dropped"
        else:
            at = next((x for x in p if tax.rank.get(x) in target_ranks), None)
            h = at if at else "above_rank"
        cache[t] = h
    return h
rows, partners = [], collections.defaultdict(collections.Counter)
samples = [l.rstrip("\n").split("\t") for l in open(reports)][1:]
for s, path in samples:
    uniq, shared, nset = collections.Counter(), collections.Counter(), collections.Counter()
    for _, status, taxids in centrifuge_reads(path, 64, 0):
        if status != 2: continue
        ks = {tok(t) for t in taxids}
        hit = ks.intersection(TARGETS)
        for t in hit:
            if len(ks) == 1: uniq[t] += 1
            else:
                shared[t] += 1; nset[t] += len(ks)
                for o in ks - {t}: partners[(t, s)][o] += 1
    for t in TARGETS:
        rows.append([s, t, TARGETS[t], uniq[t], shared[t], round(nset[t] / shared[t], 2) if shared[t] else 0,
                     ";".join(f"{o}:{n}" for o, n in partners[(t, s)].most_common(4))])
    print(s, {TARGETS[t]: (uniq[t], shared[t]) for t in TARGETS}, flush=True)
with open(out, "w") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(["sample", "taxid", "name", "unique_reads", "shared_reads", "mean_candidates_when_shared", "top_partners"])
    w.writerows(rows)
