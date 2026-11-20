"""Stage timing table from run_em.py result pickles: agg.py <results dir>"""
import glob
import pickle
import sys

rows = []
for f in sorted(glob.glob(sys.argv[1] + "/*.pkl"), key=lambda p: [int(x) if x.isdigit() else x
                                                                    for x in p.rsplit("/", 1)[1].split("_")]):
    r = pickle.load(open(f, "rb"))
    e, i = r["em_stage"], r["imp_stage"]
    rows.append((r["sample"], r["parse_s"], r["t_em"], r["iterations"], r["t_imp"], e.get("squarem", 0),
                 e.get("estep", 0), e.get("loglik", 0), e.get("gap", 0), e.get("tree_flow", 0),
                 e.get("dominated", 0) + e.get("class_components", 0), i.get("chain", 0), i.get("tree_theta", 0),
                 i.get("draw_assignments", 0), i.get("tree_flow", 0), r["maxrss_kb"] / 1024))
cols = ["sample", "parse", "em", "iter", "imp", "em:squarem", "em:estep", "em:loglik", "em:gap", "em:laminar",
        "em:setup", "imp:chain", "imp:SIR_tree_theta", "imp:gibbs_draw", "imp:laminar_flow", "rss_MB"]
print("\t".join(cols))
for r in rows:
    print("\t".join(str(round(x, 1)) if isinstance(x, float) else str(x) for x in r))
tot = [sum(r[k] for r in rows) for k in range(1, 15)]
print("TOTAL\t" + "\t".join(str(round(x, 1)) for x in tot))
