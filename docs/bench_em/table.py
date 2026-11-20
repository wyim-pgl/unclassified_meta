"""Per-sample markdown table from stats/agg/cmp_counts TSVs (run in the benchmark folder)."""
import pandas as pd

st = pd.read_csv("stats_RNA_S.tsv", sep="\t").set_index("sample")
b = pd.read_csv("base_RNA_S_stages.tsv", sep="\t").set_index("sample")
o = pd.read_csv("opt3_RNA_S_stages.tsv", sep="\t").set_index("sample")
c = pd.read_csv("cmp_counts_RNA_S.tsv", sep="\t").set_index("sample")
rows = []
for s in st.index:
    rows.append(
        f"| {s.split('_')[0]} | {st.loc[s, 'reads']:,} | {st.loc[s, 'targets']:,} | {st.loc[s, 'ambiguous_reads']:,} "
        f"| {st.loc[s, 'ambiguous_classes']:,} | {st.loc[s, 'mean_candidates_ambiguous']:.2f} | {st.loc[s, 'crossing']} "
        f"| {int(b.loc[s, 'iter']):,} / {int(o.loc[s, 'iter']):,} | {b.loc[s, 'em']:.0f} / {o.loc[s, 'em']:.0f} "
        f"| {b.loc[s, 'imp']:.0f} / {o.loc[s, 'imp']:.0f} "
        f"| {(b.loc[s, 'em'] + b.loc[s, 'imp']) / (o.loc[s, 'em'] + o.loc[s, 'imp']):.1f} "
        f"| {b.loc[s, 'rss_MB']:.0f} / {o.loc[s, 'rss_MB']:.0f} | {c.loc[s, 'll_new_minus_base']:+.2f} "
        f"| {c.loc[s, 'max_abs_diff_taxa']:.2f} |")
tb, to = b.loc["TOTAL"], o.loc["TOTAL"]
rows.append(f"| **total** | {st['reads'].sum():,} | | {st['ambiguous_reads'].sum():,} | | | "
            f"| {int(tb['iter']):,} / {int(to['iter']):,} | {tb['em']:.0f} / {to['em']:.0f} "
            f"| {tb['imp']:.0f} / {to['imp']:.0f} | {(tb['em'] + tb['imp']) / (to['em'] + to['imp']):.1f} | | | |")
open("per_sample_table.md", "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
print("base stages:", {k: round(tb[k]) for k in ["em:estep", "em:loglik", "imp:SIR_tree_theta", "imp:gibbs_draw"]})
print("opt stages:", {k: round(to[k]) for k in ["em:estep", "em:loglik", "imp:gibbs_draw"]})
