"""Render the simulation metrics, timing and audit tables for REPORT.md."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
ORDER = ["LCA", "EM", "MI"]


def mdtable(headers, rows):
    return "\n".join("| " + " | ".join(map(str, row)) + " |" for row in [headers, ["---"] * len(headers), *rows])


def wilson(k, n):
    z = 1.959964
    p = k / n
    center = (p + z*z/(2*n))/(1+z*z/n)
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n))/(1+z*z/n)
    return f"{100*max(0,center-half):.1f}–{100*min(1,center+half):.1f}%"


def render():
    f = pd.read_csv(ROOT / "results.tsv", sep="\t")
    rep = pd.read_csv(ROOT / "replicate_metrics.tsv", sep="\t")
    rt = pd.read_csv(ROOT / "runtimes.tsv", sep="\t")
    m = pd.read_csv(ROOT / "metrics.tsv", sep="\t").set_index(["scenario", "method", "subset"])
    lines = []
    def header(title):
        lines.extend(["", "### " + title, ""])
    def metric(scenario, method, subset, stat):
        return m.loc[(scenario, method, subset), stat]
    def measure(scenario, method, subset, stat, pct=False):
        scale = 100 if pct else 1
        value = scale * metric(scenario, method, subset, stat)
        se = scale * metric(scenario, method, subset, stat + "_mcse")
        digits = 2 if pct else 3
        return f"{value:.{digits}f} ± {se:.{digits}f}"

    strict = pd.read_csv(ROOT/"strict_metrics.tsv",sep="\t").set_index(["scenario","method","subset"])
    exact = f[f.strict_revision]
    counts = exact.groupby("scenario").dataset.nunique()
    header("Primary comparison at the requested revision")
    lines.append(f"This subset contains {counts['alt']} alternative and {counts['global_null']} global-null datasets. It includes completed fits from before the shared checkout's first later commit and full CLI runs using the verified local export of cabe340. The extended tables below include {f.dataset.nunique()} distinct datasets in all. The source-cohort field in results.tsv identifies every row.\n")
    rows = []
    for method in ORDER:
        a, n, z, unchanged = (strict.loc["alt",method,"changed"],strict.loc["global_null",method,"all"],
                              strict.loc["alt",method,"all"],strict.loc["alt",method,"unchanged"])
        rows.append([method, f"{a.bias:+.3f} ± {a.bias_mcse:.3f}", f"{a.rmse:.3f}",
                     f"{unchanged.bias:+.3f} / {unchanged.rmse:.3f}", f"{100*a.coverage:.2f}",
                     f"{100*a.bh:.2f}", f"{100*z.fdp:.2f}", f"{100*n.p_lt05:.2f}", f"{100*n.bh:.2f}"])
    lines.append(mdtable(["Method","Changed bias ± MCSE","Changed RMSE","Unchanged bias / RMSE",
                          "Changed coverage (%)","BH power (%)","Alternative FDR (%)","Global-null raw FPR (%)",
                          "Global-null BH FPR (%)"],rows))
    rows=[]
    ss = exact[exact.scenario.eq("global_null")]
    for method in ORDER:
        z=ss[ss.method.eq(method)].groupby("dataset").bh_call.any()
        a=exact[exact.scenario.eq("alt") & exact.method.eq(method)]
        rows.append([method, f"{int(z.sum())}/{len(z)}",wilson(int(z.sum()),len(z)),
                     f"{100*ss[ss.method.eq(method)].covered_centered.mean():.2f}",
                     f"{a[a.prob_ambig.gt(0)].fmi.mean():.4f}" if method=="MI" else "—",
                     f"{a[a.prob_ambig.eq(0)].fmi.mean():.4f}" if method=="MI" else "—"])
    lines.append("\n"+mdtable(["Method","Null datasets with any BH discovery","95% Wilson CI","Global-null coverage (%)","Mean FMI, ambiguous (alt)","Mean FMI, unique (alt)"], rows))
    header("Extended comparison across all audited datasets")
    lines.append("The remaining numerical tables use all accepted datasets, with each seed counted once. The revision audit records agreement of the pooled outputs with cabe340 applied to the saved component fits.")
    header("Effect estimates under the alternative")
    lines.append("Bias is estimated minus true log fold change after centering. Each ± value is one dataset-level Monte Carlo standard error. RMSE is in natural-log units. Up and down rows expose cancellation in pooled signed bias.\n")
    rows = []
    for method in ORDER:
        for subset in ["changed", "up", "down", "unchanged", "ambiguous_unchanged"]:
            a = m.loc["alt", method, subset]
            rows.append([method, subset.replace("_", " "), int(a.estimable), measure("alt", method, subset, "bias"),
                         f"{a.rmse:.3f}", measure("alt", method, subset, "coverage", True),
                         f"{100*a.coverage_raw:.2f}", f"{2*a.halfwidth:.3f}"])
    lines.append(mdtable(["Method", "Species", "Intervals", "Bias ± MCSE", "RMSE", "Centered coverage (%) ± MCSE",
                          "Raw coverage (%)", "Mean CI width"], rows))
    lines.append("\nEach changed species has an unchanged sister. The next table gives the mean centered error per changed species across alternative datasets. This distinguishes species-specific bias from cancellation between species.\n")
    rows = []
    changed = f[(f.scenario == "alt") & f.true_lfc.ne(0)]
    for fid, z in changed.groupby("feature_id"):
        taxid = fid.split("|")[1]
        paired = z.pivot(index="dataset", columns="method", values="error")
        delta = paired.MI - paired.LCA
        rows.append([taxid, f"{z.true_lfc.iloc[0]:+.3f}", f"{z.prob_ambig.iloc[0]:.2f}",
                     *[f"{z.loc[z.method==method,'error'].mean():+.3f} ± {z.loc[z.method==method,'error'].sem():.3f}" for method in ORDER],
                     f"{delta.mean():+.3f} ± {delta.sem():.3f}"])
    lines.append(mdtable(["TaxID", "True LFC", "Ambiguity probability", "LCA bias ± MCSE", "EM bias ± MCSE",
                          "MI bias ± MCSE", "MI − LCA, paired ± MCSE"], rows))
    rows = []
    for method in ORDER:
        sisters = f[(f.scenario == "alt") & f.method.eq(method) & f.changed_sister_lfc.notna()].copy()
        aligned = (sisters.error * np.sign(sisters.changed_sister_lfc)).groupby(sisters.dataset).mean()
        rows.append([method, len(sisters), f"{aligned.mean():+.3f} ± {aligned.sem():.3f}",
                     f"{100*sisters.bh_call.mean():.2f}", f"{100*sisters.covered_centered.mean():.2f}"])
    lines.append("\nAmong the 14 unchanged species per dataset that share reads with a changed sister, positive signed error below means an apparent shift in the direction of that sister's true effect.\n")
    lines.append(mdtable(["Method", "Species-datasets", "Sister-direction error ± MCSE", "BH FPR (%)", "Centered coverage (%)"], rows))
    header("Discoveries and global-null calibration")
    rows = []
    for method in ORDER:
        g = rep.query('scenario == "global_null" and subset == "all" and method == @method')
        k, n = int(g.any_bh.sum()), len(g)
        rows.append([method, measure("alt", method, "changed", "bh", True),
                     measure("alt", method, "all", "fdp", True), measure("global_null", method, "all", "p_lt05", True),
                     measure("global_null", method, "all", "bh", True),
                     f"{100*k/n:.2f} ({k}/{n}, 95% CI {wilson(k,n)})",
                     measure("global_null", method, "all", "coverage", True)])
    lines.append(mdtable(["Method", "Power (%)", "Alternative FDR (%)", "Global-null raw FPR (%)",
                          "Global-null BH FPR (%)", "Global-null FDR, any BH discovery", "Global-null centered coverage (%)"], rows))
    lines.append("\nGlobal-null FDR uses a Wilson binomial interval across datasets. Other ± values are one Monte Carlo standard error in percentage points. The BH FPR denominator is all true-null species. Its value can be small even when several datasets have a discovery.")
    rows = []
    for method in ORDER:
        a, d, n = (m.loc["alt", method, "all"], m.loc["alt", method, "changed"], m.loc["global_null", method, "all"])
        rows.append([method, f"{100*d.robust:.2f}", f"{100*a.robust_fdp:.2f}", f"{100*n.robust_fdp:.2f}", f"{100*n.coverage_raw:.2f}"])
    lines.append("\n" + mdtable(["Method", "Robust power (%)", "Robust alternative FDR (%)", "Robust global-null FDR (%)",
                                 "Raw global-null CI coverage (%)"], rows))
    rows = []
    for method in ORDER:
        z = f[(f.scenario == "alt") & f.method.eq(method) & f.true_lfc.ne(0)]
        rows.append([method, f"{100*z.p.lt(.05).mean():.2f}", f"{z.p.median():.4f}", f"{z.q.median():.3f}",
                     f"{z.se.median():.3f}"])
    lines.append("\n" + mdtable(["Method", "Changed species with raw p < 0.05 (%)", "Median changed-species p",
                                 "Median changed-species BH q", "Median changed-species SE"], rows))
    mi_alt = f[(f.scenario == "alt") & f.method.eq("MI") & f.true_lfc.ne(0)].dropna(subset=["error", "se", "df"])
    matched = (mi_alt.error.abs() <= 1.96 * mi_alt.se).mean()
    lines.append(f"\nUsing the same 1.96 cutoff on the pooled MI SE gives {matched:.2%} centered coverage for changed species, compared with {mi_alt.covered_centered.mean():.2%} using the specified MI t interval. This diagnostic isolates the contribution of the critical value to coverage.")
    header("Assignment uncertainty and paired comparisons")
    rows = []
    for arm in ["alt", "global_null"]:
        subset = f[(f.scenario == arm) & (f.method == "MI")]
        for label, mask in [("all ambiguous", subset.prob_ambig.gt(0)), ("quartet", subset.class_type.eq("quartet")),
                            ("pair", subset.class_type.eq("pair")), ("unique only", subset.class_type.eq("unique"))]:
            z = subset[mask]
            rows.append([arm, label, len(z), f"{z.fmi.mean():.4f}", f"{z.fmi.median():.4f}",
                         f"{z.fmi.quantile(.9):.4f}", f"{z.df.median():.2f}"])
    lines.append(mdtable(["Scenario", "Species", "Species-datasets", "Mean FMI", "Median FMI", "90th pct FMI", "Median MI df"], rows))
    lines.append("\nFMI is the share of pooled variance attributable to between-imputation differences. A uniquely identifiable species can have nonzero FMI because ANCOM-BC2 estimates its coefficient using bias corrections learned across species.")
    rows = []
    for other in ["LCA", "EM"]:
        for arm, sub, stat in [("alt", "changed", "bh"), ("alt", "changed", "coverage"),
                               ("alt", "all", "fdp"), ("global_null", "all", "any_bh")]:
            g = rep[(rep.scenario == arm) & (rep.subset == sub)].pivot(index="dataset", columns="method", values=stat)
            delta = (g.MI - g[other])*100
            rows.append([f"MI minus {other}", arm, {"bh": "power", "coverage": "changed CI coverage",
                         "fdp": "FDR", "any_bh": "global-null FDR"}[stat],
                         f"{delta.mean():+.2f} ± {delta.sem():.2f}"])
    lines.append("\n" + mdtable(["Paired comparison", "Scenario", "Metric", "Difference, percentage points ± MCSE"], rows))
    z = f[(f.scenario == "alt") & f.true_lfc.ne(0)].pivot(index=["dataset", "feature_id"], columns="method",
                                                           values=["lfc_centered", "se"])
    gap = (z.lfc_centered.MI - z.lfc_centered.EM).abs()
    lines.append(f"\nFor changed species the median absolute difference between the centered EM and MI coefficients was {gap.median():.4f} log units. The median ratio of MI SE to EM SE was {(z.se.MI/z.se.EM).median():.4f}.")
    header("Runtime and data checks")
    ds = f[["dataset", "shared_reads"]].drop_duplicates()
    shared = ds.shared_reads.sum()
    total = 2_400_000 * len(ds)
    ok = rt[rt.status == "ok"]
    complete_mi = int(f.loc[f.method == "MI", "n_fits"].eq(5).sum())
    lines.append(f"The {len(ds)} completed datasets contained {total:,} reads. Of these, {shared:,} ({100*shared/total:.2f}%) had multiple candidate species. The accepted datasets represent {7*len(ds)} R fits, excluding repeated and failed attempts. Across methods {f.covered_centered.count():,} of {len(f):,} species-method rows had usable intervals. All {complete_mi:,} MI rows with five estimable fits can be identified in results.tsv.\n")
    rows = []
    for col, label in [("generation_s", "Generate per-read lines"), ("table_lca_s", "LCA table"),
                       ("table_em_s", "EM table + 5 draws"), ("da_lca_s", "LCA DA"),
                       ("da_em_s", "EM expected DA + 5 imputed DA fits")]:
        v = ok[col] / 60
        rows.append([label, f"{v.median():.2f}", f"{v.min():.2f}–{v.max():.2f}", f"{v.sum():.1f}"])
    lines.append(mdtable(["Stage", "Median min/dataset", "Range (min)", "Cumulative min"], rows))
    start = pd.to_datetime(rt.started_utc, utc=True)
    finish = (start + pd.to_timedelta(rt.elapsed_s, unit="s")).max()
    first = min([start.min(), *[pd.Timestamp(json.loads(p.read_text())["started_utc"]) for p in (ROOT/"runs").glob("*/initial_run.json")]])
    lines.append(f"\nElapsed wall time from first generation to last dataset result was {(finish-first).total_seconds()/60:.1f} min ({first:%H:%M:%S} to {finish:%H:%M:%S} UTC). Dataset jobs overlapped, with up to eight slots and one thread per fit. Stage times for the accepted executions include subprocess startup and concurrent machine load. Superseded and failed attempts are recorded separately in the dataset JSON logs.")
    lines.append("\nThe analysis checked every LCA cell against singleton truth, expected-table library sizes against 200,000 (rounding tolerance 0.00004), integer-draw totals, unique-species counts and between-imputation variance against the five saved fits. The unrounded bias-corrected coefficients and original p/q values are in results.tsv.")
    return "\n".join(lines).strip() + "\n"


def update_report():
    path = ROOT / "REPORT.md"
    doc = path.read_text()
    begin = "<!-- generated-results -->"
    end = "<!-- end-generated-results -->"
    if begin in doc:
        before, rest = doc.split(begin, 1)
        _, after = rest.split(end, 1)
    else:
        before, after = doc.split("\n", 1)
        before += "\n\n"
    runs = pd.read_csv(ROOT / "runtimes.tsv", sep="\t")
    counts = runs.loc[runs.status.eq("ok"), "scenario"].value_counts()
    nalt, nnull = counts["alt"], counts["global_null"]
    section = (f"## Results for {nalt} alternative and {nnull} global-null datasets\n\n"
               + render() + "\nA fixed-size replay of these seeds uses "
               + f"`python simulate.py --alt 1 {nalt} --null 1 {nnull} --workers 8` in a fresh study directory.\n")
    tmp = path.with_suffix(".md.tmp")
    tmp.write_text(before + begin + "\n" + section + end + after)
    tmp.replace(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--update", action="store_true", help="refresh the results section in REPORT.md")
    if parser.parse_args().update:
        update_report()
    else:
        print(render(), end="")
