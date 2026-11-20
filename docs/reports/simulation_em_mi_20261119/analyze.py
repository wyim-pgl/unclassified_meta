"""Audit read allocation, extract DA outputs and measure performance against known truth."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import datetime as dt

import numpy as np
import pandas as pd
from unclassified_meta.da import read_results

ROOT = Path(__file__).resolve().parent
METHODS = {"LCA": "da_lca/ancombc2_results.tsv", "EM": "da_em/ancombc2_results.tsv",
           "MI": "da_em/mi_results.tsv"}
TERM = "grouptreated"


def read_counts(path):
    return pd.read_csv(path, sep="\t", index_col="feature_id").sort_index()


def audit(folder, products):
    origin = read_counts(folder / "origin.tsv")
    unique = read_counts(folder / "singleton.tsv")
    lca = read_counts(products / "lca/counts.tsv")
    em = read_counts(products / "em/counts.tsv")
    assert origin.shape == (40, 12) and (origin.sum() == 200_000).all()
    pd.testing.assert_frame_equal(unique, lca, check_dtype=False)
    pure = unique.index[(unique == origin).all(axis=1)]
    assert len(pure) == 10
    led = pd.read_csv(products / "em/ledger.tsv", sep="\t").set_index("sample")
    shared = origin.sum() - unique.sum()
    np.testing.assert_array_equal(shared, led.loc[shared.index, "reads_shared_by_em"])
    np.testing.assert_allclose(em.sum(), origin.sum(), atol=0.00004, rtol=0)
    np.testing.assert_array_equal(em.loc[pure], origin.loc[pure])
    draws = sorted((products / "em/imputations").glob("counts_*.tsv"))
    assert len(draws) == 5
    for draw in draws:
        im = read_counts(draw)
        assert (im >= unique).all().all() and (im % 1 == 0).all().all()
        np.testing.assert_array_equal(im.sum(), origin.sum())
        np.testing.assert_array_equal(im.loc[pure], origin.loc[pure])
    return int(shared.sum())


def extract(folder, rec, products):
    shared = audit(folder, products)
    truth = pd.read_csv(folder / "truth.tsv", sep="\t")
    truth["class_type"] = np.select([truth.prob_ambig == 0, truth.prob_ambig == 0.5], ["unique", "pair"], default="quartet")
    i = truth.feature_id.str.split("|", regex=False).str[1].astype(int) - 100
    component = (i // 8) * 8 + np.select([i % 8 < 4, i % 8 < 6], [0, 4], default=i % 8)
    nonzero_component = truth.true_lfc.where(truth.true_lfc.ne(0)).groupby(component).transform("max")
    truth["changed_sister_lfc"] = nonzero_component.where(truth.true_lfc.eq(0))
    frames = []
    for method, path in METHODS.items():
        source = read_results(products / path)
        assert source.feature_id.is_unique and set(source.feature_id) <= set(truth.feature_id)
        f = truth.merge(source, how="left", validate="1:1", on="feature_id")
        df = truth.copy()
        df["dataset"], df["scenario"], df["rep"], df["method"] = folder.name, rec["scenario"], rec["rep"], method
        df["source_cohort"] = rec["source_cohort"]
        df["strict_revision"] = rec["strict_revision"]
        df["shared_reads"] = shared
        for col in ["lfc", "se", "p", "q", "diff", "diff_robust", "df", "fmi", "n_fits"]:
            k = col + "_" + TERM
            if k not in f and col not in {"df", "fmi", "n_fits"}:
                raise ValueError(f"{folder}/{path}: missing {k}")
            df[col] = f[k] if k in f else np.nan
        if method == "MI":
            draws = sorted((products / "da_em/imputations").glob("ancombc2_results_*.tsv"))
            assert len(draws) == 5
            fits = [read_results(path).set_index("feature_id").reindex(truth.feature_id) for path in draws]
            est = np.stack([r["lfc_" + TERM] for r in fits])
            se = np.stack([r["se_" + TERM] for r in fits])
            mean = est.mean(axis=0)
            within = (se ** 2).mean(axis=0)
            between = est.var(axis=0, ddof=1)
            between[between <= 1e-12 * within] = 0
            total = within + 1.2 * between
            ok = df.lfc.notna()
            assert df.loc[ok, "n_fits"].eq(5).all()
            np.testing.assert_allclose(df.loc[ok, "lfc"], mean[ok], atol=1e-12)
            np.testing.assert_allclose(df.loc[ok, "se"] ** 2, total[ok], atol=1e-12)
            np.testing.assert_allclose(df.loc[ok, "fmi"], (1.2 * between / total)[ok], atol=1e-12)
        df["offset"] = df.lfc.where(df.true_lfc == 0).median()
        df["lfc_centered"] = df.lfc - df.offset
        df["error"] = df.lfc_centered - df.true_lfc
        df["bh_call"] = df.q.lt(.05)
        df["robust_call"] = df.diff_robust.astype(str).str.lower().isin(["true", "1", "1.0"])
        assert (df.robust_call <= df.bh_call).all()
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def add_intervals(f):
    critical = np.full(len(f), 1.96)
    mi = f.method.eq("MI")
    df = f.loc[mi, "df"]
    # R's Student quantile supports noninteger and infinite degrees of freedom.
    values = df.to_csv(index=False, header=False, na_rep="NA").replace("inf", "Inf")
    out = subprocess.run([os.environ["UNCLASSIFIED_META_RSCRIPT"], "--vanilla", "-e",
                          'options(digits=17); d <- scan(file("stdin"), quiet=TRUE); write(qt(.975,d),stdout(),ncolumns=1)'],
                         input=values, capture_output=True, text=True, check=True)
    critical[mi] = pd.to_numeric(pd.Series(out.stdout.split()), errors="coerce")
    f["halfwidth"] = critical * f.se
    for col, est in [("covered_centered", f.lfc_centered), ("covered_raw", f.lfc)]:
        valid = np.isfinite(est) & np.isfinite(f.halfwidth) & (f.se > 0)
        f[col] = ((est - f.true_lfc).abs() <= f.halfwidth).astype(float).where(valid)


def by_dataset(f):
    rows = []
    for (dataset, scenario, method), x in f.groupby(["dataset", "scenario", "method"]):
        signal = x.true_lfc != 0
        for name, mask in {"all": x.feature_id.notna(), "changed": signal, "unchanged": ~signal,
                           "up": x.true_lfc > 0, "down": x.true_lfc < 0,
                           "ambiguous": x.prob_ambig > 0, "unique": x.prob_ambig == 0,
                           "ambiguous_unchanged": ~signal & (x.prob_ambig > 0),
                           "null_sister": x.changed_sister_lfc.notna(),
                           "quartet_changed": signal & x.class_type.eq("quartet"),
                           "pair_changed": signal & x.class_type.eq("pair")}.items():
            s = x[mask]
            if s.empty:
                continue
            rows.append(dict(dataset=dataset, scenario=scenario, method=method, subset=name, n=len(s),
                             n_estimable=int(s.covered_centered.count()), bias=s.error.mean(), mse=s.error.pow(2).mean(),
                             coverage=s.covered_centered.mean(), coverage_raw=s.covered_raw.mean(),
                             halfwidth=s.halfwidth.mean(), fmi=s.fmi.mean(), mean_df=s.df.mean(),
                             p_lt05=s.p.lt(.05).mean(), bh=s.bh_call.mean(), robust=s.robust_call.mean(),
                             fdp=(s.bh_call & s.true_lfc.eq(0)).sum() / max(1, s.bh_call.sum()),
                             robust_fdp=(s.robust_call & s.true_lfc.eq(0)).sum() / max(1, s.robust_call.sum()),
                             any_bh=int(s.bh_call.any()), offset=s.offset.iloc[0]))
    return pd.DataFrame(rows)


def summarize(rep):
    rows = []
    for key, g in rep.groupby(["scenario", "method", "subset"], sort=False):
        row = dict(zip(["scenario", "method", "subset"], key))
        row.update(datasets=len(g), species_datasets=g.n.sum(), estimable=g.n_estimable.sum())
        for stat in ["bias", "mse", "coverage", "coverage_raw", "halfwidth", "fmi", "p_lt05", "bh",
                     "robust", "fdp", "robust_fdp", "any_bh", "offset"]:
            row[stat] = g[stat].mean()
            row[stat + "_mcse"] = g[stat].sem()
        row["rmse"] = np.sqrt(row["mse"])
        rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--completed-only", action="store_true", help="interim summary while datasets are running")
    parser.add_argument("--reuse-checked", action="store_true", help="reuse already audited rows from the current study checkpoint")
    opt = parser.parse_args()
    prior_rows, prior_reps = {}, {}
    if opt.reuse_checked:
        for name, cache in [("results",prior_rows), ("replicate_metrics",prior_reps)]:
            old = pd.read_csv(ROOT/(name+".tsv"),sep="\t",float_precision="round_trip")
            old["scenario"] = old.scenario.replace({"global_null":"null"})
            cache.update(dict(tuple(old.groupby("dataset"))))
    retained = set()
    records, frames = [], []
    for folder in sorted((ROOT / "runs").iterdir()):
        if opt.completed_only and not (folder / "run.json").exists():
            continue
        rec = json.loads((folder / "run.json").read_text())
        base = folder
        rec["original_status"] = rec["status"]
        original_end = dt.datetime.fromisoformat(rec["started_utc"]) + dt.timedelta(seconds=rec["elapsed_s"])
        # First later commit on the shared branch, with only documentation changes.
        initial_checkout = original_end < dt.datetime.fromisoformat("2026-11-19T18:57:02+00:00")
        pinned = str(ROOT / "pinned/src")
        fresh_pinned = rec.get("pythonpath") == pinned or (folder / "pinned_launch.txt").is_file()
        replay_record = folder / "pinned_fit/run.json"
        if replay_record.exists():
            replay = json.loads(replay_record.read_text())
            if replay["status"] == "ok":
                base = folder / "pinned_fit"
                rec["status"] = "ok"
                rec["commands"] = replay["commands"]
                replay_end = dt.datetime.fromisoformat(replay["started_utc"]) + dt.timedelta(seconds=replay["elapsed_s"])
                rec["elapsed_s"] = (replay_end - dt.datetime.fromisoformat(rec["started_utc"])).total_seconds()
        rec["strict_revision"] = initial_checkout or fresh_pinned or base != folder
        if base != folder:
            cohort = "pinned_replay"
        elif fresh_pinned:
            cohort = "pinned_fresh"
        elif initial_checkout:
            cohort = "initial_checkout"
        else:
            cohort = "later_checkout"
        rec["source_cohort"] = cohort
        records.append({"dataset": folder.name, **{k: v for k, v in rec.items() if k != "commands"},
                        **{c["label"] + "_s": c["seconds"] for c in rec["commands"]}})
        if rec["status"] == "ok":
            prior = prior_rows.get(folder.name)
            if (prior is not None and prior.source_cohort.eq(cohort).all()
                    and prior.strict_revision.eq(rec["strict_revision"]).all()
                    and (base/"repool_check.json").exists() and len(prior)==120):
                retained.add(folder.name)
                frames.append(prior)
            else:
                frames.append(extract(folder, rec, base))
    f = pd.concat(frames, ignore_index=True)
    add_intervals(f)
    new_rows = f[~f.dataset.isin(retained)]
    rep = pd.concat([by_dataset(new_rows), *[prior_reps[k] for k in sorted(retained)]], ignore_index=True)
    primary = rep[rep.dataset.isin(f.loc[f.strict_revision,"dataset"].unique())]
    for name, data in [("strict_metrics", summarize(primary)), ("results", f), ("replicate_metrics", rep), ("metrics", summarize(rep)),
                       ("runtimes", pd.DataFrame(records))]:
        data = data.replace({"scenario": {"null": "global_null"}})
        data.to_csv(ROOT / (name + ".tsv"), sep="\t", index=False, na_rep="NA")
    print(f"{len(records)} dataset records, {len(frames)} complete, {len(f)} species/method results")
