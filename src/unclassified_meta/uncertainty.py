"""em-uncertainty: how much of a panel feature's count is observed, how much is assigned, and how stable both are.

Input: a sealed `panel-quant` folder (any --ambiguous mode). Its `equivalence_classes.tsv` holds every fragment of
every sample once, grouped by its set of candidate features/decoys (empty = unmapped); `unique_counts.tsv` holds
the fragments whose best alignments fall in one feature only.

Per sample:
  - EM point: `em.expected_counts` of the mapped classes, as `panel-quant --ambiguous em` computes it.
  - Read bootstrap (`--bootstraps` B): the sample's fragments are resampled with replacement, i.e. one multinomial
    draw of (input fragments) over the equivalence classes in proportion to their fragments, the unmapped class
    included, so the replicate keeps the input fragment total. A fragment is a read pair, so mates stay together.
    The EM is refitted on every replicate. RNG: `numpy.random.default_rng([seed, 3, sample_index])`, one
    multinomial draw per replicate in replicate order over the classes sorted by candidate tuple (sample_index =
    the sample's column position in the quant's unique_counts.tsv). Each sample is one task, so the results do
    not depend on --threads.
  - Realizations (`--realizations` R): integer draws of the posterior assignment, `em.imputations(mapped, fit,
    R, default_rng([seed, 4, sample_index]))` (stream 4; panel-quant's imputations use stream 2).
  - Alternative mapping criteria per feature: unique (single-candidate classes), even split (each shared fragment
    split equally among its candidates, decoys included), EM, any-candidate (each shared fragment counted for every
    one of its candidates: an upper bound). unique <= EM <= any_candidate always; even split lies between unique
    and any_candidate.

Limits: the bootstrap and the realizations measure read-sampling and assignment instability of ONE library. They
are not biological replicates, add no samples to a test, and their p-values must never be pooled as if they were.
The confidence classes below are descriptive labels with a fixed, printed rule, not probabilities.

Output folder layout (all count tables: `feature_id` first, then one column per sample, features in the order of
features.tsv, decoys excluded; MANIFEST.json is written last and seals every file listed here):
  unique_counts.tsv                byte copy of the quant's table (observed; never replaced by EM values)
  em_expected_counts.tsv           EM expected fragments, 6 decimals
  normalized_abundance.tsv         EM fragments per million input fragments of the sample (input = all classes,
                                   unmapped included), 6 decimals
  aggregated_counts.tsv            EM counts summed by the aggregation key (first column `feature_id` = the key)
  aggregated_unique_counts.tsv     unique counts summed by the same key
  aggregated_groups.tsv            feature_id, group (the key each feature was summed into)
  bootstrap/counts_NNN.tsv         EM counts of bootstrap replicate NNN (001..B), 6 decimals
  realizations/counts_NNN.tsv      integer posterior realization NNN (001..R)
  em_bootstrap_abundance.tsv       long: feature_id, label, sample, em, boot_mean, boot_sd, boot_p2_5, boot_p97_5,
                                   boot_cv (NaN when boot_mean < --min-mean-for-cv), n_bootstraps
  em_assignment_uncertainty.tsv    long: feature_id, label, sample, unique, shared, em, unique_fraction,
                                   em_from_shared, even_split, any_candidate, realization_mean, realization_sd,
                                   realization_p2_5, realization_p97_5
  viral_unique_read_support.tsv    one row per feature, after '#' header lines with the rule (read it with
                                   `read_support`)
  figures/unique_vs_shared_support.png, figures/bootstrap_abundance.png
  MANIFEST.json                    {"step": "em-uncertainty", ...}; consumers call `verify_uncertainty(folder)`.

assignment_confidence (per feature, totals over samples; unique_fraction = unique_total / em_total):
  low       em_driven (unique_fraction < --unique-fraction-cut, default 0.5) or unique_total < 10 or
            boot_cv_median > 0.5 or em_total > 2 * unique_total (EM more than doubles the observed count; at the
            default cut this coincides with em_driven, it stays in force when the cut is changed)
  high      not low, and unique_fraction >= 0.8 and boot_cv_median <= 0.2 and identity_class in (same, reference)
            and not mixed_identity
  moderate  otherwise (also when the bootstrap CV is undefined, i.e. no sample has boot_mean >= 1)
"""
from __future__ import annotations

import argparse
import logging
import os
import platform
import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__, em, manifest, panel
from .fragments import DECOY_PREFIX
from .io_utils import make_staging_dir, refuse_overlap, remove_staging_dir

log = logging.getLogger("unclassified_meta")

STEP = "em-uncertainty"
RNG_DESCRIPTION = ("bootstrap: numpy.random.default_rng([seed, 3, sample_index]), one multinomial(input fragments, "
                   "class proportions) per replicate in replicate order, classes sorted by candidate tuple; "
                   "realizations: em.imputations(mapped, fit, R, numpy.random.default_rng([seed, 4, sample_index])); "
                   "sample_index = column position in the quant's unique_counts.tsv")
CORE_FILES = ("unique_counts.tsv", "em_expected_counts.tsv", "normalized_abundance.tsv", "aggregated_counts.tsv",
              "aggregated_unique_counts.tsv", "aggregated_groups.tsv", "em_bootstrap_abundance.tsv",
              "em_assignment_uncertainty.tsv", "viral_unique_read_support.tsv")
FIGURES = ("figures/unique_vs_shared_support.png", "figures/bootstrap_abundance.png")
HIGH_UNIQUE_FRACTION, HIGH_CV, LOW_CV, MIN_UNIQUE, MAX_FOLD = 0.8, 0.2, 0.5, 10, 2.0
CONFIDENT_IDENTITY = ("same", "reference")


def confidence_rule(cut: float) -> str:
    return (f"assignment_confidence: low = unique_fraction < {cut:g} (em_driven) or unique_total < {MIN_UNIQUE} or "
            f"boot_cv_median > {LOW_CV:g} or em_total > {MAX_FOLD:g} * unique_total; high = not low and "
            f"unique_fraction >= {HIGH_UNIQUE_FRACTION:g} and boot_cv_median <= {HIGH_CV:g} and identity_class in "
            f"({', '.join(CONFIDENT_IDENTITY)}) and not mixed_identity; moderate = otherwise. unique_fraction = "
            "unique_total / em_total (totals over samples). Descriptive classes, not probabilities; bootstrap and "
            "realizations measure read sampling and assignment instability, not biological replication.")


def bootstrap_name(b: int) -> str:
    return f"bootstrap/counts_{b:03d}.tsv"


def realization_name(r: int) -> str:
    return f"realizations/counts_{r:03d}.tsv"


# ---------------------------------------------------------------- per-sample computation (one task per sample)

def resample_classes(classes: dict, rng: np.random.Generator) -> dict:
    """One read-bootstrap replicate of a sample: (input fragments) fragments drawn with replacement, i.e. a
    multinomial over the classes (sorted by candidate tuple, unmapped () included) in proportion to their
    fragments. The replicate keeps the input fragment total."""
    keys = sorted(classes)
    n = np.array([classes[k] for k in keys], dtype=np.int64)
    total = int(n.sum())
    if not total:
        return {}
    draw = rng.multinomial(total, n / total)
    return {k: int(d) for k, d in zip(keys, draw) if d}


def _sample_task(task: tuple) -> dict:
    """EM point, B bootstrap refits and R realizations of one sample. Its own RNG streams make the result
    independent of which process runs it."""
    name, index, classes, seed, bootstraps, realizations = task
    total = int(sum(classes.values()))
    mapped = {k: int(v) for k, v in classes.items() if k and v}
    fit = em.expected_counts(mapped)
    if not fit.converged:
        raise ValueError(f"{name}: the EM assignment did not converge in {fit.iterations} iterations")
    rng = np.random.default_rng([int(seed), 3, int(index)])
    boots, not_converged = [], 0
    for _ in range(bootstraps):
        rep = resample_classes(classes, rng)
        bfit = em.expected_counts({k: d for k, d in rep.items() if k})
        not_converged += not bfit.converged
        boots.append(em.round_preserving(bfit.counts, 6))
    reals = em.imputations(mapped, fit, realizations, np.random.default_rng([int(seed), 4, int(index)]))
    for r, d in enumerate(reals):
        if sum(d.values()) != sum(mapped.values()):
            raise ValueError(f"{name}: realization {r + 1} does not keep every mapped fragment")
    return {"sample": name, "total": total, "em": em.round_preserving(fit.counts, 6), "boot": boots,
            "real": reals, "not_converged": not_converged}


def criteria(classes: dict, features: list) -> dict:
    """unique, shared, even_split and any_candidate fragments per feature (numpy arrays in `features` order)."""
    idx = {f: i for i, f in enumerate(features)}
    out = {k: np.zeros(len(features)) for k in ("unique", "shared", "even_split", "any_candidate")}
    for cls, v in classes.items():
        if not cls or not v:
            continue
        for t in cls:
            i = idx.get(t)
            if i is None:
                continue                                  # a decoy candidate: not a feature row
            if len(cls) == 1:
                out["unique"][i] += v
            else:
                out["shared"][i] += v
                out["even_split"][i] += v / len(cls)
    out["even_split"] += out["unique"]
    out["any_candidate"] = out["unique"] + out["shared"]
    return out


# ---------------------------------------------------------------- inputs

def _truthy(v) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes")


def aggregation_groups(features: pd.DataFrame, aggregate_by: str | None, taxonomy: Path | None,
                       quant_rec: dict) -> tuple[pd.Series, str, dict]:
    """(group key per feature_id, description, {taxonomy path: sha256}). A feature with an empty key stays a group
    of its own (its feature_id)."""
    fid = features["feature_id"].astype(str)
    if taxonomy is not None:
        taxonomy = Path(taxonomy)
        if not (taxonomy.parent / manifest.NAME).exists():
            raise manifest.ManifestError(f"{taxonomy} has no {manifest.NAME} next to it: use the panel_taxonomy.tsv "
                                         "written by `panel-taxonomy`")
        rec = manifest.verify_file(taxonomy, what="panel taxonomy")
        want = (quant_rec.get("panel") or {}).get("panel_id")
        if rec.get("panel_id") and want and rec["panel_id"] != want:
            raise ValueError(f"{taxonomy} was made for another panel (panel_id {rec['panel_id']}) than the quant "
                             f"tables (panel_id {want}); run panel-taxonomy on that panel")
        tx = pd.read_csv(taxonomy, sep="\t", dtype=str, keep_default_na=False)
        if not {"feature_id", "taxid"} <= set(tx.columns):
            raise ValueError(f"{taxonomy}: needs the columns feature_id and taxid")
        taxid = dict(zip(tx["feature_id"], tx["taxid"].str.strip()))
        key = [f"taxid:{taxid[f]}" if taxid.get(f) else f for f in fid]
        return (pd.Series(key, index=fid.values), "panel_taxonomy:taxid",
                {str(taxonomy.resolve()): rec["artifacts"][taxonomy.name]})
    if aggregate_by is None:
        has_parent = "parent_feature" in features.columns and \
            features["parent_feature"].fillna("").astype(str).str.strip().ne("").any()
        aggregate_by = "parent_feature" if has_parent else "label"
    if aggregate_by not in features.columns:
        raise ValueError(f"--aggregate-by {aggregate_by}: not a column of features.tsv "
                         f"({', '.join(features.columns)})")
    vals = features[aggregate_by].fillna("").astype(str).str.strip()
    key = [v if v else f for v, f in zip(vals, fid)]
    return pd.Series(key, index=fid.values), aggregate_by, {}


# ---------------------------------------------------------------- tables

def _matrix(features: list, samples: list, values: dict, decimals: int | None = 6) -> pd.DataFrame:
    df = pd.DataFrame({"feature_id": features})
    for s in samples:
        v = values[s]
        df[s] = [v.get(f, 0) for f in features] if isinstance(v, dict) else v
        if decimals is None:
            df[s] = df[s].astype(np.int64)
    return df


def _write(df: pd.DataFrame, path: Path, decimals: int | None = 6, header: list | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for line in header or []:
            fh.write(f"# {line}\n")
        df.to_csv(fh, sep="\t", index=False, float_format=f"%.{decimals}f" if decimals is not None else None,
                  lineterminator="\n")


def read_support(path: str | Path) -> pd.DataFrame:
    """viral_unique_read_support.tsv without its leading '#' rule lines (labels may contain '#', so the file is
    not read with comment='#')."""
    with open(path) as fh:
        skip = 0
        for line in fh:
            if not line.startswith("#"):
                break
            skip += 1
    return pd.read_csv(path, sep="\t", skiprows=skip, keep_default_na=False, na_values=["", "NaN", "nan"])


def _cv(mean: np.ndarray, sd: np.ndarray, min_mean: float) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(mean >= min_mean, sd / np.where(mean > 0, mean, 1.0), np.nan)


def _median(values: list) -> float:
    a = np.array(values, dtype=float)
    a = a[~np.isnan(a)]
    return float(np.median(a)) if a.size else float("nan")


def classify(unique_fraction: float, unique_total: float, em_total: float, cv: float, identity: str, mixed: bool,
             cut: float) -> tuple[bool, str]:
    """(em_driven, assignment_confidence) by the rule of the module notes."""
    em_driven = bool(not np.isnan(unique_fraction) and unique_fraction < cut)
    if em_driven or unique_total < MIN_UNIQUE or (not np.isnan(cv) and cv > LOW_CV) or em_total > MAX_FOLD * unique_total:
        return em_driven, "low"
    if unique_fraction >= HIGH_UNIQUE_FRACTION and not np.isnan(cv) and cv <= HIGH_CV \
            and identity in CONFIDENT_IDENTITY and not mixed:
        return em_driven, "high"
    return em_driven, "moderate"


def _figures(support: pd.DataFrame, boot_totals: np.ndarray, features: list, folder: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    (folder / "figures").mkdir(parents=True, exist_ok=True)
    order = support.sort_values("em_total", ascending=False, kind="mergesort")
    top = order.head(40).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, max(2.0, 0.22 * len(top) + 1)))
    names = [str(x)[:50] for x in top["label"]]
    ax.barh(names, top["unique_total"], color="#2b6cb0", label="unique fragments")
    ax.barh(names, (top["em_total"] - top["unique_total"]).clip(lower=0), left=top["unique_total"],
            color="#ed8936", label="EM share of shared fragments")
    ax.set_xlabel("fragments (all samples)")
    ax.tick_params(axis="y", labelsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(folder / FIGURES[0], dpi=200)
    plt.close(fig)
    pos = {f: i for i, f in enumerate(features)}
    top = order.head(30).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, max(2.0, 0.24 * len(top) + 1)))
    data = [boot_totals[:, pos[f]] for f in top["feature_id"]]
    if data and boot_totals.shape[0]:
        ax.boxplot(data, vert=False, widths=0.6, showfliers=False)
        for i, d in enumerate(data, start=1):
            ax.plot(d, np.full(len(d), i), "o", ms=1.5, alpha=0.4, color="#2b6cb0")
        ax.set_yticks(range(1, len(data) + 1), [str(x)[:50] for x in top["label"]], fontsize=6)
    ax.set_xlabel("EM fragments per bootstrap replicate (all samples)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(folder / FIGURES[1], dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- run

def run(quant_dir, outdir, bootstraps: int = 100, realizations: int = 20, seed: int = 1, aggregate_by: str | None = None,
        taxonomy=None, threads: int = 1, min_mean_for_cv: float = 1.0, unique_fraction_cut: float = 0.5) -> Path:
    """Write the em-uncertainty folder (layout in the module notes) and return it."""
    quant_dir, outdir = Path(quant_dir), Path(outdir)
    bootstraps, realizations, seed, threads = int(bootstraps), int(realizations), int(seed), max(int(threads), 1)
    if bootstraps < 2:
        raise ValueError("--bootstraps must be at least 2 (a standard deviation needs two replicates)")
    if realizations < 1:
        raise ValueError("--realizations must be at least 1")
    if seed < 0:
        raise ValueError("--seed must be 0 or a positive whole number")
    if min_mean_for_cv < 0 or not 0 <= unique_fraction_cut <= 1:
        raise ValueError("--min-mean-for-cv must be >= 0 and --unique-fraction-cut within [0, 1]")
    if outdir.exists() and outdir.resolve() == quant_dir.resolve():
        raise ValueError(f"{STEP}: the output folder is the quant folder {quant_dir}; write to another folder")
    classes_by_sample, _ = panel.read_class_tables(quant_dir)              # verifies the quant seal
    qrec = panel.verify_quant(quant_dir)
    if "features.tsv" not in qrec["artifacts"]:
        raise manifest.ManifestError(f"{quant_dir}: features.tsv is not sealed by provenance.json; run panel-quant again")
    inputs = [quant_dir / n for n in ("features.tsv", *panel.CLASS_ARTIFACTS, "provenance.json")]
    outputs = [outdir / n for n in (*CORE_FILES, *FIGURES, manifest.NAME, manifest.MARKER)]
    refuse_overlap(outputs, inputs + ([Path(taxonomy)] if taxonomy else []), STEP)

    features = pd.read_csv(quant_dir / "features.tsv", sep="\t", dtype={"feature_id": str}, keep_default_na=False)
    feats = list(features["feature_id"])
    uniq_q = pd.read_csv(quant_dir / "unique_counts.tsv", sep="\t", dtype={"feature_id": str})
    samples = [str(c) for c in uniq_q.columns if c != "feature_id"]
    if list(uniq_q["feature_id"]) != feats:
        raise ValueError(f"{quant_dir}: unique_counts.tsv and features.tsv list different features")
    unknown_samples = sorted(set(classes_by_sample) - set(samples))
    if unknown_samples:
        raise ValueError(f"{quant_dir}: equivalence classes of samples not in unique_counts.tsv: {unknown_samples[:3]}")
    known = set(feats)
    for s, cl in classes_by_sample.items():
        bad = sorted({t for k in cl for t in k if t not in known and not t.startswith(DECOY_PREFIX)})
        if bad:
            raise ValueError(f"{s}: candidates missing from features.tsv: {bad[:3]}")
    groups, agg_desc, tax_sha = aggregation_groups(features, aggregate_by, taxonomy, qrec)

    # ---- compute (each sample one task with its own RNG streams: identical results for any --threads)
    tasks = [(s, i, classes_by_sample.get(s, {}), seed, bootstraps, realizations) for i, s in enumerate(samples)]
    if threads > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=min(threads, len(tasks))) as ex:
            results = list(ex.map(_sample_task, tasks))
    else:
        results = [_sample_task(t) for t in tasks]
    res = {r["sample"]: r for r in results}
    not_conv = sum(r["not_converged"] for r in results)
    if not_conv:
        log.warning("%d bootstrap EM refits did not converge (their certified-gap stop was not reached); they are "
                    "kept and counted in the manifest", not_conv)
    crit = {s: criteria(classes_by_sample.get(s, {}), feats) for s in samples}
    for s in samples:                                       # the observed table and the classes must agree
        if not np.array_equal(crit[s]["unique"], uniq_q[s].to_numpy(dtype=float)):
            raise ValueError(f"{quant_dir}: unique_counts.tsv of {s} does not match equivalence_classes.tsv")
    pos = {f: i for i, f in enumerate(feats)}
    vec = lambda d: np.array([d.get(f, 0.0) for f in feats], dtype=float)      # noqa: E731
    em_v = {s: vec(res[s]["em"]) for s in samples}
    boot = {s: np.array([vec(b) for b in res[s]["boot"]]) for s in samples}     # B x features
    real = {s: np.array([vec(r) for r in res[s]["real"]]) for s in samples}     # R x features
    labels = list(features["label"]) if "label" in features.columns else feats

    # ---- long tables
    ab_rows, as_rows = [], []
    per_sample = {}
    for s in samples:
        b, r, e, c = boot[s], real[s], em_v[s], crit[s]
        bm, bsd = b.mean(axis=0), b.std(axis=0, ddof=1)
        bp = np.percentile(b, [2.5, 97.5], axis=0)
        rm = r.mean(axis=0)
        rsd = r.std(axis=0, ddof=1) if r.shape[0] > 1 else np.full(len(feats), np.nan)
        rp = np.percentile(r, [2.5, 97.5], axis=0)
        bcv, rcv = _cv(bm, bsd, min_mean_for_cv), _cv(rm, rsd, min_mean_for_cv)
        per_sample[s] = (bcv, rcv)
        with np.errstate(divide="ignore", invalid="ignore"):
            ufrac = np.where(e > 0, c["unique"] / np.where(e > 0, e, 1.0), np.nan)
        for i, f in enumerate(feats):
            ab_rows.append((f, labels[i], s, e[i], bm[i], bsd[i], bp[0, i], bp[1, i], bcv[i], b.shape[0]))
            as_rows.append((f, labels[i], s, c["unique"][i], c["shared"][i], e[i], ufrac[i], e[i] - c["unique"][i],
                            c["even_split"][i], c["any_candidate"][i], rm[i], rsd[i], rp[0, i], rp[1, i]))
    abundance = pd.DataFrame(ab_rows, columns=["feature_id", "label", "sample", "em", "boot_mean", "boot_sd",
                                               "boot_p2_5", "boot_p97_5", "boot_cv", "n_bootstraps"])
    assignment = pd.DataFrame(as_rows, columns=["feature_id", "label", "sample", "unique", "shared", "em",
                                                "unique_fraction", "em_from_shared", "even_split", "any_candidate",
                                                "realization_mean", "realization_sd", "realization_p2_5",
                                                "realization_p97_5"])

    # ---- per-feature support
    tot = {k: sum(crit[s][k] for s in samples) for k in ("unique", "shared", "even_split", "any_candidate")}
    em_tot = sum(em_v[s] for s in samples) if samples else np.zeros(len(feats))
    srows = []
    for i, f in enumerate(feats):
        u, e = float(tot["unique"][i]), float(em_tot[i])
        ufrac = u / e if e > 0 else float("nan")
        bcv = _median([per_sample[s][0][i] for s in samples])
        rcv = _median([per_sample[s][1][i] for s in samples])
        ident = str(features["identity_class"].iloc[i]) if "identity_class" in features.columns else ""
        mixed = _truthy(features["mixed_identity"].iloc[i]) if "mixed_identity" in features.columns else False
        driven, conf = classify(ufrac, u, e, bcv, ident, mixed, unique_fraction_cut)
        four = [u, float(tot["even_split"][i]), e, float(tot["any_candidate"][i])]
        row = {"feature_id": f, "label": labels[i], "identity_class": ident, "mixed_identity": mixed}
        if "parent_feature" in features.columns:
            row["parent_feature"] = features["parent_feature"].iloc[i]
        row.update({"unique_total": u, "shared_total": float(tot["shared"][i]), "em_total": e,
                    "even_split_total": four[1], "any_candidate_total": four[3], "unique_fraction": ufrac,
                    "em_driven": driven, "criteria_ratio": four[3] / max(u, 1.0), "criteria_min": min(four),
                    "criteria_max": max(four), "boot_cv_median": bcv, "realization_cv_median": rcv,
                    "samples_with_unique_support": int(sum(crit[s]["unique"][i] > 0 for s in samples)),
                    "assignment_confidence": conf})
        srows.append(row)
    support = pd.DataFrame(srows)

    # ---- matrices
    norm = {s: (em_v[s] / res[s]["total"] * 1e6 if res[s]["total"] else np.zeros(len(feats))) for s in samples}
    em_mat = _matrix(feats, samples, em_v)
    uniq_mat = _matrix(feats, samples, {s: crit[s]["unique"] for s in samples}, None)
    key = [groups[f] for f in feats]
    agg = em_mat.drop(columns="feature_id").groupby(pd.Series(key, name="feature_id"), sort=False).sum().reset_index()
    agg_u = uniq_mat.drop(columns="feature_id").groupby(pd.Series(key, name="feature_id"), sort=False).sum().reset_index()
    boot_totals = np.array([sum(boot[s][b] for s in samples) for b in range(bootstraps)]) if samples \
        else np.zeros((bootstraps, len(feats)))

    # ---- publish: unpublish first, stage everything, move into place, seal last
    outdir.mkdir(parents=True, exist_ok=True)
    for sub in ("bootstrap", "realizations", "figures"):
        if (outdir / sub).is_symlink():
            raise ValueError(f"{outdir / sub} is a link; {STEP} writes its tables inside the output folder itself")
    previous = manifest.previous_artifacts(outdir)
    manifest.begin(outdir, STEP)
    stage = make_staging_dir(outdir, STEP)
    try:
        shutil.copyfile(quant_dir / "unique_counts.tsv", stage / "unique_counts.tsv")
        if manifest.sha256(stage / "unique_counts.tsv") != qrec["artifacts"]["unique_counts.tsv"]:
            raise manifest.ManifestError(f"{quant_dir / 'unique_counts.tsv'} changed after it was verified; "
                                         "run panel-quant again")
        _write(em_mat, stage / "em_expected_counts.tsv")
        _write(_matrix(feats, samples, norm), stage / "normalized_abundance.tsv")
        _write(agg, stage / "aggregated_counts.tsv")
        _write(agg_u, stage / "aggregated_unique_counts.tsv", None)
        _write(pd.DataFrame({"feature_id": feats, "group": key}), stage / "aggregated_groups.tsv", None)
        _write(abundance, stage / "em_bootstrap_abundance.tsv")
        _write(assignment, stage / "em_assignment_uncertainty.tsv")
        _write(support, stage / "viral_unique_read_support.tsv",
               header=[confidence_rule(unique_fraction_cut),
                       f"em_driven: unique_fraction < {unique_fraction_cut:g}; criteria_ratio: any_candidate_total / "
                       "max(unique_total, 1); criteria_min/max over unique, even split, EM and any-candidate totals; "
                       f"CV medians over samples whose mean is >= {min_mean_for_cv:g} fragments"])
        names = list(CORE_FILES)
        for b in range(bootstraps):
            names.append(bootstrap_name(b + 1))
            _write(_matrix(feats, samples, {s: boot[s][b] for s in samples}), stage / names[-1])
        for r in range(realizations):
            names.append(realization_name(r + 1))
            _write(_matrix(feats, samples, {s: real[s][r] for s in samples}, None), stage / names[-1], None)
        _figures(support, boot_totals, feats, stage)
        names += list(FIGURES)
        manifest.remove_stale(outdir, previous, names)
        for sub in ("bootstrap", "realizations", "figures"):
            (outdir / sub).mkdir(exist_ok=True)
        for n in names:
            os.replace(stage / n, outdir / n)
    finally:
        remove_staging_dir(stage)
    inputs_sha = {n: qrec["artifacts"][n] for n in ("features.tsv", *panel.CLASS_ARTIFACTS)}
    inputs_sha["provenance.json"] = manifest.sha256(quant_dir / "provenance.json")
    inputs_sha.update(tax_sha)
    import matplotlib
    manifest.seal(outdir, names, {
        "step": STEP, "quant": str(quant_dir.resolve()), "quant_provenance_sha256": inputs_sha["provenance.json"],
        "panel_id": (qrec.get("panel") or {}).get("panel_id"), "samples": samples,
        "bootstraps": bootstraps, "realizations": realizations, "seed": seed, "rng": RNG_DESCRIPTION,
        "aggregate_by": agg_desc, "min_mean_for_cv": min_mean_for_cv, "unique_fraction_cut": unique_fraction_cut,
        "confidence_rule": confidence_rule(unique_fraction_cut), "bootstrap_not_converged": int(not_conv),
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                     "matplotlib": matplotlib.__version__, "unclassified_meta": __version__},
        "inputs_sha256": inputs_sha})
    counts = support["assignment_confidence"].value_counts()
    log.info("%s: %d features x %d samples, %d bootstraps, %d realizations; confidence high %d, moderate %d, low %d "
             "-> %s", STEP, len(feats), len(samples), bootstraps, realizations, counts.get("high", 0),
             counts.get("moderate", 0), counts.get("low", 0), outdir)
    return outdir


def verify_uncertainty(folder) -> dict:
    """The manifest record of a complete em-uncertainty folder whose files are unchanged (ManifestError otherwise).
    The record lists `bootstraps` and `realizations`; the tables are bootstrap/counts_NNN.tsv and
    realizations/counts_NNN.tsv for NNN = 001..that number (layout in the module notes)."""
    folder = Path(folder)
    rec = manifest.verify(folder, required=CORE_FILES, what="em-uncertainty tables")
    if rec.get("step") != STEP:
        raise manifest.ManifestError(f"{folder / manifest.NAME} was not written by {STEP} (step {rec.get('step')!r})")
    want = [bootstrap_name(b + 1) for b in range(int(rec.get("bootstraps", 0)))] + \
        [realization_name(r + 1) for r in range(int(rec.get("realizations", 0)))]
    missing = [n for n in want if n not in rec["artifacts"]]
    if missing:
        raise manifest.ManifestError(f"{folder / manifest.NAME} does not cover {', '.join(missing[:3])}; run {STEP} again")
    return rec


# ---------------------------------------------------------------- CLI

def add_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument("--quant", required=True, help="sealed panel-quant folder (with equivalence_classes.tsv)")
    p.add_argument("--bootstraps", type=int, default=100, help="read-bootstrap replicates per sample (default 100)")
    p.add_argument("--realizations", type=int, default=20,
                   help="integer posterior realizations of the assignment (default 20)")
    p.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
    p.add_argument("--aggregate-by", default=None, metavar="COLUMN",
                   help="features.tsv column to sum EM counts by (default parent_feature when filled, else label)")
    p.add_argument("--taxonomy", default=None, metavar="panel_taxonomy.tsv",
                   help="sum by the taxID of a sealed `panel-taxonomy` table instead of --aggregate-by")
    p.add_argument("--min-mean-for-cv", type=float, default=1.0,
                   help="CV is reported only when the mean is at least this many fragments (default 1)")
    p.add_argument("--unique-fraction-cut", type=float, default=0.5,
                   help="em_driven when unique / EM fragments is below this (default 0.5)")
    p.add_argument("-t", "--threads", type=int, default=1, help="processes, one sample each (default 1)")
    p.add_argument("-o", "--outdir", required=True, help="output folder")


def main(args) -> int:
    if args.taxonomy and args.aggregate_by:
        raise ValueError("use either --aggregate-by or --taxonomy, not both")
    run(args.quant, args.outdir, bootstraps=args.bootstraps, realizations=args.realizations, seed=args.seed,
        aggregate_by=args.aggregate_by, taxonomy=args.taxonomy, threads=args.threads,
        min_mean_for_cv=args.min_mean_for_cv, unique_fraction_cut=args.unique_fraction_cut)
    return 0
