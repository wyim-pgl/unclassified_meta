"""da-robustness and da-resolutions: how far a differential-abundance call survives other methods, resampled read
tables, other prevalence filters, dropped samples and other clustering resolutions (docs/DESIGN_robustness.md, 3).

Statistics and limits
---------------------
* Methods. ANCOM-BC2 (r/ancombc2.R; bias-corrected natural-log fold change), ALDEx2 (r/aldex2.R; aldex.glm
  coefficient on the log2 CLR scale, mean over Dirichlet Monte-Carlo instances; counts rounded) and LinDA
  (r/linda.R; bias-corrected log2 CLR coefficient, fixed or mixed model). Every method sees the same samples,
  formula, reference level and the SAME features: one prevalence filter (--prv-cut, share of samples with a count
  > 0 in the primary table) is applied here, and the methods' own filters are off. Effects are on different
  scales; they are compared by direction and rank only, never averaged.
* Analyses. primary (EM expected counts, or --counts), unique-only, bootstrap tables (reads resampled within each
  sample, EM refitted), posterior realizations (integer draws of the read assignment), leave-one-sample-out (or
  leave-one-unit-out with --unit) and alternative prevalence cuts. In every analysis table a feature without a
  single read cannot be tested by any method and is left out of that fit for all methods alike.
* Resampling is not replication. Bootstrap and realization fits measure how much the call depends on which reads
  were sampled and where ambiguous reads were assigned. They are not biological replicates: they are summarized
  as the share of fits that reproduce the primary direction / significance, and their p-values are NEVER pooled
  (Rubin's rules stay with `da --imputations`, the one justified pooling, for ANCOM-BC2).
* Stability arithmetic. For a method, the denominator is the number of SUCCESSFUL fits of that kind; a feature not
  tested in a successful fit (no reads in that table, or dropped by the method) counts as neither same-direction
  nor significant. A failed fit is listed in fit_status.tsv and left out of every denominator.
* Agreement is not validation. The methods share the reads, the samples and the model; their agreement says the
  call does not hinge on one framework's assumptions. It is reported as "robust across three complementary
  differential abundance frameworks", never as validation by independent methods.
* Dimensions stay separate: assignment_confidence (em-uncertainty), abundance_stability (bootstrap CV),
  da_support, cross_method_agreement, replicate_consistency and resolution_stability are separate columns; the
  only combination is the printed rule behind `conclusion`.
* Power. With few samples per group a non-significant result is not evidence of no effect; mde_80 = (1.96 +
  0.8416) x SE is the effect a two-sided 5 % test detects with 80 % power, on each method's own scale.
"""
from __future__ import annotations

import logging
import os
import platform
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import da, manifest
from .runner import StepError

log = logging.getLogger("unclassified_meta")

METHODS = ("ancombc2", "aldex2", "linda")
R_DIR = Path(__file__).parent / "r"
SCRIPTS = {"ancombc2": R_DIR / "ancombc2.R", "aldex2": R_DIR / "aldex2.R", "linda": R_DIR / "linda.R"}
R_PACKAGES = {"ancombc2": "ANCOMBC", "aldex2": "ALDEx2", "linda": "MicrobiomeStat"}
METHOD_NAMES = {"ancombc2": "ANCOM-BC2", "aldex2": "ALDEx2", "linda": "LinDA"}
SCALE = {"ancombc2": "ln", "aldex2": "log2 CLR", "linda": "log2 CLR"}
ANCOMBC2_DEFINITION = "ANCOM-BC2 bias-corrected log fold change (natural log)"
STANDARD_COLUMNS = ("feature_id", "term", "effect", "se", "ci_low", "ci_high", "stat", "p", "q", "effect_definition")
Z_MDE80 = 1.96 + 0.8416                   # two-sided alpha 0.05, power 0.80
DIRECTION_MIN, SIGNIFICANCE_MIN, LOO_DIRECTION_MIN = 0.9, 0.8, 0.9
CV_STABLE, CV_MODERATE = 0.2, 0.5
CONCORDANCE_CLASSES = ("all_significant_same_direction", "majority_significant_same_direction", "single_method",
                       "conflicting_direction", "none")
NUMBER_WORDS = {1: "one", 2: "two", 3: "three", 4: "four"}

RESAMPLING_NOTE = ("bootstrap tables and posterior realizations resample the reads / read assignments of the SAME "
                   "samples: they measure read-sampling and assignment instability, not biological variation, are "
                   "never counted as replicates, and their p-values are not pooled")
NOT_TESTED_NOTE = ("a feature not tested in a successful fit (no reads in that table, or dropped by the method) "
                   "counts as not-same-direction and not-significant; failed fits (fit_status.tsv) are excluded from "
                   "the denominators")
AGREEMENT_NOTE = ("agreement between the methods is not independent validation: they analyse the same reads, "
                  "samples and model, so they share every upstream error (read assignment, contamination, design); "
                  "agreement shows only that the call does not hinge on one framework's assumptions")
MI_NOTE = "resampled fits summarized by shares, p-values not pooled (Rubin's rules: da --imputations, ANCOM-BC2)"
CONCLUSION_RULE = (
    "conclusion = high_confidence only when ALL hold: (1) q < alpha at the primary analysis in at least 2 methods, "
    "all with the same effect direction; (2) in each of those methods the bootstrap AND realization direction "
    f"stability >= {DIRECTION_MIN} and significance stability >= {SIGNIFICANCE_MIN} (not assessed = not met); "
    "(3) assignment_confidence is not low; (4) abundance_stability is not unstable; (5) replicate_consistency is "
    "true. not_significant: q >= alpha in every method at the primary analysis. Otherwise exploratory.")


def preferred_wording(n_methods: int = 3) -> str:
    return (f"Robust across {NUMBER_WORDS.get(n_methods, str(n_methods))} complementary differential abundance "
            "frameworks")


# ------------------------------------------------------------------------------------------------ table I/O

def _read_tsv(path, str_cols=("feature_id",)) -> pd.DataFrame:
    """A TSV written by this program: leading '#' lines are notes; a '#' inside a row is data; identifiers stay
    text ("01" stays "01")."""
    skip = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            skip += 1
    head = pd.read_csv(path, sep="\t", skiprows=skip, nrows=0, quoting=3).columns
    return pd.read_csv(path, sep="\t", skiprows=skip, quoting=3,
                       dtype={c: str for c in str_cols if c in head}, keep_default_na=False,
                       na_values=["", "NA", "NaN", "nan"])


def read_counts(path) -> pd.DataFrame:
    """feature x sample table (first column feature_id) as floats indexed by feature_id."""
    df = _read_tsv(path)
    if not len(df.columns) or df.columns[0] != "feature_id":
        raise StepError(f"{path}: the first column must be 'feature_id'")
    df = df.set_index("feature_id")
    if df.index.duplicated().any():
        raise StepError(f"{path}: duplicate feature_id {df.index[df.index.duplicated()][0]}")
    df.columns = [str(c) for c in df.columns]
    vals = df.apply(pd.to_numeric, errors="coerce")
    if vals.isna().any().any() or (vals < 0).any().any():
        raise StepError(f"{path}: counts must be finite, non-negative numbers")
    return vals.astype(float)


def write_tsv(path, df: pd.DataFrame, header=()) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fo:
        for h in header:
            fo.write("# " + " ".join(str(h).split()) + "\n")
        df.to_csv(fo, sep="\t", index=False, na_rep="NA", float_format="%.6g")


def read_standard(path) -> pd.DataFrame:
    df = _read_tsv(path, ("feature_id", "term", "effect_definition"))
    missing = [c for c in STANDARD_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path} lacks column(s) {', '.join(missing)}")
    for c in ("effect", "se", "ci_low", "ci_high", "stat", "p", "q"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df[list(STANDARD_COLUMNS)]


def ancombc2_to_standard(res: pd.DataFrame) -> pd.DataFrame:
    """ANCOM-BC2's wide table (lfc_/se_/W_/p_/q_<term>) as the standard long table; the intercept is left out.
    stat = W (lfc / se); ci = lfc +- 1.96 se (ANCOM-BC2's own test uses the same normal reference)."""
    terms = [c[4:] for c in res.columns if c.startswith("lfc_") and c != "lfc_(Intercept)"]
    parts = []
    num = lambda c: pd.to_numeric(res[c], errors="coerce") if c in res.columns else pd.Series(np.nan, index=res.index)
    for t in terms:
        eff, se = num(f"lfc_{t}"), num(f"se_{t}")
        parts.append(pd.DataFrame({"feature_id": res["feature_id"].astype(str), "term": t, "effect": eff, "se": se,
                                   "ci_low": eff - 1.959964 * se, "ci_high": eff + 1.959964 * se,
                                   "stat": num(f"W_{t}"), "p": num(f"p_{t}"), "q": num(f"q_{t}"),
                                   "effect_definition": ANCOMBC2_DEFINITION}))
    if not parts:
        return pd.DataFrame(columns=list(STANDARD_COLUMNS))
    return pd.concat(parts, ignore_index=True)[list(STANDARD_COLUMNS)]


# ------------------------------------------------------------------------------------------------ inputs

@dataclass
class Inputs:
    primary: str
    features: str
    metadata: str
    unique: str | None = None
    bootstraps: list = field(default_factory=list)
    realizations: list = field(default_factory=list)
    support: str | None = None
    boot_abundance: str | None = None
    normalized: str | None = None
    resolution: str | None = None
    expected: dict = field(default_factory=dict)        # path -> sha256 vouched for by a manifest (None: unsealed)
    em_record: dict = field(default_factory=dict)

    def paths(self) -> list[str]:
        return [p for p in (self.primary, self.features, self.metadata, self.unique, self.support,
                            self.boot_abundance, self.normalized, self.resolution, *self.bootstraps,
                            *self.realizations) if p]


def _sealed_sha(path, unsealed_ok: bool, what: str) -> str | None:
    """sha256 recorded for `path` by its folder's manifest; None for a hand-made table (only with --unsealed-input).
    The same rule as da.verify_inputs."""
    folder = Path(path).parent
    if any((folder / n).exists() for n in ("provenance.json", manifest.NAME)):
        return manifest.verify_file(path, what=what)["artifacts"][Path(path).name]
    if manifest.is_generated(folder) or (folder / "ledger.tsv").exists():
        raise StepError(f"{folder} was written by this program but has no complete manifest; run that step again")
    if not unsealed_ok:
        raise StepError(f"{path} has no manifest; pass --unsealed-input to use a table made by hand as it is")
    log.warning("using %s without a manifest (--unsealed-input): no completeness or checksum check", path)
    return None


def verify_em_folder(folder) -> dict:
    """The sealed record of an em-uncertainty folder (uncertainty.verify_uncertainty when that module provides it,
    else the generic manifest check)."""
    rec = None
    try:
        from . import uncertainty                                  # written by the em-uncertainty step
        fn = getattr(uncertainty, "verify_uncertainty", None)
    except ImportError:
        fn = None
    if fn is not None:
        rec = fn(folder)
    if not isinstance(rec, dict) or not rec.get("artifacts"):
        rec = manifest.verify(folder, what="em-uncertainty tables")
    if rec.get("step") not in (None, "em-uncertainty"):
        raise StepError(f"{folder} is a '{rec.get('step')}' folder, not an em-uncertainty folder")
    if "em_expected_counts.tsv" not in rec["artifacts"]:
        raise StepError(f"{folder}: the manifest does not seal em_expected_counts.tsv")
    return rec


def _numbered(arts: dict, sub: str) -> list[str]:
    pat = re.compile(rf"{sub}/counts_(\d+)\.tsv")
    found = [(int(m.group(1)), k) for k in arts if (m := pat.fullmatch(k))]
    return [k for _, k in sorted(found)]


def collect_inputs(args) -> Inputs:
    unsealed = bool(getattr(args, "unsealed_input", False))
    if bool(args.em_uncertainty) == bool(args.counts):
        raise StepError("give --em-uncertainty DIR (primary = EM expected counts) or --counts FILE, not both")
    if args.em_uncertainty:
        folder = Path(args.em_uncertainty)
        rec = verify_em_folder(folder)
        arts = rec["artifacts"]
        opt = lambda n: str(folder / n) if n in arts else None
        boots, reals = _numbered(arts, "bootstrap"), _numbered(arts, "realizations")
        nb, nr = min(args.bootstrap_fits, len(boots)), min(args.realization_fits, len(reals))
        if nb < args.bootstrap_fits:
            log.warning("da-robustness: %d bootstrap table(s) available, %d requested", len(boots), args.bootstrap_fits)
        if nr < args.realization_fits:
            log.warning("da-robustness: %d realization table(s) available, %d requested", len(reals),
                        args.realization_fits)
        inp = Inputs(primary=str(folder / "em_expected_counts.tsv"), features=str(args.features),
                     metadata=str(args.metadata), unique=opt("unique_counts.tsv"),
                     bootstraps=[str(folder / k) for k in boots[:nb]], realizations=[str(folder / k) for k in reals[:nr]],
                     support=opt("viral_unique_read_support.tsv"), boot_abundance=opt("em_bootstrap_abundance.tsv"),
                     normalized=opt("normalized_abundance.tsv"), em_record=rec)
        for p in inp.paths():
            if Path(p).parent == folder or folder in Path(p).parents:
                inp.expected[p] = arts[str(Path(p).relative_to(folder))]
    else:
        inp = Inputs(primary=str(args.counts), features=str(args.features), metadata=str(args.metadata))
        inp.expected[inp.primary] = _sealed_sha(inp.primary, unsealed, "count tables")
        if args.bootstrap_fits or args.realization_fits:
            log.info("da-robustness: --counts has no bootstrap or realization tables; give --em-uncertainty for "
                     "the read-resampling analyses")
    inp.expected[inp.features] = _sealed_sha(inp.features, unsealed, "feature tables")
    inp.expected[inp.metadata] = None
    if getattr(args, "resolution_table", None):
        inp.resolution = str(args.resolution_table)
        inp.expected[inp.resolution] = _sealed_sha(inp.resolution, True, "resolution tables")
    return inp


# ------------------------------------------------------------------------------------------------ design

def prevalence_filter(counts: pd.DataFrame, cut: float) -> list[str]:
    """Features present (count > 0) in at least `cut` of the samples - applied ONCE, then given to every method."""
    prev = (counts > 0).mean(axis=1)
    return list(counts.index[prev >= cut - 1e-12])


def analysis_table(table: pd.DataFrame, features: list, samples: list, what: str) -> pd.DataFrame:
    """The analysis' table on the shared feature set; a feature without a read in it cannot be tested by any
    method and is left out for all of them."""
    miss_f = [f for f in features if f not in table.index]
    miss_s = [s for s in samples if s not in table.columns]
    if miss_f or miss_s:
        raise StepError(f"{what}: missing {'features ' + ', '.join(miss_f[:3]) if miss_f else ''}"
                        f"{' samples ' + ', '.join(miss_s[:3]) if miss_s else ''} of the primary table")
    t = table.loc[features, samples]
    return t[(t > 0).any(axis=1)]


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(name)) or "_"


@dataclass(frozen=True)
class Analysis:
    name: str
    kind: str                     # primary, unique, boot, real, loo, prv
    counts: Path
    metadata: Path
    left_out: str = ""
    n_features: int = 0


@dataclass(frozen=True)
class Job:
    method: str
    analysis: Analysis
    outdir: Path


def group_sizes(meta: pd.DataFrame, group: str) -> dict:
    g = meta[group].astype(str)
    g = g[~g.str.strip().isin(["", "NA"])]
    return {k: int(v) for k, v in g.value_counts().sort_index().items()}


def _command(rscript: str, job: Job, features: Path, args) -> list[str]:
    m, a = job.method, job.analysis
    cmd = [rscript, str(SCRIPTS[m]), "--counts", str(a.counts), "--metadata", str(a.metadata), "--group", args.group,
           "--outdir", str(job.outdir), "--alpha", str(args.alpha), "--p-adjust", args.p_adjust, "--threads", "1"]
    if m == "ancombc2":           # the shared filter is done; ANCOM-BC2's own prevalence / library filters are off
        cmd += ["--features", str(features), "--prv-cut", "0", "--lib-cut", "0"]
    else:
        cmd += ["--seed", str(args.seed)]
    if m == "aldex2":
        cmd += ["--mc-samples", str(args.mc_samples)]
    if args.formula:
        cmd += ["--formula", args.formula]
    if args.reference_level:
        cmd += ["--reference-level", args.reference_level]
    if args.rand_formula:         # aldex2.R notes that it ignores it (fixed effects only)
        cmd += ["--rand-formula", args.rand_formula]
    return cmd


def run_job(job: Job, rscript: str, features: Path, args, stage: Path, outdir: Path) -> tuple[dict, pd.DataFrame | None]:
    """One method on one analysis table. Returns its fit_status row and its standard table (None when it failed)."""
    job.outdir.mkdir(parents=True, exist_ok=True)
    status = {"method": job.method, "analysis": job.analysis.name, "kind": job.analysis.kind,
              "left_out": job.analysis.left_out, "features_in": job.analysis.n_features, "status": "failed",
              "rows": 0, "message": ""}
    try:
        p = subprocess.run(_command(rscript, job, features, args), capture_output=True, text=True)
    except OSError as e:
        status["message"] = f"could not start Rscript: {e}"
        return status, None
    if p.returncode != 0:
        lines = [x.strip() for x in (p.stderr or "").splitlines() if x.strip()]
        status["message"] = (f"exit {p.returncode}: " + " | ".join(lines[-3:])).replace(str(stage), str(outdir))
        log.warning("da-robustness: %s on %s failed (%s)", job.method, job.analysis.name, status["message"])
        return status, None
    try:
        if job.method == "ancombc2":
            res = ancombc2_to_standard(da.read_results(job.outdir / "ancombc2_results.tsv"))
        else:
            res = read_standard(job.outdir / "results.tsv")
    except (OSError, ValueError, KeyError) as e:
        status["message"] = f"unreadable result: {e}".replace(str(stage), str(outdir))
        return status, None
    write_tsv(job.outdir / "results.tsv", res)
    if not (job.outdir / "warnings.txt").exists():
        (job.outdir / "warnings.txt").write_text("none\n")
    status.update(status="ok", rows=len(res))
    return status, res


def missing_r_packages(rscript: str, methods) -> list[str]:
    pk = ", ".join(f'"{R_PACKAGES[m]}"' for m in methods)
    expr = f"x <- c({pk}); cat(x[!vapply(x, requireNamespace, logical(1), quietly = TRUE)], sep = '\\n')"
    p = subprocess.run([rscript, "-e", expr], capture_output=True, text=True)
    if p.returncode != 0:
        raise StepError(f"{rscript} could not check the R packages: {p.stderr.strip()[-300:]}")
    return [x for x in p.stdout.split() if x]


# ------------------------------------------------------------------------------------------------ statistics

def _index(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop_duplicates(["feature_id", "term"]).set_index(["feature_id", "term"])


def resample_stability(primary: pd.DataFrame, fits: list, alpha: float) -> pd.DataFrame:
    """Per (feature, term) of `primary`: n successful fits, median and 2.5/97.5 % of their effects, direction
    stability (share of fits whose effect has the primary's sign) and significance stability (share with
    q < alpha). Not tested in a fit = neither. `fits` holds the standard tables of the SUCCESSFUL fits only."""
    idx = primary.index
    n = len(fits)
    out = pd.DataFrame(index=idx)
    out["n"] = n
    if n == 0:
        for c in ("effect_median", "effect_p2_5", "effect_p97_5", "direction_stability", "significance_stability"):
            out[c] = np.nan
        return out
    eff = pd.concat([_index(f)["effect"].rename(i) for i, f in enumerate(fits)], axis=1).reindex(idx)
    q = pd.concat([_index(f)["q"].rename(i) for i, f in enumerate(fits)], axis=1).reindex(idx)
    s0 = np.sign(primary["effect"])
    same = np.sign(eff).eq(s0, axis=0) & eff.notna() & (s0 != 0).to_numpy()[:, None]
    sig = q.lt(alpha) & q.notna()
    known = primary["effect"].notna()
    out["effect_median"] = eff.median(axis=1)
    out["effect_p2_5"] = eff.quantile(0.025, axis=1)
    out["effect_p97_5"] = eff.quantile(0.975, axis=1)
    out["direction_stability"] = (same.sum(axis=1) / n).where(known)
    out["significance_stability"] = sig.sum(axis=1) / n
    return out


def classify_concordance(signs: dict, significant: set, n_methods: int) -> str:
    """signs: {method: sign of its primary effect} of the methods that tested the feature; significant: methods
    with q < alpha. all_* = significant in every requested method."""
    if not significant:
        return "none"
    if len({signs[m] for m in significant}) > 1:
        return "conflicting_direction"
    if len(significant) == 1:
        return "single_method"
    if len(significant) >= n_methods:
        return "all_significant_same_direction"
    return "majority_significant_same_direction"


def spearman(a: pd.Series, b: pd.Series) -> tuple[float, int]:
    both = a.notna() & b.notna()
    n = int(both.sum())
    if n < 3:
        return float("nan"), n
    ra, rb = a[both].rank().to_numpy(), b[both].rank().to_numpy()
    if ra.std() == 0 or rb.std() == 0:
        return float("nan"), n
    return float(np.corrcoef(ra, rb)[0, 1]), n


def concordance_table(primary: dict, methods: list, alpha: float, labels: dict) -> tuple[pd.DataFrame, list]:
    """cross_method_concordance.tsv rows and the pairwise Spearman header lines. primary: {method: standard table
    of its primary fit} (methods whose primary fit failed are absent: not tested)."""
    keys = list(dict.fromkeys(k for m in methods if m in primary for k in _index(primary[m]).index))
    idx = pd.MultiIndex.from_tuples(keys, names=["feature_id", "term"]) if keys else \
        pd.MultiIndex.from_tuples([], names=["feature_id", "term"])
    eff = {m: (_index(primary[m])["effect"].reindex(idx) if m in primary else pd.Series(np.nan, index=idx))
           for m in methods}
    qv = {m: (_index(primary[m])["q"].reindex(idx) if m in primary else pd.Series(np.nan, index=idx))
          for m in methods}
    rows = []
    for k in idx:
        signs = {m: float(np.sign(eff[m][k])) for m in methods if np.isfinite(eff[m][k])}
        sig = {m for m in signs if np.isfinite(qv[m][k]) and qv[m][k] < alpha}
        row = {"feature_id": k[0], "label": labels.get(k[0], k[0]), "term": k[1], "n_methods_tested": len(signs),
               "n_methods_significant": len(sig), "methods_significant": ",".join(m for m in methods if m in sig)}
        for m in methods:
            row[f"sign_{m}"] = signs.get(m, np.nan)
            row[f"effect_{m}"] = eff[m][k]
            row[f"q_{m}"] = qv[m][k]
        row["direction_agreement"] = bool(signs) and len(set(signs.values())) == 1 and 0.0 not in signs.values()
        row["concordance"] = classify_concordance(signs, sig, len(methods))
        rows.append(row)
    cols = ["feature_id", "label", "term", "n_methods_tested", "n_methods_significant", "methods_significant",
            *[f"{p}_{m}" for m in methods for p in ("sign", "effect", "q")], "direction_agreement", "concordance"]
    table = pd.DataFrame(rows, columns=cols)
    header = []
    for i, a in enumerate(methods):
        for b in methods[i + 1:]:
            r, n = spearman(eff[a], eff[b])
            header.append(f"Spearman rank correlation of primary effects, {a} vs {b}: "
                          f"{'NA' if np.isnan(r) else f'{r:.3f}'} over {n} feature x term pairs tested by both")
    return table, header


def max_single_sample_share(abund: pd.DataFrame, groups: pd.Series) -> pd.Series:
    """Per feature: the largest share of one group's total (normalized) abundance contributed by a single sample
    (1.0 = everything in that group comes from one sample)."""
    out = pd.Series(np.nan, index=abund.index)
    for g in groups.unique():
        cols = [s for s in groups.index[groups == g] if s in abund.columns]
        if not cols:
            continue
        sub = abund[cols]
        tot = sub.sum(axis=1)
        share = (sub.max(axis=1) / tot).where(tot > 0)
        out = pd.concat([out, share], axis=1).max(axis=1)
    return out


def replicate_rows(method: str, primary: pd.DataFrame, loo: list, alpha: float, labels: dict) -> pd.DataFrame:
    """Long leave-one-out rows plus one SUMMARY row per (feature, term). loo: [(left_out, standard table)] of the
    successful leave-one-out fits."""
    p = _index(primary)
    s0 = np.sign(p["effect"])
    psig = p["q"].lt(alpha)
    long, same_cols, sig_cols = [], [], []
    for left, res in loo:
        r = _index(res).reindex(p.index)
        same = (np.sign(r["effect"]) == s0) & r["effect"].notna() & (s0 != 0)
        sig = r["q"].lt(alpha) & r["q"].notna()
        same_cols.append(same.rename(left))
        sig_cols.append(sig.rename(left))
        long.append(pd.DataFrame({"method": method, "feature_id": p.index.get_level_values(0),
                                  "term": p.index.get_level_values(1), "left_out": left,
                                  "effect": r["effect"].to_numpy(), "q": r["q"].to_numpy(),
                                  "significant": sig.to_numpy(), "same_direction": same.to_numpy()}))
    n = len(loo)
    summ = pd.DataFrame({"method": method, "feature_id": p.index.get_level_values(0),
                         "term": p.index.get_level_values(1), "left_out": "SUMMARY",
                         "effect": p["effect"].to_numpy(), "q": p["q"].to_numpy(), "significant": psig.to_numpy(),
                         "same_direction": np.nan, "loo_n": n})
    if n:
        same_m, sig_m = pd.concat(same_cols, axis=1), pd.concat(sig_cols, axis=1)
        summ["loo_direction_stability"] = (same_m.sum(axis=1) / n).where(p["effect"].notna()).to_numpy()
        summ["loo_significance_stability"] = (sig_m.sum(axis=1) / n).to_numpy()
        drivers = [",".join(c for c in sig_m.columns if not sig_m.loc[k, c]) if psig[k] else "" for k in p.index]
        summ["single_sample_driven"] = [bool(psig[k]) and bool(d) for k, d in zip(p.index, drivers)]
        summ["drivers"] = drivers
    else:
        summ["loo_direction_stability"] = np.nan
        summ["loo_significance_stability"] = np.nan
        summ["single_sample_driven"] = False
        summ["drivers"] = ""
    summ["mde_80"] = (Z_MDE80 * p["se"]).to_numpy()
    out = pd.concat([*long, summ], ignore_index=True)
    out.insert(2, "label", out["feature_id"].map(lambda f: labels.get(f, f)))
    return out


def classify_cv(cv: float) -> str:
    if cv is None or not np.isfinite(cv):
        return "not_assessed"
    return "stable" if cv <= CV_STABLE else ("moderate" if cv <= CV_MODERATE else "unstable")


def conclude(sig_signs: dict, stability: dict, assignment: str, abundance: str,
             replicate_consistent: bool) -> tuple[str, list]:
    """The printed rule (CONCLUSION_RULE). sig_signs: {method: sign} of the methods with q < alpha at the primary
    analysis; stability: {method: {boot_direction_stability, real_direction_stability,
    boot_significance_stability, real_significance_stability}}. Returns (conclusion, unmet conditions)."""
    if not sig_signs:
        return "not_significant", ["no method has q < alpha at the primary analysis"]
    why = []
    if len(sig_signs) < 2:
        why.append("significant in fewer than 2 methods")
    if len(set(sig_signs.values())) > 1:
        why.append("the significant methods disagree in direction")
    for m in sig_signs:
        st = stability.get(m, {})
        for key, lim in (("boot_direction_stability", DIRECTION_MIN), ("real_direction_stability", DIRECTION_MIN),
                         ("boot_significance_stability", SIGNIFICANCE_MIN),
                         ("real_significance_stability", SIGNIFICANCE_MIN)):
            v = st.get(key, np.nan)
            if v is None or not np.isfinite(v):
                why.append(f"{m} {key} not assessed")
            elif v < lim:
                why.append(f"{m} {key} {v:.2f} < {lim}")
    if assignment == "low":
        why.append("assignment_confidence low")
    if abundance == "unstable":
        why.append("abundance_stability unstable")
    if not replicate_consistent:
        why.append("replicate_consistency false")
    return ("high_confidence" if not why else "exploratory"), why


def _fmt(x, nd=3) -> str:
    try:
        x = float(x)
    except (TypeError, ValueError):
        return str(x)
    return f"{x:.{nd}g}" if np.isfinite(x) else "NA"


# ------------------------------------------------------------------------------------------------ CLI

def add_arguments_robustness(p) -> None:
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--em-uncertainty", metavar="DIR",
                     help="sealed em-uncertainty folder: primary = em_expected_counts.tsv, plus unique-only, "
                          "bootstrap and realization tables")
    src.add_argument("--counts", help="a single count table (primary only; no read-resampling analyses)")
    p.add_argument("--features", required=True, help="features.tsv of the panel (feature_id, label, parent_feature)")
    p.add_argument("--metadata", required=True, help="TSV with a 'sample' column, the group column and covariates")
    p.add_argument("--group", required=True, help="metadata column that holds the groups to compare")
    p.add_argument("--formula", help="fixed-effects formula, e.g. 'group + batch' (default: the group column)")
    p.add_argument("--reference-level", help="group used as the baseline")
    p.add_argument("--rand-formula", metavar="'(1|unit)'",
                   help="random effects (ANCOM-BC2 and LinDA; ALDEx2 has none and fits the fixed effects only)")
    p.add_argument("--unit", help="metadata column (e.g. plant): leave out one unit at a time instead of one sample")
    p.add_argument("--methods", default=",".join(METHODS), help="comma list of ancombc2, aldex2, linda")
    p.add_argument("--prv-cut", type=float, default=0.10,
                   help="shared prevalence filter: keep features with a count > 0 in at least this share of samples "
                        "of the primary table (the methods' own filters are off)")
    p.add_argument("--prv-cuts", default="0.25,0.5", help="alternative prevalence cuts refitted on the primary table "
                                                          "('' for none)")
    p.add_argument("--bootstrap-fits", type=int, default=20, help="bootstrap tables refitted (<= available)")
    p.add_argument("--realization-fits", type=int, default=10, help="posterior realizations refitted (<= available)")
    p.add_argument("--resolution-table", help="multi_resolution_robustness.tsv from da-resolutions")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--p-adjust", default="holm", help="multiple-testing adjustment of every method (per term)")
    p.add_argument("--mc-samples", type=int, default=128, help="ALDEx2 Monte-Carlo Dirichlet instances")
    p.add_argument("--seed", type=int, default=1, help="seed of ALDEx2 and LinDA (ANCOM-BC2 uses its fixed seed)")
    p.add_argument("-t", "--threads", type=int, default=1, help="fits run in parallel (one R process each)")
    p.add_argument("-o", "--outdir", default="da_robustness")
    p.add_argument("--rscript", help="Rscript of the R environment (default: $UNCLASSIFIED_META_RSCRIPT, then PATH)")
    p.add_argument("--unsealed-input", action="store_true",
                   help="use tables that carry no manifest (made by hand) as they are")


def add_arguments_resolutions(p) -> None:
    p.add_argument("--run", action="append", required=True, metavar="LABEL=DIR",
                   help="a da-robustness output folder at one clustering resolution, e.g. 98=da_rob_98 (repeat)")
    p.add_argument("--parent-column", default="parent_feature",
                   help="column of each run's features.used.tsv that names the parent (feature_id when empty)")
    p.add_argument("-o", "--outdir", default="da_resolutions")


def _parse_methods(text: str) -> list[str]:
    ms = [m.strip().lower() for m in str(text).split(",") if m.strip()]
    bad = [m for m in ms if m not in METHODS]
    if bad or not ms:
        raise StepError(f"--methods: unknown {', '.join(bad) or '(none)'}; choose from {', '.join(METHODS)}")
    return list(dict.fromkeys(ms))


def _parse_cuts(text: str) -> list[float]:
    out = []
    for x in str(text or "").split(","):
        if x.strip():
            v = float(x)
            if not 0 <= v <= 1:
                raise StepError(f"--prv-cuts: {v} is not a share between 0 and 1")
            out.append(v)
    return list(dict.fromkeys(out))


def _guard(outdir: Path, inputs, command: str):
    from .io_utils import folder_lock, preflight
    preflight(command, inputs, owned_dirs=[outdir])
    for sub in ("fits", "figures"):
        if (outdir / sub).is_symlink():
            raise StepError(f"{outdir / sub} is a link; {command} writes its files inside the output folder itself")
    return folder_lock(outdir, command)


def main_robustness(args) -> int:
    """`da-robustness`: checks the output folder against the inputs, locks it, runs everything, seals it."""
    outdir = Path(args.outdir)
    ins = [x for x in (args.em_uncertainty, args.counts, args.features, args.metadata,
                       getattr(args, "resolution_table", None)) if x]
    with _guard(outdir, ins, "da-robustness"):
        run_robustness(args)
    return 0


def main_resolutions(args) -> int:
    outdir = Path(args.outdir)
    runs = _parse_runs(args.run)
    with _guard(outdir, [d for _, d in runs], "da-resolutions"):
        run_resolutions(args)
    return 0


# ------------------------------------------------------------------------------------------------ da-robustness

def run_robustness(args) -> Path:
    from .io_utils import make_staging_dir, remove_staging_dir
    methods = _parse_methods(args.methods)
    cuts = _parse_cuts(args.prv_cuts)
    if not 0 <= args.prv_cut <= 1:
        raise StepError("--prv-cut must be a share between 0 and 1")
    if args.threads < 1 or args.mc_samples < 1 or args.bootstrap_fits < 0 or args.realization_fits < 0:
        raise StepError("--threads and --mc-samples must be >= 1; --bootstrap-fits / --realization-fits >= 0")
    inp = collect_inputs(args)
    rscript = da.find_rscript(args.rscript)
    missing = missing_r_packages(rscript, methods)
    if missing:
        raise StepError(f"{rscript} lacks the R package(s) {', '.join(missing)}; install environment-r.yml or "
                        "drop the method from --methods")
    outdir = Path(args.outdir)
    previous = manifest.previous_artifacts(outdir)
    manifest.begin(outdir, "da-robustness")
    stage = make_staging_dir(outdir, "da-robustness")
    try:
        cp = manifest.private_copies({p: inp.expected.get(p) for p in inp.paths()}, stage / "inputs")
        made, record = _analyse(args, inp, cp, methods, cuts, rscript, stage, outdir)
        _publish(stage, outdir, made, previous, record, inp.paths(), "da-robustness")
    finally:
        remove_staging_dir(stage)
    return outdir / "robustness_summary.tsv"


def _analyse(args, inp: Inputs, cp: dict, methods: list, cuts: list, rscript: str, stage: Path, outdir: Path):
    alpha = float(args.alpha)
    work = stage / "work"
    work.mkdir()
    primary = read_counts(cp[inp.primary])
    feats = _read_tsv(cp[inp.features], ("feature_id", "label", "parent_feature"))
    if "feature_id" not in feats.columns:
        raise StepError(f"{inp.features} has no feature_id column")
    meta = pd.read_csv(cp[inp.metadata], sep="\t", dtype=str, keep_default_na=False)   # every column as text
    for col in ("sample", args.group, *([args.unit] if args.unit else [])):
        if col not in meta.columns:
            raise StepError(f"{inp.metadata} has no column '{col}'")
    if meta["sample"].duplicated().any():
        raise StepError(f"{inp.metadata}: duplicate sample names")
    no_counts = [s for s in meta["sample"] if s not in primary.columns]
    if no_counts:
        log.warning("da-robustness: metadata samples without counts left out: %s", ", ".join(no_counts))
    meta = meta[meta["sample"].isin(primary.columns)].reset_index(drop=True)
    samples = list(meta["sample"])
    if len(samples) < 4:
        raise StepError(f"{len(samples)} sample(s) with counts and metadata; differential abundance needs more")
    extra = [s for s in primary.columns if s not in set(samples)]
    if extra:
        log.warning("da-robustness: samples with counts but not in the metadata left out: %s", ", ".join(extra))
    sizes = group_sizes(meta, args.group)
    selected = prevalence_filter(primary[samples], args.prv_cut)
    if len(selected) < 2:
        raise StepError(f"{len(selected)} feature(s) pass --prv-cut {args.prv_cut}; nothing to test")
    log.info("da-robustness: %d of %d features pass the shared prevalence filter (%s)", len(selected), len(primary),
             args.prv_cut)
    fi = feats.drop_duplicates("feature_id").set_index("feature_id")
    unknown = [f for f in selected if f not in fi.index]
    if unknown:
        log.warning("da-robustness: %d tested feature(s) are not in %s (no label): %s", len(unknown), inp.features,
                    ", ".join(unknown[:5]))
    used = fi.reindex(selected).reset_index()
    if "label" not in used.columns:
        used["label"] = np.nan
    used["label"] = used["label"].where(used["label"].notna() & (used["label"].astype(str) != ""), used["feature_id"])
    write_tsv(stage / "features.used.tsv", used)
    labels = dict(zip(used["feature_id"], used["label"].astype(str)))
    meta_path = stage / "metadata.used.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)

    # ---- analysis tables (all on the shared feature set)
    analyses: list[Analysis] = []

    def add(name, kind, table, mpath=meta_path, smp=samples, left=""):
        t = analysis_table(table, selected, smp, name)
        f = work / f"{name}.tsv"
        t.reset_index().to_csv(f, sep="\t", index=False, float_format="%.6f")
        analyses.append(Analysis(name, kind, f, mpath, left, len(t)))

    add("primary", "primary", primary)
    if inp.unique:
        add("unique", "unique", read_counts(cp[inp.unique]))
    for i, b in enumerate(inp.bootstraps, 1):
        add(f"boot_{i:03d}", "boot", read_counts(cp[b]))
    for i, r in enumerate(inp.realizations, 1):
        add(f"real_{i:03d}", "real", read_counts(cp[r]))
    if args.unit:
        units = [u for u in dict.fromkeys(meta[args.unit]) if u.strip() not in ("", "NA")]
        drop = [(u, list(meta["sample"][meta[args.unit] == u])) for u in units]
    else:
        drop = [(s, [s]) for s in samples]
    seen = set()
    for left, gone in drop:
        name = f"loo_{_safe(left)}"
        while name in seen:
            name += "_"
        seen.add(name)
        keep = meta[~meta["sample"].isin(gone)]
        mp = work / f"{name}.metadata.tsv"
        keep.to_csv(mp, sep="\t", index=False)
        add(name, "loo", primary, mp, list(keep["sample"]), left)
    for c in cuts:
        if c <= args.prv_cut:
            log.warning("da-robustness: --prv-cuts %s is not above --prv-cut %s (same features as primary)", c,
                        args.prv_cut)
        keep_f = [f for f in prevalence_filter(primary[samples], c) if f in set(selected)]
        if len(keep_f) < 2:
            log.warning("da-robustness: prevalence cut %s leaves %d feature(s); not fitted", c, len(keep_f))
            continue
        t = primary.loc[keep_f, samples]
        f = work / f"prv_{c:g}.tsv"
        t.reset_index().to_csv(f, sep="\t", index=False, float_format="%.6f")
        analyses.append(Analysis(f"prv_{c:g}", "prv", f, meta_path, "", len(t)))

    # ---- fits
    jobs = [Job(m, a, stage / "fits" / f"{m}__{a.name}") for m in methods for a in analyses]
    log.info("da-robustness: %d fits (%d methods x %d analyses) on %d thread(s)", len(jobs), len(methods),
             len(analyses), args.threads)
    features_path = stage / "features.used.tsv"
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        done = list(ex.map(lambda j: run_job(j, rscript, features_path, args, stage, outdir), jobs))
    status = pd.DataFrame([s for s, _ in done])
    results = {(j.method, j.analysis.name): r for j, (s, r) in zip(jobs, done) if r is not None}
    terms_of = lambda r: r[r["term"].astype(str).str.startswith(args.group)]     # group terms (and interactions)
    results = {k: terms_of(v) for k, v in results.items()}
    by_kind = {}
    for a in analyses:
        by_kind.setdefault(a.kind, []).append(a)
    prim = {m: results[(m, "primary")] for m in methods if (m, "primary") in results}
    if not prim:
        raise StepError("every method failed on the primary table; see the warnings above")

    # ---- per-method stability
    stab_rows = []
    stability_lookup = {}
    for m in methods:
        if m not in prim:
            continue
        p = _index(prim[m])
        ok = lambda kind: [results[(m, a.name)] for a in by_kind.get(kind, []) if (m, a.name) in results]
        bs, rs = resample_stability(p, ok("boot"), alpha), resample_stability(p, ok("real"), alpha)
        un = _index(results[(m, "unique")]).reindex(p.index) if (m, "unique") in results else None
        prv_ok = [_index(results[(m, a.name)]).reindex(p.index) for a in by_kind.get("prv", [])
                  if (m, a.name) in results]
        if prv_ok:
            filt = pd.concat([(r["q"] < alpha).rename(i) for i, r in enumerate(prv_ok)], axis=1).all(axis=1)
        else:
            filt = pd.Series(np.nan, index=p.index, dtype=object)
        defn = p["effect_definition"]
        df = pd.DataFrame({
            "method": m, "feature_id": p.index.get_level_values(0),
            "label": [labels.get(f, f) for f in p.index.get_level_values(0)], "term": p.index.get_level_values(1),
            "effect_definition": defn.to_numpy(), "primary_effect": p["effect"].to_numpy(),
            "primary_se": p["se"].to_numpy(), "primary_q": p["q"].to_numpy(),
            "primary_significant": (p["q"] < alpha).to_numpy(),
            "boot_n": bs["n"].to_numpy(), "boot_effect_median": bs["effect_median"].to_numpy(),
            "boot_effect_p2_5": bs["effect_p2_5"].to_numpy(), "boot_effect_p97_5": bs["effect_p97_5"].to_numpy(),
            "boot_direction_stability": bs["direction_stability"].to_numpy(),
            "boot_significance_stability": bs["significance_stability"].to_numpy(),
            "real_n": rs["n"].to_numpy(), "real_direction_stability": rs["direction_stability"].to_numpy(),
            "real_significance_stability": rs["significance_stability"].to_numpy(),
            "unique_effect": un["effect"].to_numpy() if un is not None else np.nan,
            "unique_q": un["q"].to_numpy() if un is not None else np.nan,
            "filter_significant_all": filt.to_numpy(), "mi_note": MI_NOTE})
        stab_rows.append(df)
        for _, r in df.iterrows():
            stability_lookup[(m, r["feature_id"], r["term"])] = r
    stability = pd.concat(stab_rows, ignore_index=True)

    concord, rho = concordance_table(prim, methods, alpha, labels)

    # ---- replicate sensitivity
    groups = meta.set_index("sample")[args.group]
    abund = (read_counts(cp[inp.normalized]) if inp.normalized else primary / primary.sum(axis=0) * 1e6)
    abund = abund.reindex(index=selected, columns=[s for s in samples if s in abund.columns])
    share = max_single_sample_share(abund, groups)
    n_per_group = ";".join(f"{k}={v}" for k, v in sizes.items())
    rep_parts = []
    for m in methods:
        if m not in prim:
            continue
        loo = [(a.left_out, results[(m, a.name)]) for a in by_kind.get("loo", []) if (m, a.name) in results]
        rep_parts.append(replicate_rows(m, prim[m], loo, alpha, labels))
    replicate = pd.concat(rep_parts, ignore_index=True)
    replicate["n_per_group"] = np.where(replicate["left_out"] == "SUMMARY", n_per_group, "")
    replicate["max_single_sample_share"] = np.where(replicate["left_out"] == "SUMMARY",
                                                    replicate["feature_id"].map(share), np.nan)
    summ_rep = replicate[replicate["left_out"] == "SUMMARY"].set_index(["method", "feature_id", "term"])

    # ---- the separate dimensions and the rule
    conf = {}
    if inp.support:
        sup = _read_tsv(cp[inp.support])
        if "assignment_confidence" in sup.columns:
            conf = dict(zip(sup["feature_id"], sup["assignment_confidence"].astype(str)))
    cv = {}
    if inp.boot_abundance:
        ba = _read_tsv(cp[inp.boot_abundance], ("feature_id", "sample"))
        if {"boot_cv", "sample"} <= set(ba.columns):
            ba = ba[ba["sample"].isin(samples)]
            cv = pd.to_numeric(ba["boot_cv"], errors="coerce").groupby(ba["feature_id"]).median().to_dict()
    res_stab = {}
    if inp.resolution:
        rt = _read_tsv(cp[inp.resolution], ("parent", "term"))
        if {"parent", "term", "resolution_stability"} <= set(rt.columns):
            res_stab = {(a, b): c for a, b, c in zip(rt["parent"], rt["term"], rt["resolution_stability"])}
    parent = _parents(used)
    summary_rows = []
    for _, c in concord.iterrows():
        f, t = c["feature_id"], c["term"]
        sig = {m: c[f"sign_{m}"] for m in methods if isinstance(c["methods_significant"], str)
               and m in c["methods_significant"].split(",")}
        stab = {m: stability_lookup[(m, f, t)] for m in methods if (m, f, t) in stability_lookup}
        stab_d = {m: {k: s[k] for k in ("boot_direction_stability", "real_direction_stability",
                                         "boot_significance_stability", "real_significance_stability")}
                  for m, s in stab.items()}
        reps = [summ_rep.loc[(m, f, t)] for m in (sig or stab) if (m, f, t) in summ_rep.index]
        rep_ok = bool(reps) and all(not bool(r["single_sample_driven"]) and np.isfinite(r["loo_direction_stability"])
                                    and r["loo_direction_stability"] >= LOO_DIRECTION_MIN for r in reps)
        assign = conf.get(f, "not_assessed")
        med_cv = cv.get(f, np.nan)
        ab = classify_cv(med_cv)
        rs = res_stab.get((parent.get(f, f), t), "not_assessed")
        conclusion, why = conclude(sig, stab_d, assign, ab, rep_ok)
        basis = []
        if c["concordance"] in CONCORDANCE_CLASSES[:2]:
            basis.append("method agreement")
        if sig and all(np.isfinite(stab_d.get(m, {}).get(k, np.nan)) and stab_d[m][k] >= DIRECTION_MIN
                       for m in sig for k in ("boot_direction_stability", "real_direction_stability")):
            basis.append("effect-direction stability")
        if assign in ("high", "moderate"):
            basis.append("read-assignment stability")
        if rs == "stable":
            basis.append("clustering-resolution stability")
        if rep_ok:
            basis.append("biological replicate consistency")
        support = "; ".join(f"{m}: effect={_fmt(c[f'effect_{m}'])} ({SCALE[m]}) q={_fmt(c[f'q_{m}'], 2)}"
                            for m in sig) or "none"
        summary_rows.append({"feature_id": f, "label": c["label"], "term": t, "assignment_confidence": assign,
                             "abundance_stability": ab, "abundance_boot_cv_median": med_cv, "da_support": support,
                             "n_methods_significant": len(sig), "cross_method_agreement": c["concordance"],
                             "replicate_consistency": rep_ok, "resolution_stability": rs,
                             "robustness_basis": ";".join(basis), "conclusion": conclusion,
                             "unmet_conditions": "; ".join(why) if conclusion == "exploratory" else ""})
    summary = pd.DataFrame(summary_rows)

    # ---- write
    write_tsv(stage / "fit_status.tsv", status, ["one row per method x analysis; failed fits are excluded from "
                                                 "every stability denominator"])
    write_tsv(stage / "differential_abundance_stability.tsv", stability, [
        "method x feature x term; effects stay on each method's own scale (effect_definition)",
        f"direction stability = share of successful fits with the primary's sign; significance stability = share "
        f"with q < {alpha}; {NOT_TESTED_NOTE}", RESAMPLING_NOTE,
        f"filter_significant_all: q < {alpha} at every alternative prevalence cut ({','.join(map(str, cuts)) or 'none'})"])
    write_tsv(stage / "cross_method_concordance.tsv", concord, [
        f"methods: {', '.join(methods)}; significant = q < {alpha} at the primary analysis; no p-value pooling",
        "concordance: all_significant_same_direction (every method) / majority_significant_same_direction / "
        "single_method / conflicting_direction (significant methods disagree) / none", *rho, AGREEMENT_NOTE,
        f"wording: '{preferred_wording(len(methods))}', not validation by independent methods"])
    write_tsv(stage / "replicate_sensitivity.tsv", replicate, [
        f"left-out {'unit (' + args.unit + ')' if args.unit else 'sample'} per row; SUMMARY rows: leave-one-out "
        "direction / significance stability over the successful fits, single_sample_driven = significant at the "
        "primary analysis but not after leaving out some single sample/unit (listed in drivers)",
        "max_single_sample_share = largest share of one group's normalized abundance contributed by one sample; "
        "mde_80 = (1.96 + 0.8416) x primary SE: the effect detectable with 80 % power (two-sided alpha 0.05), "
        "on the method's own scale"])
    write_tsv(stage / "robustness_summary.tsv", summary, [
        "one row per feature x term; the dimensions are separate columns", CONCLUSION_RULE,
        f"replicate_consistency: no method single_sample_driven and leave-one-out direction stability >= "
        f"{LOO_DIRECTION_MIN} in the significant methods (all tested methods when none is significant)",
        f"abundance_stability: median bootstrap CV over samples <= {CV_STABLE} stable, <= {CV_MODERATE} moderate, "
        "else unstable", AGREEMENT_NOTE])
    figs = make_figures(stage / "figures", stability, concord, summary, results, by_kind, methods, abund, groups,
                        alpha)
    versions = software_versions(stage, methods)
    record = {"step": "da-robustness", "methods": methods, "seed": args.seed, "mc_samples": args.mc_samples,
              "parameters": {k: (v if isinstance(v, (int, float, str, bool, type(None))) else str(v))
                             for k, v in sorted(vars(args).items()) if not callable(v)},
              "inputs_sha256": {k: manifest.sha256(v) for k, v in cp.items()},
              "em_uncertainty": ({k: inp.em_record.get(k) for k in ("bootstraps", "realizations", "seed")}
                                 if inp.em_record else None),
              "versions": versions, "features_tested": len(selected), "samples": len(samples),
              "fits": status[["method", "analysis", "status"]].to_dict("records")}
    (stage / "robustness_report.md").write_text(
        render_report(args, methods, sizes, selected, primary, analyses, status, stability, concord, summary,
                      replicate, rho, versions, record, cuts))
    made = ["features.used.tsv", "metadata.used.tsv", "fit_status.tsv", "differential_abundance_stability.tsv",
            "cross_method_concordance.tsv", "replicate_sensitivity.tsv", "robustness_summary.tsv",
            "robustness_report.md", *figs]
    for j, (s, r) in zip(jobs, done):
        rel = Path("fits") / j.outdir.name
        for name in ("results.tsv", "warnings.txt") + (("session_info.txt",) if j.analysis.kind == "primary" else ()):
            if (j.outdir / name).is_file() and (name != "results.tsv" or r is not None):
                made.append(str(rel / name))
    n_ok = int((status["status"] == "ok").sum())
    log.info("da-robustness: %d of %d fits succeeded; %d high-confidence, %d exploratory feature x term",
             n_ok, len(status), int((summary["conclusion"] == "high_confidence").sum()),
             int((summary["conclusion"] == "exploratory").sum()))
    return made, record


def _parents(used: pd.DataFrame, column: str = "parent_feature") -> dict:
    if column not in used.columns:
        return dict(zip(used["feature_id"], used["feature_id"]))
    par = used[column].astype(str).where(used[column].notna() & ~used[column].astype(str).isin(["", "NA", "nan"]),
                                         used["feature_id"])
    return dict(zip(used["feature_id"], par))


def software_versions(stage: Path, methods: list) -> dict:
    import matplotlib
    from . import __version__
    v = {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
         "matplotlib": matplotlib.__version__, "unclassified_meta": __version__}
    pat = re.compile(r"\b(ANCOMBC|ALDEx2|MicrobiomeStat|lme4|lmerTest)_([0-9][0-9.\-]*)")
    for m in methods:
        f = stage / "fits" / f"{m}__primary" / "session_info.txt"
        if not f.is_file():
            continue
        txt = f.read_text(errors="replace")
        r = re.search(r"R version ([0-9.]+)", txt)
        pk = dict(pat.findall(txt))
        v[f"R ({m})"] = ", ".join(([f"R {r.group(1)}"] if r else []) + [f"{k} {x}" for k, x in sorted(pk.items())])
    return v


# ------------------------------------------------------------------------------------------------ figures

def make_figures(figdir: Path, stability, concord, summary, results, by_kind, methods, abund, groups,
                 alpha) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figdir.mkdir(parents=True, exist_ok=True)
    made = []
    # top features: the smallest primary q over the methods
    best = stability.groupby(["feature_id", "term"])["primary_q"].min().sort_values()
    top = list(best.index[:12])
    lab = dict(zip(zip(stability["feature_id"], stability["term"]), stability["label"]))
    colors = {"ancombc2": "#1f77b4", "aldex2": "#ff7f0e", "linda": "#2ca02c"}

    def save(fig, name):
        fig.tight_layout()
        fig.savefig(figdir / name, dpi=150)
        plt.close(fig)
        made.append(f"figures/{name}")

    # 1. bootstrap effect distributions
    n = max(len(top), 1)
    ncol = min(4, n)
    fig, axes = plt.subplots(int(np.ceil(n / ncol)), ncol, figsize=(3.2 * ncol, 2.6 * np.ceil(n / ncol)),
                             squeeze=False)
    boot_eff = {m: [_index(results[(m, a.name)])["effect"] for a in by_kind.get("boot", []) if (m, a.name) in results]
                for m in methods}
    prim_eff = {m: _index(results[(m, "primary")])["effect"] for m in methods if (m, "primary") in results}
    for ax, k in zip(axes.flat, top):
        data, pos, prim_pts = [], [], []
        for i, m in enumerate(methods):
            effs = [e for e in (b.get(k, np.nan) for b in boot_eff[m]) if np.isfinite(e)]
            if effs:
                data.append(effs)
                pos.append(i)
            if m in prim_eff:
                prim_pts.append((i, prim_eff[m].get(k, np.nan)))
        if data:
            ax.boxplot(data, positions=pos, widths=0.5, showfliers=False)
        for i, e in prim_pts:
            ax.plot(i, e, "D", color=colors.get(methods[i], "k"), ms=5)
        ax.axhline(0, color="grey", lw=0.6)
        ax.set_xticks(range(len(methods)), methods, fontsize=7)
        ax.set_title(f"{lab.get(k, k[0])}\n{k[1]}", fontsize=7)
    for ax in list(axes.flat)[len(top):]:
        ax.axis("off")
    if not by_kind.get("boot"):
        axes.flat[0].text(0.5, 0.5, "no bootstrap fits", ha="center", transform=axes.flat[0].transAxes)
    fig.suptitle("bootstrap effects per method (box) and primary effect (diamond); scales differ per method",
                 fontsize=8)
    save(fig, "bootstrap_effects.png")

    # 2. direction-stability heatmap: features x method/analysis
    keys = list(best.index[:30])
    cols = [(m, c) for m in methods for c in ("boot_direction_stability", "real_direction_stability")]
    mat = np.full((len(keys), len(cols)), np.nan)
    st = stability.set_index(["method", "feature_id", "term"])
    for i, k in enumerate(keys):
        for j, (m, c) in enumerate(cols):
            if (m, *k) in st.index:
                mat[i, j] = st.loc[(m, *k), c]
    fig, ax = plt.subplots(figsize=(1.0 + 0.6 * len(cols), 1.2 + 0.25 * max(len(keys), 1)))
    im = ax.imshow(mat if len(keys) else np.full((1, len(cols)), np.nan), vmin=0, vmax=1, cmap="viridis",
                   aspect="auto")
    ax.set_xticks(range(len(cols)), [f"{m}\n{c.split('_')[0]}" for m, c in cols], fontsize=6)
    ax.set_yticks(range(len(keys)), [f"{lab.get(k, k[0])} [{k[1]}]" for k in keys], fontsize=6)
    fig.colorbar(im, ax=ax, label="direction stability")
    ax.set_title("share of resampled fits with the primary direction", fontsize=8)
    save(fig, "direction_stability_heatmap.png")

    # 3. pairwise concordance of primary effects
    pairs = [(a, b) for i, a in enumerate(methods) for b in methods[i + 1:]]
    fig, axes = plt.subplots(1, max(len(pairs), 1), figsize=(3.6 * max(len(pairs), 1), 3.4), squeeze=False)
    for ax, (a, b) in zip(axes.flat, pairs):
        x, y = concord[f"effect_{a}"].astype(float), concord[f"effect_{b}"].astype(float)
        sa = concord[f"q_{a}"].astype(float) < alpha
        sb = concord[f"q_{b}"].astype(float) < alpha
        for mask, col, name in ((sa & sb, "#d62728", "both"), (sa ^ sb, "#ff9896", "one"),
                                (~sa & ~sb, "#bbbbbb", "neither")):
            ax.scatter(x[mask], y[mask], s=12, c=col, label=f"significant in {name}")
        ax.axhline(0, color="grey", lw=0.5)
        ax.axvline(0, color="grey", lw=0.5)
        ax.set_xlabel(f"{a} effect ({SCALE[a]})", fontsize=8)
        ax.set_ylabel(f"{b} effect ({SCALE[b]})", fontsize=8)
        ax.legend(fontsize=6)
    if not pairs:
        axes.flat[0].text(0.5, 0.5, "one method only", ha="center", transform=axes.flat[0].transAxes)
    save(fig, "concordance.png")

    # 4. per-sample normalized abundance by group for the top features
    feats = list(dict.fromkeys(k[0] for k in top))[:12]
    n = max(len(feats), 1)
    ncol = min(4, n)
    fig, axes = plt.subplots(int(np.ceil(n / ncol)), ncol, figsize=(3.0 * ncol, 2.5 * np.ceil(n / ncol)),
                             squeeze=False)
    levels = list(dict.fromkeys(groups.astype(str)))
    for ax, f in zip(axes.flat, feats):
        for i, g in enumerate(levels):
            smp = [s for s in groups.index[groups.astype(str) == g] if s in abund.columns]
            v = abund.loc[f, smp].to_numpy(dtype=float) if f in abund.index else np.array([])
            v = np.log10(v + 1)              # log10(x + 1) on a linear axis (no mathtext tick labels: fast)
            ax.scatter(np.full(len(v), i) + np.linspace(-0.12, 0.12, max(len(v), 1))[:len(v)], v, s=12)
            if len(v) > 1:                   # name the largest sample of each group: the one a call may rest on
                ax.annotate(smp[int(np.argmax(v))], (i, v.max()), fontsize=5)
        ax.set_ylabel("log10(abundance + 1)", fontsize=6)
        ax.set_xticks(range(len(levels)), levels, fontsize=7)
        ax.set_title(str(lab.get(next((k for k in top if k[0] == f), (f, "")), f)), fontsize=7)
    for ax in list(axes.flat)[len(feats):]:
        ax.axis("off")
    fig.suptitle("per-sample normalized abundance (each point is one biological sample)", fontsize=8)
    save(fig, "replicate_abundance.png")
    return made


# ------------------------------------------------------------------------------------------------ report

def _md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_none_\n"
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]).replace("|", "/") for c in cols) + " |")
    return "\n".join(lines) + "\n"


def render_report(args, methods, sizes, selected, primary, analyses, status, stability, concord, summary,
                  replicate, rho, versions, record, cuts) -> str:
    alpha = float(args.alpha)
    nm = len(methods)
    low = {g: n for g, n in sizes.items() if n <= 4}
    ok = status[status["status"] == "ok"]
    st = stability.set_index(["method", "feature_id", "term"])
    rep = replicate[replicate["left_out"] == "SUMMARY"].set_index(["method", "feature_id", "term"])
    out = ["# Differential-abundance robustness", ""]
    out += ["## Design", "",
            f"- samples per group ({args.group}): " + ", ".join(f"{g} = {n}" for g, n in sizes.items()),
            f"- formula: ~ {args.formula or args.group}; reference level: {args.reference_level or 'first level'}; "
            f"random effects: {args.rand_formula or 'none'} (ANCOM-BC2 and LinDA"
            f"{'; ALDEx2 fits the fixed effects only and treats repeated samples of one unit as independent' if args.rand_formula and 'aldex2' in methods else ''})",
            f"- leave-one-out unit: {args.unit or 'sample'}",
            f"- methods: {', '.join(METHOD_NAMES[m] for m in methods)}",
            f"- features tested: {len(selected)} of {len(primary)} (one shared prevalence filter, count > 0 in >= "
            f"{args.prv_cut} of the samples; the methods' own filters off)",
            "- analyses: " + ", ".join(f"{k} {sum(1 for a in analyses if a.kind == k)}"
                                       for k in ("primary", "unique", "boot", "real", "loo", "prv")
                                       if any(a.kind == k for a in analyses)),
            f"- fits: {len(ok)} of {len(status)} succeeded", ""]
    if low:
        mde = st["primary_se"].groupby(level=0).median() * Z_MDE80
        out += ["## Low replication", "",
                f"At least one group has 4 or fewer samples ({', '.join(f'{g} = {n}' for g, n in low.items())}). "
                "With so few biological replicates the tests have little power: only large effects can reach "
                f"q < {alpha}, variance estimates are unstable, and one sample can carry a result "
                "(see replicate_sensitivity.tsv).",
                "",
                "- A non-significant result is not evidence of no effect.",
                "- Minimal detectable effect at 80 % power (median over features, (1.96 + 0.8416) x SE, each "
                "method's own scale): " + ", ".join(f"{m} {_fmt(v)} ({SCALE[m]})" for m, v in mde.items()),
                "- Bootstrap tables and posterior realizations are resamples of the same samples, not "
                "replicates; they cannot add power.", ""]
    out += ["## Wording", "",
            f"Agreement is reported as \"{preferred_wording(nm)}\". The methods analyse the same reads, samples "
            "and model, so agreement is not independent validation and must not be described as such.", "",
            "## Conclusion rule", "", CONCLUSION_RULE, ""]

    def rows(kind):
        sub = summary[summary["conclusion"] == kind]
        recs = []
        for _, s in sub.iterrows():
            k = (s["feature_id"], s["term"])
            per = []
            for m in methods:
                if (m, *k) not in st.index:
                    continue
                r = st.loc[(m, *k)]
                per.append(f"{m}: {_fmt(r['primary_effect'])} {SCALE[m]}, q={_fmt(r['primary_q'], 2)}, "
                           f"dir boot/real {_fmt(r['boot_direction_stability'], 2)}/"
                           f"{_fmt(r['real_direction_stability'], 2)}")
            reps = [f"{m} loo dir {_fmt(rep.loc[(m, *k), 'loo_direction_stability'], 2)}"
                    + (f" (drivers {rep.loc[(m, *k), 'drivers']})" if rep.loc[(m, *k), "drivers"] else "")
                    for m in methods if (m, *k) in rep.index]
            recs.append({"feature": s["label"], "term": s["term"], "methods": "<br>".join(per),
                         "abundance": f"{s['abundance_stability']} (CV {_fmt(s['abundance_boot_cv_median'], 2)})",
                         "assignment": s["assignment_confidence"], "replicates": "; ".join(reps),
                         "agreement": s["cross_method_agreement"],
                         **({"unmet": s["unmet_conditions"]} if kind == "exploratory" else {})})
        return pd.DataFrame(recs)

    hc, ex = rows("high_confidence"), rows("exploratory")
    out += ["## High-confidence findings", "",
            f"{len(hc)} feature x term. These are \"{preferred_wording(nm)}\" and stable to read resampling, "
            "read assignment and leaving out single samples." if nm == 3 else f"{len(hc)} feature x term.", "",
            _md_table(hc), "## Exploratory findings", "",
            f"{len(ex)} feature x term significant in at least one method but not meeting the rule above.", "",
            _md_table(ex),
            f"{int((summary['conclusion'] == 'not_significant').sum())} feature x term were not significant in any "
            "method.", ""]
    out += ["## What each robustness claim refers to", "",
            "- **method agreement**: q < alpha in at least two methods, same direction, on the same data; not "
            "independent validation.",
            "- **effect-direction stability**: share of bootstrap (reads resampled within samples, EM refitted) and "
            "posterior-realization (ambiguous reads re-drawn) fits with the primary direction. "
            + RESAMPLING_NOTE + ".",
            "- **read-assignment stability**: assignment_confidence of em-uncertainty (unique-read support, "
            "bootstrap CV, agreement of the mapping criteria).",
            "- **clustering-resolution stability**: the same parent feature significant in the same direction at "
            "every clustering resolution (da-resolutions).",
            "- **biological replicate consistency**: no single sample (or unit) carries the call, judged by "
            "leave-one-out refits.", "",
            "## Cross-method rank correlation", "", *[f"- {x}" for x in rho], "", AGREEMENT_NOTE + ".", ""]
    failed = status[status["status"] != "ok"]
    out += ["## Failed fits", "", "Failed fits are excluded from every stability denominator.", "",
            _md_table(failed[["method", "analysis", "message"]]) if len(failed) else "_none_\n"]
    out += ["## Software versions", "", *[f"- {k}: {v}" for k, v in versions.items()], ""]
    p = record["parameters"]
    out += ["## Parameters and seeds", "",
            *[f"- {k}: {p[k]}" for k in ("alpha", "p_adjust", "prv_cut", "prv_cuts", "bootstrap_fits",
                                         "realization_fits", "mc_samples", "seed") if k in p],
            "- ANCOM-BC2 seed: fixed in r/ancombc2.R (set.seed(11))",
            f"- em-uncertainty: {record['em_uncertainty'] or 'not used (--counts)'}", ""]
    text = "\n".join(out)
    if "validated by" in text.lower():          # the wording rule (docs/DESIGN_robustness.md, principle 4)
        raise StepError("internal: the report text would describe method agreement as validation")
    return text


# ------------------------------------------------------------------------------------------------ publishing

def _publish(stage: Path, outdir: Path, made: list, previous: dict, record: dict, inputs, command: str) -> None:
    from .io_utils import refuse_overlap
    refuse_overlap([outdir / f for f in made], inputs, command)
    manifest.remove_stale(outdir, previous, made)
    for f in made:
        dst = outdir / f
        dst.parent.mkdir(parents=True, exist_ok=True)
        os.replace(stage / f, dst)
    # top-level files of the previous generation that this one did not write again (unchanged ones only)
    for name in sorted(set(previous) - set(made)):
        f = outdir / name
        if "/" in name or manifest.STALE_NAME.fullmatch(name) or f.is_symlink() or not f.is_file():
            continue
        if manifest.sha256(f) == previous[name]:
            f.unlink()
    fits = outdir / "fits"
    if fits.is_dir() and not fits.is_symlink():
        for d in fits.iterdir():
            if d.is_dir() and not d.is_symlink() and not any(d.iterdir()):
                d.rmdir()
    manifest.seal(outdir, made, record)


# ------------------------------------------------------------------------------------------------ da-resolutions

def _parse_runs(values) -> list[tuple[str, str]]:
    runs = []
    for v in values or []:
        if "=" not in v:
            raise StepError(f"--run {v}: give LABEL=DIR, e.g. 98=da_rob_98")
        lab, d = v.split("=", 1)
        if not lab.strip() or not d.strip():
            raise StepError(f"--run {v}: give LABEL=DIR")
        runs.append((lab.strip(), d.strip()))
    labs = [r[0] for r in runs]
    if len(set(labs)) != len(labs):
        raise StepError("--run: duplicate labels")
    if len(runs) < 2:
        raise StepError("da-resolutions needs at least two --run LABEL=DIR")
    return runs


def _feature_sign(row, methods) -> str:
    sig = [m for m in str(row.get("methods_significant") or "").split(",") if m]
    s = {float(row.get(f"sign_{m}")) for m in sig if pd.notna(row.get(f"sign_{m}"))}
    if not s:
        return ""
    return "+" if s == {1.0} else "-" if s == {-1.0} else "+-"


def resolution_table(runs: dict, parent_column: str = "parent_feature") -> pd.DataFrame:
    """runs: {label: (cross_method_concordance DataFrame, features.used DataFrame)} -> one row per parent x term."""
    order = {c: i for i, c in enumerate(CONCORDANCE_CLASSES)}
    per = {}
    for lab, (conc, used) in runs.items():
        parent = _parents(used, parent_column)
        methods = [c[5:] for c in conc.columns if c.startswith("sign_")]
        conc = conc.assign(parent=conc["feature_id"].map(lambda f: parent.get(f, f)),
                           sign=[_feature_sign(r, methods) for _, r in conc.iterrows()])
        for (par, term), g in conc.groupby(["parent", "term"], sort=False):
            best = min(g["concordance"], key=lambda c: order.get(c, len(order)))
            signs = sorted(set("".join(g.loc[g["n_methods_significant"] > 0, "sign"])))
            per.setdefault((par, term), {})[lab] = {
                "n": len(g), "concordance": best, "significant": bool((g["n_methods_significant"] > 0).any()),
                "signs": "".join(signs)}
    labels = list(runs)
    rows = []
    for (par, term), d in per.items():
        row = {"parent": par, "term": term}
        for lab in labels:
            x = d.get(lab, {"n": 0, "concordance": "not_present", "significant": False, "signs": ""})
            row[f"n_features_{lab}"] = x["n"]
            row[f"concordance_{lab}"] = x["concordance"]
            row[f"significant_{lab}"] = x["significant"]
            row[f"signs_{lab}"] = x["signs"]
        sig_at = [lab for lab in labels if d.get(lab, {}).get("significant")]
        all_signs = {c for lab in sig_at for c in d[lab]["signs"]}
        row["consistent_direction"] = bool(sig_at) and len(all_signs) == 1
        row["significant_at"] = ",".join(sig_at)
        row["resolution_stability"] = ("not_significant" if not sig_at else
                                       "stable" if len(sig_at) == len(labels) and row["consistent_direction"]
                                       else "resolution_dependent")
        rows.append(row)
    return pd.DataFrame(rows)


def run_resolutions(args) -> Path:
    from .io_utils import make_staging_dir, remove_staging_dir
    runs = _parse_runs(args.run)
    tables, shas = {}, {}
    for lab, d in runs:
        rec = manifest.verify(d, what="da-robustness results")
        if rec.get("step") != "da-robustness":
            raise StepError(f"--run {lab}={d}: not a da-robustness folder (step {rec.get('step')})")
        for n in ("cross_method_concordance.tsv", "features.used.tsv"):
            if n not in rec["artifacts"]:
                raise StepError(f"{d}: the manifest does not seal {n}")
            shas[str(Path(d) / n)] = rec["artifacts"][n]
    outdir = Path(args.outdir)
    previous = manifest.previous_artifacts(outdir)
    manifest.begin(outdir, "da-resolutions")
    stage = make_staging_dir(outdir, "da-resolutions")
    try:
        cp = manifest.private_copies(shas, stage / "inputs")
        for lab, d in runs:
            conc = _read_tsv(cp[str(Path(d) / "cross_method_concordance.tsv")], ("feature_id", "term"))
            used = _read_tsv(cp[str(Path(d) / "features.used.tsv")], ("feature_id", args.parent_column))
            tables[lab] = (conc, used)
        out = resolution_table(tables, args.parent_column)
        write_tsv(stage / "multi_resolution_robustness.tsv", out, [
            "one row per parent feature x term; per resolution: features of that parent, their best cross-method "
            "concordance, whether any is significant (q < alpha in at least one method), and the directions of "
            "the significant ones",
            "resolution_stability: stable = significant at every resolution in one direction; resolution_dependent "
            "= significant at some resolutions only, or in different directions; not_significant = at none",
            f"parent column: {args.parent_column} of each run's features.used.tsv (feature_id when empty)"])
        record = {"step": "da-resolutions", "runs": {lab: str(Path(d).resolve()) for lab, d in runs},
                  "parent_column": args.parent_column,
                  "inputs_sha256": {k: manifest.sha256(v) for k, v in cp.items()}}
        _publish(stage, outdir, ["multi_resolution_robustness.tsv"], previous, record, list(shas), "da-resolutions")
    finally:
        remove_staging_dir(stage)
    return outdir / "multi_resolution_robustness.tsv"
