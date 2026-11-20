"""Thin wrapper around the ANCOM-BC2 R script (r/ancombc2.R)."""
from __future__ import annotations

import os
import shutil
import json
import logging
import subprocess
from pathlib import Path

from .runner import StepError

log = logging.getLogger("unclassified_meta")
SCRIPT = Path(__file__).parent / "r" / "ancombc2.R"


def find_rscript(explicit: str | None = None) -> str:
    cand = explicit or os.environ.get("UNCLASSIFIED_META_RSCRIPT") or shutil.which("Rscript")
    if not cand or not (Path(cand).exists() or shutil.which(cand)):
        raise StepError("Rscript not found. Create the R environment (micromamba create -f environment-r.yml) and "
                        "pass --rscript <env>/bin/Rscript or set UNCLASSIFIED_META_RSCRIPT.")
    return cand


def has_ancombc(rscript: str) -> bool:
    p = subprocess.run([rscript, "-e", "suppressPackageStartupMessages(library(ANCOMBC))"],
                       capture_output=True, text=True)
    return p.returncode == 0


def verify_inputs(counts: str | Path, features: str | Path, unsealed_ok: bool = False) -> dict:
    """Tables made by panel-quant or taxa-table are used only as their sealed generation: each consumed file must
    be listed, with an unchanged checksum, in its folder's manifest (provenance.json of panel-quant, MANIFEST.json
    of taxa-table). A folder this program wrote (producer marker, ledger.tsv) without a complete seal is refused.
    Tables made by hand carry no manifest and are used only with --unsealed-input.
    Returns {file: sha256 recorded by its manifest, or None for a table used unsealed}."""
    from . import manifest
    expected = {}
    for f in (counts, features):
        folder = Path(f).parent
        expected[str(f)] = None
        if any((folder / n).exists() for n in ("provenance.json", manifest.NAME)):
            rec = manifest.verify_file(f, what="count tables")
            expected[str(f)] = rec["artifacts"][Path(f).name]
        elif manifest.is_generated(folder) or (folder / "ledger.tsv").exists():
            raise StepError(f"{folder} was written by panel-quant or taxa-table but has no complete manifest "
                            "(provenance.json / MANIFEST.json): the tables were not published completely; run that "
                            "step again")
        elif not unsealed_ok:
            raise StepError(f"{f} has no manifest (not made by panel-quant or taxa-table); pass --unsealed-input to "
                            "use a table made by hand as it is")
        else:
            log.warning("using %s without a manifest (--unsealed-input): no completeness or checksum check", f)
    return expected


# every file name the R script can write; a rerun removes the ones it did not write again
DA_OUTPUT_PATTERNS = ("ancombc2_results.tsv", "structural_zeros.tsv", "session_info.txt", "warnings.txt", "*.png",
                      "*.pdf")             # mi_results.tsv is removed only as a sealed output (manifest.remove_stale)


P_ADJUST_METHODS = ("holm", "bh", "fdr", "by", "bonferroni", "none")      # what p_adjust implements (lower case)


def check_p_adjust(method: str) -> None:
    """Refuse, before any fit runs, a --p-adjust that the Python pooling (p_adjust) cannot apply."""
    if str(method).lower() not in P_ADJUST_METHODS:
        raise StepError(f"--p-adjust {method} is not available with --imputations (use holm, BH, BY, bonferroni, "
                        "none)")


def p_adjust(p, method: str):
    """R's p.adjust for holm, BH/fdr, BY, bonferroni and none (NaN entries are left out)."""
    import numpy as np
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    x = p[ok]
    n = len(x)
    if n == 0:
        return out
    m = method.lower()
    if m == "none":
        adj = x
    elif m == "bonferroni":
        adj = np.minimum(1, x * n)
    elif m == "holm":
        o = np.argsort(x)
        a = np.minimum(1, np.maximum.accumulate((n - np.arange(n)) * x[o]))
        adj = np.empty(n); adj[o] = a
    elif m in ("bh", "fdr", "by"):
        o = np.argsort(x)[::-1]
        q = (np.sum(1.0 / np.arange(1, n + 1)) if m == "by" else 1.0)
        a = np.minimum(1, np.minimum.accumulate(q * n / np.arange(n, 0, -1) * x[o]))
        adj = np.empty(n); adj[o] = a
    else:
        raise StepError(f"--p-adjust {method} is not available with --imputations (use holm, BH, BY, bonferroni, none)")
    out[ok] = adj
    return out


def read_results(path) -> "pd.DataFrame":
    """An ANCOM-BC2 result table: leading '#' lines are notes; a '#' inside a row (a label such as 'strain #2')
    is data."""
    import pandas as pd
    skip = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            skip += 1
    return pd.read_csv(path, sep="\t", skiprows=skip, quoting=3)        # the R script writes unquoted fields


# |W| at or below this is a zero statistic: a two-sided t test then has p = 1 legitimately. Above it a valid test has
# p < 1 (|W| = 1e-12 gives p = 1 - 8e-13, which R writes as < 1), so p = 1 there is ANCOM-BC2's fill for no test.
_ZERO_W = 1e-12
MI_STATUS = ("pooled", "untestable", "partly_untestable", "not_estimated", "incomplete")


def _untestable(p, w, lfc, se):
    """ANCOM-BC2 fits it could not test: p_hat = 2 pt(|W|, dof) is NaN when the feature's residual df is 0 (lm
    saturated: the feature is nonzero in no more samples than the design has columns) and ancombc2 then sets
    p_hat[is.na(p_hat)] = 1. A real two-sided test has p = 1 only at W = 0, so p >= 1 with a nonzero (or missing)
    W marks the fill. W = lfc / se in ANCOM-BC2; it is used where the W column is missing."""
    import numpy as np
    w = np.where(np.isfinite(w), w, lfc / np.where(se > 0, se, np.nan))
    return (p >= 1) & ~(np.abs(w) <= _ZERO_W)


def combine_imputations(results: list, alpha: float, method: str):
    """Rubin's rules over the ANCOM-BC2 result tables of the imputations, term by term (intercept left out).

    Every feature of any fit is reported. A fit counts for a feature only when its estimate, its standard error
    (finite, > 0) and its p are finite and ANCOM-BC2 could test it (see _untestable). A feature is pooled only when
    every fit counts; otherwise its p is set to 1 and it stays in the multiple-testing adjustment, as ANCOM-BC2 does
    for a failed test. Pooling only the fits that could test a feature would condition on the imputations that
    happened to give it enough nonzero samples, so it is not done. status_<term> says why a feature was not pooled:
    untestable (no fit could test it), partly_untestable (some could), not_estimated (no finite estimate in any fit:
    filtered by prevalence or a failed fit), incomplete (missing, or without a finite p, in some fits).

    The complete-data degrees of freedom of each feature are recovered from the counted fits' own t tests (W and p)
    and their median is used, so identical imputations reproduce ANCOM-BC2's test; when no fit gives one (all W = 0,
    where p = 1 anyway, or all p underflowed to 0) the reference is normal. fmi_<term> is lambda = (1 + 1/M) B / T,
    the share of the total variance due to the read assignment (Rubin's gamma adds a small-df correction and is
    larger). diff_robust: diff, and ANCOM-BC2's own diff_robust (significant and stable under its pseudo-count
    sensitivity analysis) in at least half of the fits.
    """
    import numpy as np
    import pandas as pd
    from . import em
    if len(results) < 2:
        raise StepError(f"combining imputations needs at least 2 fits, got {len(results)}")
    base = results[0]
    need = {c for c in base.columns if c.split("_", 1)[0] in ("lfc", "se", "W", "p") and "(Intercept)" not in c}
    for m, r in enumerate(results[1:], 2):
        missing = sorted(need - set(r.columns))
        extra = sorted(c for c in r.columns if c.startswith("lfc_") and c not in base.columns)
        if missing or extra:
            raise StepError(f"imputation fit {m} does not have the coefficients of fit 1 "
                            f"(missing {', '.join(missing) or '-'}; extra {', '.join(extra) or '-'})")
    ids = list(dict.fromkeys(f for r in results for f in r["feature_id"]))
    rs = [r.drop_duplicates("feature_id").set_index("feature_id").reindex(ids) for r in results]
    meta_cols = [c for c in ("label", "identity_class") if c in base.columns]
    out = pd.DataFrame({"feature_id": ids})
    for c in meta_cols:
        vals = pd.concat([r[c] for r in rs], axis=1).bfill(axis=1).iloc[:, 0]
        out[c] = vals.to_numpy()
    terms = [c[4:] for c in base.columns if c.startswith("lfc_") and c != "lfc_(Intercept)"]
    true_frac = lambda col: np.column_stack([r[col].astype(str).str.upper().eq("TRUE").to_numpy() for r in rs]
                                            ).mean(axis=1) if col in base.columns else np.full(len(ids), np.nan)
    col = lambda name: np.column_stack([r[name].to_numpy(dtype=float) for r in rs])
    for t in terms:
        if f"se_{t}" not in base.columns:
            raise StepError(f"the imputation fits give lfc_{t} without se_{t}")
        est, se = col(f"lfc_{t}"), col(f"se_{t}")
        estimated = np.isfinite(est) & np.isfinite(se) & (se > 0)
        untestable = np.zeros_like(estimated)
        fin = estimated.copy()
        if f"p_{t}" in base.columns:              # a fit whose own test failed is not a successful fit
            pf = col(f"p_{t}")
            wf = col(f"W_{t}") if f"W_{t}" in base.columns else np.full_like(pf, np.nan)
            untestable = estimated & np.isfinite(pf) & _untestable(pf, wf, est, se)
            fin &= np.isfinite(pf) & ~untestable
        n_fits = fin.sum(axis=1)
        m = len(rs)
        n_est, n_unt = estimated.sum(axis=1), untestable.sum(axis=1)
        status = np.select([n_fits == m, n_est == 0, (n_unt > 0) & (n_unt == n_est), n_unt > 0],
                           ["pooled", "not_estimated", "untestable", "partly_untestable"], "incomplete")
        est = np.where(fin, est, np.nan)
        dfc = np.full(len(ids), np.inf)
        if f"W_{t}" in base.columns and f"p_{t}" in base.columns:
            for i in np.flatnonzero(n_fits == m):     # only pooled features need one, from the counted fits only
                d = [em.df_from_p(float(wf[i, j]), float(pf[i, j])) for j in range(m) if fin[i, j]]
                d = [x for x in d if not np.isnan(x)]
                dfc[i] = float(np.median(d)) if d else np.inf
        c = em.rubin(est, se, dfc)
        pv = np.where((n_fits == m) & np.isfinite(c["p"]), c["p"], 1.0)
        q = p_adjust(pv, method)
        out[f"lfc_{t}"], out[f"se_{t}"], out[f"df_{t}"] = c["estimate"], c["se"], c["df"]
        out[f"fmi_{t}"], out[f"p_{t}"], out[f"q_{t}"] = c["fmi_between"], pv, q
        out[f"n_fits_{t}"] = n_fits
        out[f"n_untestable_{t}"] = n_unt
        out[f"status_{t}"] = status
        out[f"diff_{t}"] = (q < alpha) & (n_fits == len(rs))
        out[f"diff_frac_{t}"] = true_frac(f"diff_{t}")
        out[f"robust_frac_{t}"] = true_frac(f"diff_robust_{t}")
        out[f"passed_ss_frac_{t}"] = true_frac(f"passed_ss_{t}")
        out[f"diff_robust_{t}"] = out[f"diff_{t}"] & (out[f"robust_frac_{t}"] >= 0.5)
    return out, terms


def imputation_inputs(args, expected: dict) -> dict:
    """{imputation table: sha256} of the generation that sealed --counts, read from ONE manifest record that also
    vouches for the counts table itself (so the expected table and its imputations are one generation)."""
    from . import manifest
    folder = Path(args.counts).parent
    name = "provenance.json" if (folder / "provenance.json").exists() else manifest.NAME
    rec = manifest.verify(folder, name, what="count tables")
    if expected.get(str(args.counts)) is None or rec["artifacts"].get(Path(args.counts).name) != expected[str(args.counts)]:
        raise StepError(f"--imputations: {args.counts} is not the table sealed with its imputations; run the step again")
    n = int(rec.get("imputations") or 0)
    if n < 2:
        raise StepError(f"--imputations: {folder} holds {n} imputation table(s); make at least 2 (taxa-table / "
                        "panel-quant --imputations, then merge)")
    names = [f"imputations/counts_{m + 1:03d}.tsv" for m in range(n)]
    missing = [x for x in names if x not in rec["artifacts"]]
    if missing:
        raise StepError(f"--imputations: {folder} does not seal {', '.join(missing[:3])}")
    return {str(folder / x): rec["artifacts"][x] for x in names}


def _fit_imputations(args, cmd: list, cp: dict, imp_files: list, stage: Path, outdir: Path) -> list:
    """Fit ANCOM-BC2 on every imputation table (private copies) and combine the fits (mi_results.tsv)."""
    import pandas as pd
    rows = list(pd.read_csv(cp[str(args.counts)], sep="\t")["feature_id"])
    results, made = [], []
    (stage / "imputations").mkdir()
    n = len(imp_files)
    for m, f in enumerate(imp_files):
        table = cp[f]
        from .taxa import read_draw
        if list(read_draw(table).index) != rows:              # whole reads only, same rows
            raise StepError(f"{f}: its rows differ from {args.counts}")
        sub = stage / f"fit_{m + 1:03d}"
        sub.mkdir()
        c = list(cmd)
        c[c.index("--counts") + 1] = str(table)
        c[c.index("--outdir") + 1] = str(sub)
        log.info("da --imputations: fit %d of %d", m + 1, n)
        p = subprocess.run(c, capture_output=True, text=True)
        if p.returncode != 0:
            print(p.stderr.strip().replace(str(stage), str(outdir)))
            raise StepError(f"ANCOM-BC2 failed on imputation {m + 1} (exit {p.returncode}); see the messages above")
        res = read_results(sub / "ancombc2_results.tsv")
        name = f"imputations/ancombc2_results_{m + 1:03d}.tsv"
        res.to_csv(stage / name, sep="\t", index=False)
        results.append(res)
        made.append(name)
    combined, terms = combine_imputations(results, float(args.alpha), args.p_adjust)
    with open(stage / "mi_results.tsv", "w") as fo:
        fo.write(f"# ANCOM-BC2 on {n} imputations of the ambiguous reads, combined with Rubin's rules "
                 f"(T = W + (1 + 1/M) B; t reference with Barnard-Rubin df from each feature's complete-data df); "
                 f"p adjusted with {args.p_adjust} per term, with p = 1 for a feature not estimated and tested in "
                 "every fit (status: pooled, untestable = ANCOM-BC2 p = 1 with W != 0 in every fit, "
                 "partly_untestable, not_estimated, incomplete); "
                 "fmi = lambda = (1 + 1/M) B / T, the share of the variance due to the read assignment; "
                 "diff_robust: diff and ANCOM-BC2's diff_robust in at least half of the fits\n")
        combined.to_csv(fo, sep="\t", index=False)
    for t in terms:
        log.info("da --imputations: %s: %d of %d features differ (q < %s), %d robust", t,
                 int(combined[f"diff_{t}"].sum()), len(combined), args.alpha, int(combined[f"diff_robust_{t}"].sum()))
    return made


def run_ancombc2(args, inputs=()) -> Path:
    """Run the R script into a private staging folder inside the output folder; only after it succeeded are its
    files checked against the inputs and moved into place, and the output folder is sealed with a manifest."""
    from . import manifest
    from .io_utils import make_staging_dir, refuse_overlap, remove_staging_dir
    if getattr(args, "imputations", False):
        check_p_adjust(args.p_adjust)            # before the previous result is unpublished and any fit runs
    expected = verify_inputs(args.counts, args.features, getattr(args, "unsealed_input", False))
    rscript = find_rscript(args.rscript)
    if not has_ancombc(rscript):
        raise StepError(f"{rscript} has no ANCOMBC package; use the Rscript of the environment-r.yml environment")
    outdir = Path(args.outdir)
    if (outdir / "imputations").is_symlink():
        raise StepError(f"{outdir / 'imputations'} is a link; the per-imputation results must be written inside "
                        "the output folder itself")
    previous = manifest.previous_artifacts(outdir)
    manifest.begin(outdir, "da")
    stage = make_staging_dir(outdir, "da")
    try:
        ledger = Path(args.counts).parent / "ledger.tsv"
        wanted = {**expected, str(args.metadata): None}
        if ledger.exists():
            wanted[str(ledger)] = (manifest.verify_file(ledger)["artifacts"].get("ledger.tsv")
                                   if expected[str(args.counts)] is not None else None)
        imp = imputation_inputs(args, expected) if getattr(args, "imputations", False) else {}
        wanted.update(imp)
        cp = manifest.private_copies(wanted, stage / "inputs")
        cmd = [rscript, str(SCRIPT), "--counts", cp[str(args.counts)], "--features", cp[str(args.features)],
               "--metadata", cp[str(args.metadata)], "--group", args.group, "--outdir", str(stage),
               "--prv-cut", str(args.prv_cut),
               "--lib-cut", str(args.lib_cut), "--alpha", str(args.alpha), "--p-adjust", args.p_adjust,
               "--threads", str(args.threads)]
        if ledger.exists():
            cmd += ["--ledger", str(cp[str(ledger)])]
        if args.formula:
            cmd += ["--formula", args.formula]
        if args.reference_level:
            cmd += ["--reference-level", args.reference_level]
        if args.struc_zero:
            cmd += ["--struc-zero"]
        if getattr(args, "neg_lb", False):
            cmd += ["--neg-lb"]
        if getattr(args, "rand_formula", None):
            cmd += ["--rand-formula", args.rand_formula]
        if getattr(args, "lme_relax_checks", False):
            cmd += ["--lme-relax-checks"]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.stderr.strip():
            # R writes into the staging folder, which is removed below; name the published files instead
            print(p.stderr.strip().replace(str(stage), str(outdir)))
        if p.returncode != 0:
            raise StepError(f"ANCOM-BC2 failed (exit {p.returncode}); see the messages above")
        made = sorted(f.name for f in stage.iterdir() if f.is_file() and f.name != ".owned_by_unclassified_meta")
        imp_made = []
        if getattr(args, "imputations", False):
            imp_made = _fit_imputations(args, cmd, cp, list(imp), stage, outdir)
            made.append("mi_results.tsv")
        refuse_overlap([outdir / f for f in made], [*inputs, args.counts, args.features, args.metadata], "da")
        manifest.remove_stale(outdir, previous, imp_made + made)
        if imp_made:
            (outdir / "imputations").mkdir(exist_ok=True)
        for f in made + imp_made:
            os.replace(stage / f, outdir / f)
        # outputs of an earlier run that this run did not write again (e.g. structural_zeros.tsv) are removed
        stale = set(previous) | {f.name for pat in DA_OUTPUT_PATTERNS for f in outdir.glob(pat)}
        stale = {n for n in stale if not manifest.STALE_NAME.fullmatch(n)}    # those go only through remove_stale
        for name in sorted(stale - set(made)):
            f = outdir / name
            if f.parent == outdir and f.is_file() and not f.is_symlink():
                f.unlink()
        manifest.seal(outdir, made + imp_made, {"step": "da", "counts": str(Path(args.counts).resolve()),
                                                "inputs_sha256": {k: manifest.sha256(v) for k, v in cp.items()},
                                                "imputations": len(imp_made)})
    finally:
        remove_staging_dir(stage)
    return outdir / "ancombc2_results.tsv"
