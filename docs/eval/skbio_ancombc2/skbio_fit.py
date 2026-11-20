"""Run scikit-bio ANCOM-BC2 with the settings of the project's R wrapper (r/ancombc2.R) and write the result in
the R wrapper's column layout (lfc_/se_/W_/p_/q_/diff_/passed_ss_/diff_robust_), so that it can be compared with
the R output and passed to da.combine_imputations.

R settings reproduced: tax_level NULL, prv_cut (features with nonzero share >= prv_cut), lib_cut 0, pseudo 0,
pseudo_sens TRUE (pseudo-counts 0.1, 0.5, 1; passed_ss when the q > alpha decision agrees in all four fits),
s0_perc 0.05 (var_quantile), struc_zero FALSE, p_adj_method holm (or as given), alpha.

usage: skbio_fit.py COUNTS META OUT.tsv [--formula F] [--group G] [--reference-level L] [--prv-cut X]
                    [--p-adjust holm] [--alpha 0.05] [--no-sens] [--timing OUT.json]
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")

import numpy as np
import pandas as pd


def r_name(cov: str) -> str:
    """patsy 'condition[T.macho]' -> R 'conditionmacho'; 'Intercept' -> '(Intercept)'."""
    if cov == "Intercept":
        return "(Intercept)"
    if "[T." in cov:
        var, lev = cov.split("[T.", 1)
        return var + lev.rstrip("]")
    return cov


def load(counts_path: str, meta_path: str, group: str, ref: str | None, formula: str):
    counts = pd.read_csv(counts_path, sep="\t", dtype={"feature_id": str}).set_index("feature_id")
    meta = pd.read_csv(meta_path, sep="\t", dtype=str, keep_default_na=False).set_index("sample")
    meta = meta.loc[[s for s in meta.index if s in counts.columns]]
    counts = counts[meta.index]
    vars_ = [v.strip() for v in formula.replace("*", "+").replace(":", "+").split("+")]
    for v in vars_:
        num = pd.to_numeric(meta[v], errors="coerce")
        if v != group and num.notna().all():
            meta[v] = num
        else:
            levels = sorted(meta[v].unique())          # R factor(): alphabetical levels
            if v == group and ref:
                levels = [ref] + [x for x in levels if x != ref]
            meta[v] = pd.Categorical(meta[v], categories=levels)
    return counts, meta


SM_METHOD = {"holm": "holm", "bh": "fdr_bh", "fdr": "fdr_bh", "by": "fdr_by", "bonferroni": "bonferroni"}


def r_na_convention(res: pd.DataFrame, p_adjust: str, alpha: float) -> pd.DataFrame:
    """R's .ancombc2_core: p_hat[is.na(p_hat)] = 1, then p.adjust over every feature (the p = 1 entries count)."""
    from statsmodels.stats.multitest import multipletests
    res = res.copy()
    res["pvalue"] = res["pvalue"].fillna(1.0)
    for cov in dict.fromkeys(res.index.get_level_values("Covariate")):
        idx = res.index.get_level_values("Covariate") == cov
        res.loc[idx, "qvalue"] = multipletests(res.loc[idx, "pvalue"].to_numpy(), method=SM_METHOD[p_adjust.lower()])[1]
    res["Signif"] = res["qvalue"] <= alpha
    return res


def fit(counts: pd.DataFrame, meta: pd.DataFrame, formula: str, prv_cut: float, p_adjust: str, alpha: float,
        sens: bool = True, r_na: bool = False) -> pd.DataFrame:
    from skbio.stats.composition import ancombc2
    prevalence = (counts != 0).sum(axis=1) / counts.shape[1]
    kept = counts.loc[prevalence >= prv_cut]
    table = kept.T                                      # samples x features
    pseudos = (0, 0.1, 0.5, 1) if sens else (0,)
    fits = [ancombc2(table, meta, formula, pseudocount=p, var_quantile=0.05, p_adjust=p_adjust, alpha=alpha)
            for p in pseudos]
    if r_na:
        fits = [type("R", (), {"result": r_na_convention(f.result, p_adjust, alpha)}) for f in fits]
    main = fits[0].result
    covs = list(dict.fromkeys(main.index.get_level_values("Covariate")))
    out = pd.DataFrame({"feature_id": kept.index})
    stat = {"lfc": "Log(FC)", "se": "SE", "W": "W", "p": "pvalue", "q": "qvalue"}
    for key, col in stat.items():
        for c in covs:
            out[f"{key}_{r_name(c)}"] = main.xs(c, level="Covariate")[col].reindex(kept.index).to_numpy()
    for c in covs:
        q = np.column_stack([f.result.xs(c, level="Covariate")["qvalue"].reindex(kept.index).to_numpy(dtype=float)
                             for f in fits])
        diff = main.xs(c, level="Covariate")["qvalue"].reindex(kept.index).to_numpy(dtype=float) <= alpha
        diff &= ~np.isnan(main.xs(c, level="Covariate")["qvalue"].reindex(kept.index).to_numpy(dtype=float))
        frac = (q > alpha).sum(axis=1) / len(pseudos)    # R: sum(x > alpha) / length(pseudo_list); NA -> NA
        passed = np.where(np.isnan(q).any(axis=1), np.nan, (frac == 0) | (frac == 1))
        out[f"diff_{r_name(c)}"] = np.where(diff, "TRUE", "FALSE")
        out[f"passed_ss_{r_name(c)}"] = np.where(np.isnan(passed), "NA", np.where(passed == 1, "TRUE", "FALSE"))
        out[f"diff_robust_{r_name(c)}"] = np.where(diff & (passed == 1), "TRUE", "FALSE")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("counts"); ap.add_argument("meta"); ap.add_argument("out")
    ap.add_argument("--formula", default="condition + plant"); ap.add_argument("--group", default="condition")
    ap.add_argument("--reference-level", default="healthy"); ap.add_argument("--prv-cut", type=float, default=0.10)
    ap.add_argument("--p-adjust", default="holm"); ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--no-sens", action="store_true"); ap.add_argument("--timing")
    ap.add_argument("--r-na", action="store_true", help="R convention: p = 1 where dof is NA, adjusted over all")
    a = ap.parse_args(argv)
    t0, c0 = time.perf_counter(), time.process_time()
    counts, meta = load(a.counts, a.meta, a.group, a.reference_level, a.formula)
    t1 = time.perf_counter()
    res = fit(counts, meta, a.formula, a.prv_cut, a.p_adjust, a.alpha, not a.no_sens, a.r_na)
    t2 = time.perf_counter()
    res.to_csv(a.out, sep="\t", index=False)
    import skbio
    info = {"wall_s_total": time.perf_counter() - t0, "wall_s_fit": t2 - t1, "cpu_s": time.process_time() - c0,
            "maxrss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, "features_tested": len(res),
            "skbio": skbio.__version__, "numpy": np.__version__}
    if a.timing:
        with open(a.timing, "w") as fh:
            json.dump(info, fh, indent=1)
    print(json.dumps(info))


if __name__ == "__main__":
    main()
