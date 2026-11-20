# Multiple-imputation pooling audit (`da.combine_imputations`)

Date: 2026-11-19. Code audited: `em-imputation` at e190ef0 (cabe340 EM + MI, 537f358 p = 1 fix).
Data: the ten saved ANCOM-BC2 genus-level fits of OSD18 RNA (`OSD18/RNA/da_G_em/imputations/`), re-pooled without
refitting. Baseline: `OSD18/RNA/repool_fix/mi_results_G.tsv`. The sealed `da_G_em/mi_results.tsv` was written by the
code before 537f358 and is not used as a reference.

## How ANCOM-BC2 marks a test it could not run

ANCOMBC 2.14.0 `.ancombc2_core` (fixed effects, `.iter_mle`):

```r
dof   = summary(lm(y ~ ...))$df[2]           # residual df; 999 if lm fails
W     = beta_hat / se_hat
p_hat = 2 * pt(abs(W), df = dof, lower.tail = FALSE)
p_hat[is.na(p_hat)] = 1
```

Zero counts are missing on the log scale, so `dof` = (samples where the feature is nonzero) − rank(X). With
`condition + plant` (rank 4), a feature nonzero in four or fewer samples has `dof = 0`. `pt(·, 0)` is NaN and the
fill sets p = 1. ANCOM-BC2 does not export `dof`, so the result table never flags these features directly. The flag
can still be read off exactly: a valid two-sided t test gives p = 1 only at W = 0. So p = 1 with a nonzero W is
always the fill. In the ten OSD18 fits all 12,812 rows with p = 1 have W ≠ 0 (smallest |W| is 6.9e-4). None is a
real p = 1.

An exact flag would need the R script to export `dof`, which means refitting. It is not needed for correctness.

## Findings

| # | Item | Verdict |
|---|------|---------|
| 1 | Rubin: Q̄ = mean, W = mean se², B = var(ddof 1), T = W + (1 + 1/M) B | correct |
| 2 | Barnard–Rubin df: ν_old = (M−1)/λ², ν_obs = (ν_com+1)/(ν_com+3) ν_com (1−λ), ν = (1/ν_old + 1/ν_obs)⁻¹; ν = ν_com when B = 0 | correct (checked against R) |
| 3 | ν_com recovered from each fit's own (W, p); median over fits | correct. Since 537f358 it used every fit; now it uses only the fits that count |
| 4 | `fmi_<term>` | is λ = (1+1/M)B/T, not Rubin's γ = (r + 2/(ν+3))/(r+1). Both are valid; the label is now documented |
| 5 | p = 1 with W ≠ 0 treated as failed (537f358) | correct for OSD18. **Bug:** with W = NA, or no W column, the fit counted as valid and was re-tested (synthetic case: p = 2e-9). Fixed |
| 6 | Feature testable in only some fits | p = 1, not pooled. Correct: pooling the testable subset conditions on the imputations that gave the feature enough nonzero samples. In OSD18 that would put 9 of the 320 such features at p < 1e-6 (e.g. Urinicoccus, testable in 2 of 10 fits) |
| 7 | Can an invalid test become significant? | No longer: p is forced to 1 unless every fit counts (`n_fits == M`), and `diff` also requires it |
| 8 | Expected-count fit in the pooling | never: only the `fit_###` results go to `combine_imputations` (also covered by `test_da_imputations_fits_every_sealed_table…`) |
| 9 | Multiple testing | `p_adjust` matches R `p.adjust` for BH and Holm. Unpooled features stay in the family with p = 1, as ANCOM-BC2 does. **The OSD18 run used `--p-adjust holm` (the default, as in ANCOM-BC2), so the baseline q values are Holm (FWER), not BH** |
| 10 | Malformed fit sets | before: 0 fits raised IndexError, 1 fit was silently "pooled", a missing coefficient raised KeyError. Now each is a clear StepError |
| 11 | `diff` threshold | `q < alpha`, while ANCOM-BC2 uses `q <= alpha`. No effect unless q equals alpha exactly; left as is |

## Independent numerical check

A separate implementation (scipy t distribution, brentq df recovery, R-definition `p.adjust`; scratch script, not in
the repository) re-pooled the ten fits. Against the baseline:

- features 3,486; pooled 1,656; `n_fits` identical for every feature
- largest relative difference: lfc 1.1e-13, se 4e-16, df 2.7e-13, fmi 1.1e-14, p and q 1.2e-12
- the new code reproduces the baseline to 4.4e-16 (absolute), with identical `diff` and `diff_robust`

Status of the 3,486 features: pooled 1,656; untestable in every fit 813; not estimated (prevalence filter) 594;
partly untestable 320; incomplete 103.

### The previously reported findings

| Genus | LFC | df | fmi (λ) | γ | q Holm (baseline) | q BH | n_fits | diff_robust |
|-------|-----|----|---------|---|-------------------|------|--------|-------------|
| Tobamovirus | −8.09 | 11.4 | 0.07 | 0.20 | 7.5e-8 | 7.5e-8 | 10 | TRUE |
| Ophiobolus | +3.91 | 5.3 | 0.06 | 0.29 | 0.038 | 0.019 | 10 | FALSE |
| Ampelomyces | +4.45 | 9.5 | 0.08 | 0.23 | 0.059 | 0.020 | 10 | FALSE |
| Carlavirus | −1.80 | 11.3 | 0.08 | 0.21 | 0.25 | 0.04999 | 10 | FALSE |
| Potexvirus | −5.37 | 12.2 | 0.01 | – | 0.24 | 0.04999 | 10 | FALSE |
| Parastagonospora | +3.44 | 2.7 | 0.66 | 0.78 | 1 | 1 | 10 | FALSE |
| Neocucurbitaria | +2.05 | 0.9 | 0.79 | 0.90 | 1 | 1 | 10 | FALSE |

All the reported numbers hold, but the q values are Holm. Under BH, five genera have q < 0.05: Ampelomyces at
0.020, and Carlavirus and Potexvirus right at 0.04999. Parastagonospora and Neocucurbitaria are pooled in all ten
fits, but most of their variance comes from read assignment, and Neocucurbitaria's reference df is below 1. Their
p = 1 after adjustment reflects real assignment uncertainty, not a failed test. Ten imputations are too few at these
fmi values (the README advises ≥ 20 above 0.3).

## Interpretation of non-testable imputations

A feature that ANCOM-BC2 cannot test in some imputations is one whose presence in enough samples depends on how
ambiguous reads are assigned. Its evidence comes from a handful of samples that receive reads in some draws and not
in others. Pooling only the draws in which it was testable conditions on the outcome of the imputation. Those draws
put more reads into the feature, so the subset overstates both its abundance and its precision. Reporting p = 1, with
`status = partly_untestable` and the counts `n_fits` / `n_untestable`, follows ANCOM-BC2's own treatment of a test it
cannot run and keeps the feature in the multiple-testing family. These features need more samples or a coarser rank,
not a different pooling rule.

## Changes

- `src/unclassified_meta/da.py`
  - `_untestable()`: p ≥ 1 and |W| > 1e-12 marks the fill; W comes from lfc/se when the W column is NA or missing
  - `combine_imputations`: refuses fewer than 2 fits and fits whose coefficient columns differ; recovers ν_com only
    from counted fits; forces p = 1 unless `n_fits == M`; new columns `n_untestable_<term>` and `status_<term>`
  - `mi_results.tsv` header names the status values and defines fmi as λ
- `tests/test_mi_pooling.py`: 19 regression tests (non-testable and valid p = 1, missing W, missing p,
  partly untestable, zero/inf/NaN/negative se, B = 0, absent features, missing/extra coefficients, < 2 fits,
  Rubin/Barnard–Rubin and BH/Holm against R, BH family, internal consistency)
- `README.md`: MI section describes testability, status and λ

## Decision: BH (2026-11-20)

`da --p-adjust` now defaults to BH in the CLI and in `ancombc2.R`. Both the expected-count fit and the pooled
imputations are adjusted for the false discovery rate. Holm, ANCOM-BC2's own default, is still available.

The ten genus fits re-pooled with BH are in `OSD18/RNA/repool_fix/mi_results_G_BH.tsv`. The Holm file
`mi_results_G.tsv` is kept unchanged. BH q < 0.05: Tobamovirus (7.5e-8, robust), Ophiobolus (0.019),
Ampelomyces (0.020, robust), Potexvirus (0.04999), Carlavirus (0.04999). Two are robust.

Caveat: `diff_robust` needs ANCOM-BC2's own `diff_robust` in at least half of the fits. The saved fits ran with Holm,
so those per-fit calls are Holm calls. Holm is stricter than BH, so the robust calls here are, if anything,
understated. A rerun with BH would align them.
