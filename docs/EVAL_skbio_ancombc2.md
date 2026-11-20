# scikit-bio ANCOM-BC2 as a backend for `da`: evaluation on OSD18

Date: 2026-11-19. Branch: `worktree-Ancom_pararrell`, fast-forwarded to `em-imputation` (`e190ef0`, includes the
MI pooling fix `537f358`). The R backend (`src/unclassified_meta/r/ancombc2.R`, ANCOMBC 2.14.0, R 4.5.3) is the
baseline throughout. Scripts: `docs/eval/skbio_ancombc2/`.

## Verdict

**Do not replace R, and do not add scikit-bio as a backend for OSD18 now.** scikit-bio 0.7.4 `ancombc2`
reproduces R to about 1e-8, with identical calls, on feature tables where no feature's observed samples confound
condition with plant. On OSD18 that holds for the 1,592 genera with at most two empty condition × plant cells.
It fails for the full table: 1,489 of the 3,081 tested genera are observed in three or fewer of the six cells,
and for those whose observed samples make condition and plant collinear, scikit-bio splits the effect
differently from R's `lm`. Those initial coefficients
feed the bias EM that all features share, so every coefficient shifts. The shift changes which genera are
significant, both in single fits and after Rubin pooling. Pooled over the 10 imputations, scikit-bio adds
*Ampelomyces* as a significant and robust call (q 0.035; R q 0.059).

Single-threaded on the same node, a genus fit with the sensitivity analysis takes 13 s against R's 31 min
(about 145× less CPU). A species fit (12,900 features) takes 31 s. Peak memory is 150–180 MB against R's 0.5 GB.
The structural-zero test is identical to R's. The speed is real, but on this design it does not give
the same statistical result.

## Setup

- Data: OSD18 RNA, 18 samples (healthy/macho × plants P257, P262, P353 × 3 reps), genus-level EM tables
  `taxa_G_em/merged` (expected counts and sealed imputations 001–010), species-level expected counts for the
  benchmark. Local copies matched the sealed sha256 values in the `da_G_em` manifest.
- Model as in production: `--group condition --reference-level healthy --formula 'condition + plant'`,
  `prv_cut 0.10`, `lib_cut 0`, `struc_zero` off, `p_adj_method holm`, `alpha 0.05`, `pseudo 0`,
  `pseudo_sens TRUE`, `s0_perc 0.05`, no random effects.
- R results: the sealed production fits in `OSD18/RNA/da_G_em` (expected + 10 imputations).
- Python: scikit-bio 0.7.4 (numpy 2.5.5 locally, 2.4.6 on Pronghorn), wrapped by
  `docs/eval/skbio_ancombc2/skbio_fit.py`. The wrapper copies the R wrapper's preprocessing (prevalence filter
  ≥ `prv_cut`, factor levels, reference level), runs the sensitivity analysis R's way (pseudo-counts 0, 0.1,
  0.5, 1; `passed_ss` when the q > alpha decision agrees across all four) and writes R's column layout.
  `--r-na` applies R's handling of untestable features (below).

## 1. Feature compatibility

| Feature | R (`ancombc2.R` → ANCOMBC 2.14.0) | scikit-bio 0.7.4 `ancombc2` | OSD18 needs it? |
|---|---|---|---|
| Fixed-effects formula, factor + covariate | yes | yes (patsy formula) | yes, `condition + plant` |
| Reference level | `relevel` | category order of the metadata column | yes |
| Random effects (`rand_formula`, lme4) | yes (`--rand-formula`) | **no** (random intercepts only in `dirmult_lme`, not in `ancombc2`) | no; the OSD18 runs use none, but the pipeline offers `--rand-formula` |
| Prevalence / library-size filters (`prv_cut`, `lib_cut`) | built in | **no**; the caller must filter | yes (`prv_cut 0.10`) |
| Zeros with `pseudo = 0` | log of 0 → NA, feature-specific fit on observed samples | same (masked values) | yes |
| Pseudo-count sensitivity (`pseudo_sens`) | built in (0.1, 0.5, 1 added before bias correction, ANCOMBC ≥ 2.10.1) | **not built in**; repeat the call per pseudo-count (documented recipe, same procedure) | yes (`diff_robust`) |
| Variance regularization (`s0_perc`) | 0.05 quantile, `quantile(type 7)` | `var_quantile`, `nanquantile` linear (same) | yes |
| Multiple testing | `p.adjust` incl. p = 1 entries | statsmodels `multipletests`; **NaN p left out** | yes |
| Untestable feature (dof NA or 0) | `p = 1` (`p_hat[is.na(p_hat)] = 1`) | **p = NaN**, q NaN, W finite | yes, 705 genera |
| Rank-deficient feature design | `lm` aliases the last column (NA) | `match_r=True` rebasing / estimability mask, **differs when condition and plant are collinear in the observed samples** | yes, 1,489 of 3,081 tested genera are observed in ≤ 3 of 6 cells |
| Structural zeros | inside `ancombc2(struc_zero=TRUE)` | standalone `struc_zero()` | off by default (`--struc-zero`) |
| Global / pairwise / Dunnett / trend tests | built in | result-object methods | no (two groups) |
| Labels, identity class, plots, session info | from `ancombc2.R` | none; the wrapper would have to add them | yes for output |
| Threads | `n_cl` (doParallel, `%dorng%`) | single process, BLAS threads | — |

## 2. Numerical comparison (genus, `conditionmacho`)

### Where the two agree and where they do not

| Subset of the expected-count table (same prv_cut-filtered features) | features | max \|Δ lfc\| | max \|Δ q\| | calls (diff, robust) |
|---|---|---|---|---|
| zero-free features | 184 | 6.3e-7 | 2.6e-8 | identical |
| features with zeros, every condition × plant cell observed | 714 | 3.0e-10 | 2.3e-8 | identical |
| features with at most 2 of the 6 condition × plant cells empty (incl. 121 untestable) | 1,592 | 1.8e-8 | 6.9e-8 | identical |
| all tested features | 3,081 | 4.24 | 0.46 | differ (below) |

Fitted on the whole table, the zero-free features differ from R by one exact constant per term (SD of the
difference ~1e-15; `conditionmacho` −0.1376, intercept −0.0286, `plantP262` −0.0087, `plantP353` −0.0508). Their
standard errors differ by about 1 %. That is the signature of a different sampling-bias estimate: the regression
of each feature is the same, but the bias correction applied to all of them is not.

### Root cause

The initial (pre-correction) estimates were exported from both implementations on the same table
(`init_r.R`, `init_py.py`, Pronghorn job 6282387):

- Initial sample offsets (theta): identical (max difference 4e-15).
- The same features enter the bias EM in both (R drops a coefficient with NA; scikit-bio masks it as not
  estimable; the counts agree in every class).
- Initial `conditionmacho` coefficients: identical for all 714 full-rank features (max difference 1.8e-14). They
  differ for 559 features whose observed samples make the design collinear. One example: a genus observed in
  one healthy-P353 and one macho-P262 sample. R's `lm` aliases `plantP353` (NA) and puts the whole difference on
  `conditionmacho` (−6.00). scikit-bio marks every coefficient estimable and splits it (`conditionmacho` −1.87,
  `plantP353` 4.12).
- Bias EM result for `conditionmacho`: δ_em 0.168 in R and 0.028 in scikit-bio. The difference, 0.140, is the
  LFC shift seen in the final table.

So scikit-bio's `match_r=True` mode does not reproduce R when a factor and a covariate are confounded within a
sparse feature's observed samples. This applies to any `condition + plant` model on sparse data, not only to
OSD18. (Not reported upstream; that needs your authorization.)

### Per-fit summary

`summary.py`; R = sealed production fit, Py = scikit-bio with `--r-na`. Calls given as R/Py/both.

| fit | tested | LFC Pearson | median \|Δ lfc\| | max \|Δ lfc\| | max rel. Δ se | untestable | Holm diff | Holm robust | BH q < 0.05 |
|---|---|---|---|---|---|---|---|---|---|
| expected counts | 3,081 | 0.950 | 0.139 | 4.24 | 2.2 % | 705/705/705 | 5/5/4 | 2/2/2 | 11/12/11 |
| imputation 001 | 3,079 | 0.949 | 0.054 | 3.68 | 3.6 % | 698/698/698 | 5/5/5 | 2/2/2 | 15/14/14 |
| imputation 002 | 3,077 | 0.951 | 0.043 | 2.38 | 3.6 % | 704/704/704 | 5/5/5 | 3/3/3 | 14/13/13 |
| imputation 003 | 3,093 | 0.951 | 0.134 | 3.76 | 3.5 % | 700/700/700 | 6/5/5 | 3/3/3 | 14/14/14 |

Expected-count fit, calls that differ (Holm):

| genus | lfc R | lfc Py | q R | q Py |
|---|---|---|---|---|
| Carlavirus | −1.94 | −1.80 | **0.026** | 0.054 |
| Hyaloscypha | 2.58 | 2.72 | 0.062 | **0.040** |

The robust calls (Tobamovirus, Ampelomyces) agree in the expected fit. Wald statistics follow the LFC (same
shift divided by nearly the same se). The covariate terms show the same pattern; R also reports NA estimates for
105 (`plantP262`) and 454 (`plantP353`) aliased coefficients, where scikit-bio reports a finite estimate with
p = 1.

## 3. Structural zeros

R's `ANCOMBC:::.get_struc_zero` (the function `ancombc2(struc_zero = TRUE)` calls) against scikit-bio
`struc_zero()` on all 4,849 genus features, both groups, `neg_lb` FALSE and TRUE, on the expected counts and
imputations 001–003: **identical in every cell** (0 disagreements of 4 × 2 × 2 × 4,849). The pipeline runs with
structural zeros off by default, so this does not affect the OSD18 results.

## 4. Non-testable features (R `p = 1, W ≠ 0`)

ANCOM-BC2 gives p = 1 with a finite W when a feature's residual degrees of freedom are NA (or 0). The MI pooling
fix `537f358` treats such a fit as failed.

- scikit-bio computes the degrees of freedom for these features as NaN and returns **p = NaN, q = NaN** (W finite,
  `Signif` False). The set of features is the same: 705 of 705 in the expected fit, and the same match in
  imputations 001–003.
- statsmodels' adjustment then leaves the NaN entries out. That shrinks the number of tests, so the q-values of
  the other features fall (up to 0.20 for the intercept term). R counts them as p = 1. Adjusting with p = 1, as
  `--r-na` does, restores R's adjustment.
- With NaN, the wrapper's sensitivity analysis must not treat NaN q as a decision. Without `--r-na`, 705
  `passed_ss` calls became NA instead of TRUE.

## 5. Rubin pooling (`da.combine_imputations`, fix `537f358`)

- Pooling R's 10 sealed fits here reproduces the sealed `OSD18/RNA/repool_fix/mi_results_G.tsv` exactly, which
  checks the harness.
- Testability survives the backend switch. `n_fits` agrees for all 3,486 pooled features (1,656 pooled from all 10
  fits in both), with scikit-bio fits passed either raw (p NaN; `combine_imputations` already counts a non-finite
  p as a failed fit) or with `--r-na` (p = 1 and W ≠ 0, the `537f358` rule). The fix is therefore not bypassed.
  A pass-through needs either convention to hold, plus `label` and `identity_class` columns, which the wrapper
  does not write yet.
- The pooled statistics do not survive it. The bias shift differs from one imputation to the next (0.04–0.14),
  so it enters the between-imputation variance.

| pooled, 10 imputations, `conditionmacho` | max \|Δ\| | median \|Δ\| |
|---|---|---|
| lfc | 1.36 | 0.070 |
| se | 0.65 | 0.018 |
| Barnard–Rubin df | 2.5 | 0.17 |
| FMI | 0.64 | 0.052 |
| p | 0.52 | 0.059 |
| q (Holm) | 0.51 | 0 |

| genus | lfc R / Py | FMI R / Py | df R / Py | q R / Py | call R / Py |
|---|---|---|---|---|---|
| Tobamovirus | −8.09 / −8.02 | 0.069 / 0.0004 | 11.4 / 12.3 | <1e-4 / <1e-4 | robust / robust |
| Ampelomyces | 4.45 / 4.52 | 0.082 / 0.055 | 9.5 / 9.8 | 0.059 / **0.035** | none / **robust** |
| Ophiobolus | 3.91 / 3.98 | 0.057 / 0.086 | 5.3 / 5.1 | 0.038 / 0.043 | diff / diff |

## 6. Runtime and memory

Pronghorn job 6282331 on one node (cpu-34, Xeon E5-2683 v4), fits one at a time, single-threaded
(`OMP/OPENBLAS/MKL_NUM_THREADS=1`, R `--threads 1`). Each fit includes the four pseudo-count fits, as in
production. scikit-bio runs in the lab `micromamba.sif` container: pip wheels need glibc newer than the el7
nodes have (lab wiki, `guide/pronghorn.md`, glibc section). /usr/bin/time wall time and peak RSS:

| fit | features tested | R wall | R CPU | R peak RSS | scikit-bio wall | scikit-bio CPU | scikit-bio peak RSS | speedup (CPU) |
|---|---|---|---|---|---|---|---|---|
| genus, expected counts | 3,081 | 30:59 | 1,852 s | 495 MB | 13.8 s (cold 33.3 s, 22.0 s) | 12.8 s | 149 MB | 145× |
| genus, imputation 001 | 3,079 | 31:04 | 1,857 s | 498 MB | 13.0 s | 12.9 s | 147 MB | 144× |
| genus, imputations 002 / 003 | 3,077 / 3,093 | not timed; see note | | | 13.3 s / 13.2 s | 13.2 s / 13.2 s | 149 MB | |
| species, expected counts | 12,900 | not timed; see note | | | 31.0 s | 30.6 s | 180 MB | |

scikit-bio times include container start-up and imports. The first two runs were slower while the
container image and inputs were cold on the shared filesystem; CPU time stays at about 12.5 s. The R fit run
single-threaded here is bit-identical to the sealed production fit run with `n_cl 8` (max |Δ| 0 for every
column), so R's thread count does not change its results. At the time of writing, the job was still timing R on
imputations 002–003 and the species table. Its log is
`/data/gpfs/assoc/pgl/develop/unclassified_meta_test/skbio_eval/logs/skbio_vs_r_6282331.log`.

scikit-bio results are reproducible across machines: the Pronghorn run (numpy 2.4.6, container) and the local run
(numpy 2.5.5) agree to 1.4e-10 with identical calls.

For reference, production `da --imputations` at genus level (job 6272796, R, 8 CPUs, `n_cl 8`): 11 fits in
2:00:13 (about 11 min per fit), peak RSS 0.5 GB; the species-level job runs about 45 min per fit.

## 7. Recommendation for OSD18

1. Keep R as the backend and as the reference for every reported OSD18 result.
2. Do not add a scikit-bio backend yet, not even as an option. It would give different calls on exactly the
   sparse `condition + plant` data this pipeline is built for, and the difference would not be visible to the
   user.
3. Revisit when scikit-bio reproduces R on confounded sparse designs. The test is ready: rerun
   `skbio_fit.py --r-na` on `taxa_G_em/merged/counts.tsv` and `summary.py`; equivalence means the full table
   behaves like the zero-free and full-rank subsets do now (≤ 1e-6, identical calls). A minimal reproducer for
   an upstream report is the two-sample feature in the root-cause section.
4. If a fast path is still wanted before then, the only setting where this evaluation supports equivalence is a
   table filtered so that no feature has more than two empty condition × plant cells. That filter does not
   recover R's answer on the current table: the bias estimate then comes from different features. It is a
   different analysis (more filtering), so it is a scientific decision, not a backend swap.
5. Fit-level parallelism of the R backend (`da --imputations --workers`) is a statistically neutral way to cut
   the wall time.
