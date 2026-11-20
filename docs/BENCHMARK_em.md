# taxa-table EM: profile, optimisation and validation (2026-11-19)

Baseline: `em-imputation` at `e190ef0`. Optimised: this change (`src/unclassified_meta/em.py`).
Data: OSD18 RNA, Centrifuger `-k 5` per-read output against `cfr_nr`, `--min-score 64`, rank S,
`--drop-lineage Viridiplantae,Metazoa`, `--assign em --imputations 10 --seed 11`. Machine: Pronghorn, one core per
job (`OMP/OPENBLAS/MKL_NUM_THREADS=1`), numpy 1.26.4. Scripts and commands: [bench_em/](bench_em/) (see the end).

## 1. Where the time went (baseline)

Each sample's equivalence classes were taken from the real `taxa-table` parse (`dump_classes.py`) and EM and 10
imputations were run on them alone, with timers around every stage (`run_em.py`). Seconds, sum over 18 samples:

| stage | baseline (s) | share of the 7,873 s measured |
|---|---|---|
| input parsing (Centrifuger per-read, taxonomy lookup, class building) | 274 | 3 % |
| EM (SQUAREM: E-steps 1,717, log-likelihoods 1,205, gaps 132) | 3,299 | 42 % |
| imputation: Gibbs assignment draws | 2,540 | 32 % |
| imputation: SIR chain starts (Dirichlet-tree draws) | 1,435 | 18 % |
| laminar groups (closed-form EM and exact draws), setup | 7 | 0 % |
| taxonomy, writing, merge with the panel | ~1 min (cProfile of the optimised end-to-end run) | |

The work is not the size of the data: the reads are already collapsed into equivalence classes, so a sample has
at most 35,231 ambiguous classes and 28,605 targets. Two things made it slow:

1. **The EM stop rule.** EM stopped when the tangent-plane bound max_i g_i - N (g the gradient, N the reads of the
   group) certified the log-likelihood to be within 0.5 of its maximum. A target with a tiny share whose
   maximum-likelihood share is a little larger keeps that bound open: EM multiplies its share by g_i / N per
   iteration, and N counts every read of the group (16 M in sample 20). In sample 20 one target
   (`#unresolved|72763`, 0.9 expected reads, no reads of its own) held the bound open for 147,634 iterations while
   the log-likelihood did not change in the fourth decimal after iteration 20,000.
2. **Python loops in the imputation.** Each chain start drew 256 candidate proportion vectors, each by a Python
   recursion over the laminar tree (578 ms per draw in sample 25); each of 500 Gibbs sweeps per chain drew every
   class by conditional binomials and `np.add.at` (98 ms per sweep in sample 25), one chain and one crossing group
   at a time.

Centrifuger's quantifier (`Quantifier.hpp`, `EstimateAbundanceWithEM`/`EMupdate`) was read for comparison. It
already works on weighted equivalence classes (`_readAssignment`: target set + summed weight) with integer target
indices, runs plain EM for at most 1,000 iterations with an L1 change criterion, and then redistributes abundance
down the taxonomy tree with genome-length correction. Equivalence classes and integer-indexed arrays were already
used here. The capped iteration count and L1 criterion were not adopted (they are the kind of uncertified stop
rule that QUARANTINE.md rejects), nor were the genome-length correction and the tree redistribution (this table
counts reads, not cells, at one rank).

## 2. What changed

| change | effect on results |
|---|---|
| class proportions computed once per point and shared by the E-step, log-likelihood and gap | none: with `dual_tol=None` the iterates are bitwise those of the baseline (checked on samples 23, 34) |
| stop when the group's log-likelihood is certified within 1e-4 of its maximum, by the smaller of the tangent-plane bound and a Lagrangian-dual bound (`_dual_gap`) | the same estimator (the MLE) computed more accurately: certified within 1e-4 per group instead of 0.5 |
| SIR chain starts: all 256 Dirichlet-tree draws at once as gamma variates normalised per node (`_TreeDraw`) | same distribution (Dirichlet = normalised gammas) |
| all crossing groups and all M chains advanced together; proportions drawn as one gamma vector | same distribution: only ratios within a group's classes enter the next assignment, and independent gammas normalised within a group are that group's Dirichlet |
| reads of a class drawn as one multinomial: read by read (one uniform each) for classes of up to 16 reads, conditional binomials above | same distribution (a multinomial is the sum of independent categorical draws); a candidate of weight zero is never drawn |

The dual bound: every term of the log-likelihood is n_j log(y_j) with y_j = a_j . theta the share of class j.
From log x <= x - 1, for any weights w_j > 0, max over the simplex of l <= sum_j n_j log(n_j / w_j) - N +
max_i (A^T w)_i. With w_j = n_j / y_j(theta) this is the tangent-plane bound. Lowering the weights of the terms
of the few targets above a level lam lowers their (A^T w)_i to lam at a cost of their own reads times
-log(1 - eps); the bound is evaluated exactly for those weights, so it remains a certificate. Tested in
`tests/test_em_fast.py` (never below the true gap; never above the tangent bound).

Not used: Numba (not in the environment; NumPy was enough), sparse matrices (no SciPy; `bincount` on flat index
arrays is the sparse product), class-level weights for "same targets, different likelihood" (Centrifuger's per-read
output gives every accepted candidate of a read the same weight, so a read's information is its candidate set;
reads with the same set are already one class).

## 3. Per-sample results (EM + 10 imputations, seconds; baseline / optimised)

Classified reads are those entering the EM step. Iterations: the largest over the sample's groups. ll: total
log-likelihood of the optimised expected counts minus the baseline's (higher is closer to the maximum). max diff:
largest change of an expected count of a taxon at the rank, in reads.

| sample | classified reads | targets | ambiguous reads | ambiguous classes | candidates per ambiguous read | crossing groups | iterations | EM s | imputation s | speed-up | RSS MB | ll | max diff |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 19 | 13,621,767 | 9,655 | 1,215,706 | 10,771 | 2.74 | 66 | 2,259 / 459 | 18 / 3 | 173 / 18 | 9.2 | 120 / 110 | +0.53 | 1.65 |
| 20 | 15,937,557 | 7,635 | 1,263,482 | 8,334 | 2.62 | 45 | 147,634 / 326 | 762 / 2 | 124 / 15 | 52.1 | 91 / 108 | +0.41 | 0.98 |
| 21 | 15,567,814 | 6,085 | 1,155,049 | 6,692 | 2.62 | 41 | 3,028 / 679 | 31 / 4 | 154 / 13 | 11.0 | 81 / 100 | +1.30 | 1.00 |
| 22 | 3,997,678 | 7,627 | 91,678 | 9,963 | 2.47 | 47 | 6,294 / 949 | 80 / 12 | 197 / 22 | 8.0 | 93 / 98 | +0.31 | 0.68 |
| 23 | 1,873,837 | 7,793 | 66,082 | 9,710 | 2.59 | 43 | 936 / 282 | 8 / 3 | 135 / 21 | 6.0 | 92 / 98 | +0.17 | 0.33 |
| 24 | 2,647,380 | 14,297 | 100,080 | 17,962 | 2.73 | 96 | 1,020 / 487 | 21 / 6 | 372 / 35 | 9.5 | 150 / 118 | +0.93 | 1.15 |
| 25 | 1,597,499 | 28,605 | 465,741 | 35,231 | 2.37 | 168 | 25,744 / 2,378 | 1191 / 76 | 732 / 73 | 12.9 | 262 / 168 | +1.88 | 1.48 |
| 26 | 1,021,221 | 10,020 | 304,953 | 11,784 | 2.31 | 50 | 2,210 / 843 | 21 / 7 | 162 / 19 | 7.0 | 124 / 109 | +0.54 | 0.98 |
| 27 | 1,554,056 | 8,465 | 448,044 | 9,897 | 2.33 | 61 | 603 / 492 | 5 / 6 | 155 / 22 | 5.8 | 108 / 104 | +0.23 | 0.33 |
| 28 | 109,579 | 11,452 | 46,618 | 13,552 | 3.04 | 62 | 3,974 / 1,906 | 84 / 39 | 263 / 31 | 5.0 | 134 / 113 | +0.44 | 0.43 |
| 29 | 278,143 | 16,072 | 139,300 | 26,638 | 3.10 | 73 | 8,765 / 3,864 | 208 / 76 | 276 / 55 | 3.7 | 195 / 133 | +0.47 | 0.92 |
| 30 | 345,091 | 14,520 | 127,340 | 29,495 | 3.19 | 42 | 11,590 / 2,288 | 458 / 82 | 340 / 67 | 5.4 | 193 / 129 | +0.23 | 0.41 |
| 31 | 284,934 | 11,159 | 47,911 | 14,140 | 3.08 | 69 | 788 / 490 | 11 / 4 | 260 / 26 | 9.1 | 133 / 108 | +0.64 | 0.54 |
| 32 | 218,573 | 10,591 | 33,608 | 11,386 | 3.08 | 71 | 594 / 351 | 6 / 3 | 173 / 20 | 7.6 | 124 / 114 | +1.07 | 1.00 |
| 33 | 192,183 | 10,194 | 41,195 | 13,029 | 3.09 | 46 | 365 / 439 | 4 / 8 | 161 / 26 | 4.9 | 125 / 105 | +0.66 | 0.67 |
| 34 | 119,557 | 6,639 | 65,694 | 6,350 | 2.37 | 54 | 3,182 / 1,782 | 24 / 18 | 129 / 13 | 4.8 | 82 / 111 | +0.26 | 0.34 |
| 35 | 202,234 | 9,471 | 110,018 | 10,875 | 2.47 | 75 | 2,490 / 652 | 23 / 8 | 189 / 24 | 6.6 | 111 / 107 | +0.37 | 0.40 |
| 36 | 429,495 | 15,815 | 238,038 | 22,919 | 2.52 | 94 | 14,810 / 2,552 | 345 / 74 | 304 / 51 | 5.2 | 189 / 127 | +1.46 | 1.00 |
| **total** | 59,998,598 | | 5,960,537 | | | | 236,286 / 21,219 | 3299 / 432 | 4300 / 551 | 7.7 | | | |

## 4. End-to-end `taxa-table` (the production command)

`/usr/bin/time`, one core, same CPU model (Xeon E5-2683 v4) for each pair. small: samples 28, 32, 34;
medium: 22, 25, 29, 30, 33, 36 (the two slowest EM samples included); complete: all 18 with `--merge-panel` and
`--panel-taxonomy` as in stage A.

| input | samples | input fragments | baseline s | optimised s | speed-up | CPU (user s) base / opt | peak RSS MB base / opt |
|---|---|---|---|---|---|---|---|
| small | 3 | 1,149,639 | 844 | 138 | 6.1x | 837 / 136 | 479 / 495 |
| medium | 6 | 10,509,161 | 5,166 | 803 | 6.4x | 5,156 / 770 | 906 / 906 |
| complete (+ merge-panel, panel-taxonomy) | 18 | 74,885,745 | 9,855 | 1,607 | 6.1x | 9,831 / 1,564 | 2,620 / 2,620 |

CPU use was 98-99 % of one core in every run. The medium optimised run on a Xeon Gold 6130 node took 593 s; the
table uses the run on the same CPU model as its baseline. In the production log the complete baseline took
2 h 41 min (162 min); here 2 h 44 min. The remaining time of the optimised run (cProfile): reading the Centrifuger
output and building the classes (the same as `--assign lca`, whose whole run takes 5 min), EM, imputation and
about 1 min for the taxonomy, writing and merge. Peak memory is the per-read parse and did not change.

Expected counts in the complete run: the largest change of any cell is 1.65 reads, 43 of 592,974 cells change
by more than 0.5 read, and every sample's total is unchanged (1e-5).

## 5. Validation

**Expected counts** (`cmp_counts.py`, all 18 samples). Read mass is conserved (sum of expected counts minus reads:
at most 2e-9). The optimised fits have a higher log-likelihood than the baseline in every sample (+0.17 to +1.88
nats): the baseline's 0.5-nat certificate applied to each of a sample's 41-168 crossing groups separately, so it
could leave a sample several nats below the maximum, while the optimised fits are certified to 1e-4 per group. The
largest change of a taxon's expected count is 1.65 reads (samples 19, 25: 1.65 and 1.48); for taxa with at least
10 expected reads, the largest relative change is 5 %. These changes move the estimates towards the maximum, not
away from it. Tolerance used for the tests: the certificate itself (log-likelihood within 1e-4 of an
independently computed maximum, `tests/test_em_fast.py`), because a count difference has no fixed scale in
directions the data barely identify.

**Seeded outputs** are not identical to the baseline's: the vectorised draws consume random numbers in a different
order. The same seed still gives the same tables on every run (test). Exact equivalence holds only for the EM
iterates with `dual_tol=None`.

**Imputation distribution** (`imp_draws.py`, `impval.py`). 200 independent imputations per sample with the
baseline (seed 101) and with the optimised code (seeds 202 and 303). For every target with a mean of at least one
read, z = difference of means / its standard error; equal distributions give z ~ N(0, 1). The second optimised
seed gives the reference for two runs of one sampler.

| sample | reads | M | time base / opt (s) | comparison | targets | mean z^2 | \|z\| > 3 | max \|z\| | SD ratio median (1-99 %) | mean log SD ratio (se) |
|---|---|---|---|---|---|---|---|---|---|---|
| 23 | 1,873,837 | 200 | 3,120 / 473 | base vs opt | 1,879 | 1.03 | 0.27 % | 3.50 | 1.00 (0.72-1.40) | -0.0002 (0.0026) |
| | | | | opt vs opt | 1,884 | 1.02 | 0.16 % | 3.22 | 1.00 (0.72-1.40) | -0.0015 (0.0025) |
| 31 | 284,934 | 200 | 5,201 / 686 | base vs opt | 2,692 | 1.02 | 0.37 % | 4.48 | 1.00 (0.75-1.41) | +0.0003 (0.0021) |
| | | | | opt vs opt | 2,681 | 0.99 | 0.11 % | 3.28 | 1.00 (0.73-1.41) | +0.0002 (0.0022) |
| 34 | 119,557 | 200 | 3,265 / 256 | base vs opt | 1,534 | 1.02 | 0.13 % | 3.56 | 1.00 (0.73-1.37) | +0.0000 (0.0026) |
| | | | | opt vs opt | 1,540 | 1.00 | 0.26 % | 4.27 | 1.00 (0.71-1.30) | -0.0038 (0.0027) |

Equal distributions give mean z^2 = 1, |z| > 3 in 0.27 % of targets and a mean log SD ratio of 0.

Every draw sums to the sample's reads. In sample 31 the largest |z| (4.48, `#unresolved|323620`) is between the
baseline and the seed-202 run; the baseline agrees with the seed-303 run there (means 1.9, 2.5, 2.1 reads): the
spread of a discrete, skewed count of a few reads, not a difference of samplers.

**Output compatibility** (`compat.py`, small and medium, and complete with `merged/`): manifests verify; the same
files, feature-table columns, ledger columns, seed and number of imputations; every imputation table holds whole
reads with the rows of `counts.tsv` (what `da --imputations` requires); every imputation ledger balances (input =
unclassified + below threshold + dropped + kept + unresolved) and matches its table. Features present in one run
only have at most 9e-5 expected reads (rows kept because some imputation drew them).

**Tests**: `pytest` 266 passed (full suite with R, 2026-11-19), including `tests/test_em_fast.py`: dual bound
validity, MLE within the certificate on crossing problems, dominated (zero-share) and dominant taxa, class-level
draws against the earlier conditional-binomial draw and the multinomial moments (1, 7, 16, 17 and 500 reads per
class), zero-weight and all-zero classes, batched Dirichlet-tree draws against the per-node recursion, mass
conservation, seed reproducibility, and independent groups drawn jointly against the same group drawn alone.

## 6. Reproduce

On Pronghorn (paths in the scripts; `B=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/em_opt`, with the
baseline source in `$B/baseline_src` (`git archive e190ef0 src`) and this change in `$B/opt3_src`):

```bash
sbatch dump.sbatch RNA S                                     # classes per sample (the real taxa-table parse)
sbatch --array=1-18 run_em.sbatch baseline_src RNA_S base_RNA_S
sbatch --array=1-18 run_em.sbatch opt3_src RNA_S opt3_RNA_S
python agg.py results/base_RNA_S; python agg.py results/opt3_RNA_S          # stage tables
PYTHONPATH=opt3_src python cmp_counts.py results/base_RNA_S results/opt3_RNA_S classes
for set in small medium complete; do for code in baseline_src opt3_src; do sbatch e2e.sbatch $code $set; done; done
PYTHONPATH=opt3_src python compat.py e2e/baseline_src_complete_S e2e/opt3_src_complete_S
for s in RNA_S_23_RNA_P262_H_rep2 RNA_S_31_RNA_P262_M_rep1 RNA_S_34_RNA_P353_M_rep1; do
  sbatch imp_draws.sbatch baseline_src $s 200 101; sbatch imp_draws.sbatch opt3_src $s 200 202
  sbatch imp_draws.sbatch opt3_src $s 200 303; done
python impval.py impval RNA_S_23_RNA_P262_H_rep2 RNA_S_31_RNA_P262_M_rep1 RNA_S_34_RNA_P353_M_rep1
```

Not benchmarked here: rank G and the DNA run (the same code paths); `panel-quant --ambiguous em` uses the same
functions.
