# Differential abundance with ambiguous species assignments

**EM plus multiple imputation was conservative in this high-depth simulation.** In the 28-dataset primary comparison at cabe340, MI detected 0 of 84 fourfold changes, versus 2 for the expected-count fit and 4 for LCA. Its centered 95% intervals covered 83 of 84 true changes (98.8%). Under the global null, 6 of 560 MI tests had raw p < 0.05 (1.07%), and none had BH q < 0.05. Mean changed-species bias was +0.027 log units with Monte Carlo SE 0.045. The primary data support a small average bias and leave considerable uncertainty about the dataset-level false discovery rate.

The extended analysis included 86 alternative and 87 global-null datasets. MI detected 27 of 516 changes (5.23%), compared with 59 (11.43%) for EM and 62 (12.02%) for LCA. Its coverage was 99.22%. The median MI-to-EM SE ratio was 1.001, so the different critical values explain much of the interval widening. Mean MI bias across changed species was −0.029, but individual changed species had biases as large as −0.135 log units. A small pooled mean is insufficient evidence for unbiased effects at each taxon. Empirical FDR was zero for every method. With 0 of 87 global-null datasets yielding a BH discovery, the Wilson 95% upper bound on that probability was 4.2%.

A shared-checkout update limits the revision claim for the extended analysis. All 174 successful pooled executions, including one duplicate, reproduced the cabe340 pooling decisions from their saved component fits. Their largest absolute q-value difference was 1.12 × 10⁻¹⁶ or less. The more conservative primary subset and the source audit below make the distinction explicit.

<!-- generated-results -->
## Results for 86 alternative and 87 global-null datasets

### Primary comparison at the requested revision

This subset contains 14 alternative and 14 global-null datasets. It includes completed fits from before the shared checkout's first later commit and full CLI runs using the verified local export of cabe340. The extended tables below include 173 distinct datasets in all. The source-cohort field in results.tsv identifies every row.

| Method | Changed bias ± MCSE | Changed RMSE | Unchanged bias / RMSE | Changed coverage (%) | BH power (%) | Alternative FDR (%) | Global-null raw FPR (%) | Global-null BH FPR (%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LCA | +0.049 ± 0.043 | 0.323 | +0.021 / 0.365 | 98.81 | 4.76 | 0.00 | 1.43 | 0.00 |
| EM | +0.025 ± 0.045 | 0.330 | +0.011 / 0.368 | 97.62 | 2.38 | 0.00 | 1.25 | 0.00 |
| MI | +0.027 ± 0.045 | 0.330 | +0.012 / 0.368 | 98.81 | 0.00 | 0.00 | 1.07 | 0.00 |

| Method | Null datasets with any BH discovery | 95% Wilson CI | Global-null coverage (%) | Mean FMI, ambiguous (alt) | Mean FMI, unique (alt) |
| --- | --- | --- | --- | --- | --- |
| LCA | 0/14 | 0.0–21.5% | 98.39 | — | — |
| EM | 0/14 | 0.0–21.5% | 98.21 | — | — |
| MI | 0/14 | 0.0–21.5% | 99.29 | 0.0017 | 0.0008 |

### Extended comparison across all audited datasets

The remaining numerical tables use all accepted datasets, with each seed counted once. The revision audit records agreement of the pooled outputs with cabe340 applied to the saved component fits.

### Effect estimates under the alternative

Bias is estimated minus true log fold change after centering. Each ± value is one dataset-level Monte Carlo standard error. RMSE is in natural-log units. Up and down rows expose cancellation in pooled signed bias.

| Method | Species | Intervals | Bias ± MCSE | RMSE | Centered coverage (%) ± MCSE | Raw coverage (%) | Mean CI width |
| --- | --- | --- | --- | --- | --- | --- | --- |
| LCA | changed | 516 | -0.013 ± 0.019 | 0.373 | 98.64 ± 0.49 | 97.87 | 1.860 |
| LCA | up | 258 | -0.012 ± 0.025 | 0.355 | 98.84 ± 0.66 | 98.06 | 1.876 |
| LCA | down | 258 | -0.015 ± 0.028 | 0.390 | 98.45 ± 0.76 | 97.67 | 1.845 |
| LCA | unchanged | 2924 | -0.005 ± 0.005 | 0.373 | 98.19 ± 0.24 | 97.26 | 1.848 |
| LCA | ambiguous unchanged | 2064 | -0.015 ± 0.007 | 0.377 | 98.21 ± 0.26 | 97.14 | 1.845 |
| EM | changed | 516 | -0.030 ± 0.019 | 0.384 | 98.06 ± 0.58 | 97.87 | 1.864 |
| EM | up | 258 | -0.028 ± 0.026 | 0.356 | 98.45 ± 0.76 | 98.45 | 1.865 |
| EM | down | 258 | -0.032 ± 0.028 | 0.409 | 97.67 ± 0.92 | 97.29 | 1.863 |
| EM | unchanged | 2924 | -0.007 ± 0.005 | 0.376 | 98.12 ± 0.27 | 97.40 | 1.850 |
| EM | ambiguous unchanged | 2064 | -0.019 ± 0.007 | 0.381 | 98.06 ± 0.31 | 97.34 | 1.848 |
| MI | changed | 516 | -0.029 ± 0.019 | 0.385 | 99.22 ± 0.38 | 98.84 | 2.177 |
| MI | up | 258 | -0.027 ± 0.026 | 0.356 | 100.00 ± 0.00 | 99.22 | 2.176 |
| MI | down | 258 | -0.031 ± 0.028 | 0.411 | 98.45 ± 0.76 | 98.45 | 2.179 |
| MI | unchanged | 2924 | -0.006 ± 0.005 | 0.376 | 99.21 ± 0.16 | 98.77 | 2.159 |
| MI | ambiguous unchanged | 2064 | -0.018 ± 0.007 | 0.382 | 99.18 ± 0.18 | 98.69 | 2.157 |

Each changed species has an unchanged sister. The next table gives the mean centered error per changed species across alternative datasets. This distinguishes species-specific bias from cancellation between species.

| TaxID | True LFC | Ambiguity probability | LCA bias ± MCSE | EM bias ± MCSE | MI bias ± MCSE | MI − LCA, paired ± MCSE |
| --- | --- | --- | --- | --- | --- | --- |
| 100 | +1.386 | 0.30 | +0.010 ± 0.039 | -0.102 ± 0.037 | -0.101 ± 0.037 | -0.111 ± 0.005 |
| 105 | -1.386 | 0.50 | -0.025 ± 0.042 | -0.021 ± 0.042 | -0.020 ± 0.042 | +0.005 ± 0.003 |
| 110 | +1.386 | 0.60 | -0.051 ± 0.040 | +0.010 ± 0.042 | +0.011 ± 0.042 | +0.062 ± 0.004 |
| 120 | +1.386 | 0.50 | +0.005 ± 0.035 | +0.009 ± 0.035 | +0.009 ± 0.035 | +0.004 ± 0.002 |
| 125 | -1.386 | 0.40 | -0.014 ± 0.042 | +0.062 ± 0.041 | +0.062 ± 0.041 | +0.076 ± 0.004 |
| 135 | -1.386 | 0.70 | -0.006 ± 0.043 | -0.137 ± 0.047 | -0.135 ± 0.048 | -0.130 ± 0.007 |

Among the 14 unchanged species per dataset that share reads with a changed sister, positive signed error below means an apparent shift in the direction of that sister's true effect.

| Method | Species-datasets | Sister-direction error ± MCSE | BH FPR (%) | Centered coverage (%) |
| --- | --- | --- | --- | --- |
| LCA | 1204 | -0.011 ± 0.011 | 0.00 | 97.92 |
| EM | 1204 | -0.013 ± 0.012 | 0.00 | 98.01 |
| MI | 1204 | -0.013 ± 0.012 | 0.00 | 98.92 |

### Discoveries and global-null calibration

| Method | Power (%) | Alternative FDR (%) | Global-null raw FPR (%) | Global-null BH FPR (%) | Global-null FDR, any BH discovery | Global-null centered coverage (%) |
| --- | --- | --- | --- | --- | --- | --- |
| LCA | 12.02 ± 2.25 | 0.00 ± 0.00 | 1.01 ± 0.20 | 0.00 ± 0.00 | 0.00 (0/87, 95% CI 0.0–4.2%) | 98.51 ± 0.23 |
| EM | 11.43 ± 2.18 | 0.00 ± 0.00 | 0.98 ± 0.19 | 0.00 ± 0.00 | 0.00 (0/87, 95% CI 0.0–4.2%) | 98.42 ± 0.24 |
| MI | 5.23 ± 1.39 | 0.00 ± 0.00 | 0.78 ± 0.18 | 0.00 ± 0.00 | 0.00 (0/87, 95% CI 0.0–4.2%) | 99.25 ± 0.16 |

Global-null FDR uses a Wilson binomial interval across datasets. Other ± values are one Monte Carlo standard error in percentage points. The BH FPR denominator is all true-null species. Its value can be small even when several datasets have a discovery.

| Method | Robust power (%) | Robust alternative FDR (%) | Robust global-null FDR (%) | Raw global-null CI coverage (%) |
| --- | --- | --- | --- | --- |
| LCA | 12.02 | 0.00 | 0.00 | 97.53 |
| EM | 11.43 | 0.00 | 0.00 | 97.79 |
| MI | 5.23 | 0.00 | 0.00 | 99.22 |

| Method | Changed species with raw p < 0.05 (%) | Median changed-species p | Median changed-species BH q | Median changed-species SE |
| --- | --- | --- | --- | --- |
| LCA | 79.46 | 0.0150 | 0.177 | 0.470 |
| EM | 80.43 | 0.0145 | 0.177 | 0.470 |
| MI | 77.52 | 0.0173 | 0.213 | 0.471 |

Using the same 1.96 cutoff on the pooled MI SE gives 97.87% centered coverage for changed species, compared with 99.22% using the specified MI t interval. This diagnostic isolates the contribution of the critical value to coverage.

### Assignment uncertainty and paired comparisons

| Scenario | Species | Species-datasets | Mean FMI | Median FMI | 90th pct FMI | Median MI df |
| --- | --- | --- | --- | --- | --- | --- |
| alt | all ambiguous | 2580 | 0.0024 | 0.0007 | 0.0041 | 8.46 |
| alt | quartet | 1720 | 0.0025 | 0.0007 | 0.0042 | 8.46 |
| alt | pair | 860 | 0.0023 | 0.0006 | 0.0038 | 8.46 |
| alt | unique only | 860 | 0.0014 | 0.0001 | 0.0008 | 8.46 |
| global_null | all ambiguous | 2610 | 0.0022 | 0.0008 | 0.0049 | 8.45 |
| global_null | quartet | 1740 | 0.0021 | 0.0008 | 0.0050 | 8.45 |
| global_null | pair | 870 | 0.0022 | 0.0007 | 0.0047 | 8.46 |
| global_null | unique only | 870 | 0.0012 | 0.0001 | 0.0031 | 8.46 |

FMI is the share of pooled variance attributable to between-imputation differences. A uniquely identifiable species can have nonzero FMI because ANCOM-BC2 estimates its coefficient using bias corrections learned across species.

| Paired comparison | Scenario | Metric | Difference, percentage points ± MCSE |
| --- | --- | --- | --- |
| MI minus LCA | alt | power | -6.78 ± 1.59 |
| MI minus LCA | alt | changed CI coverage | +0.58 ± 0.43 |
| MI minus LCA | alt | FDR | +0.00 ± 0.00 |
| MI minus LCA | global_null | global-null FDR | +0.00 ± 0.00 |
| MI minus EM | alt | power | -6.20 ± 1.54 |
| MI minus EM | alt | changed CI coverage | +1.16 ± 0.46 |
| MI minus EM | alt | FDR | +0.00 ± 0.00 |
| MI minus EM | global_null | global-null FDR | +0.00 ± 0.00 |

For changed species the median absolute difference between the centered EM and MI coefficients was 0.0036 log units. The median ratio of MI SE to EM SE was 1.0010.

### Runtime and data checks

The 173 completed datasets contained 415,200,000 reads. Of these, 156,702,068 (37.74%) had multiple candidate species. The accepted datasets represent 1211 R fits, excluding repeated and failed attempts. Across methods 20,760 of 20,760 species-method rows had usable intervals. All 6,920 MI rows with five estimable fits can be identified in results.tsv.

| Stage | Median min/dataset | Range (min) | Cumulative min |
| --- | --- | --- | --- |
| Generate per-read lines | 0.06 | 0.04–0.39 | 14.7 |
| LCA table | 0.32 | 0.22–0.65 | 58.8 |
| EM table + 5 draws | 2.80 | 2.14–5.17 | 520.2 |
| LCA DA | 0.93 | 0.56–2.43 | 177.8 |
| EM expected DA + 5 imputed DA fits | 4.92 | 3.07–24.66 | 965.1 |

Elapsed wall time from first generation to last dataset result was 239.6 min (18:17:58 to 22:17:37 UTC). Dataset jobs overlapped, with up to eight slots and one thread per fit. Stage times for the accepted executions include subprocess startup and concurrent machine load. Superseded and failed attempts are recorded separately in the dataset JSON logs.

The analysis checked every LCA cell against singleton truth, expected-table library sizes against 200,000 (rounding tolerance 0.00004), integer-draw totals, unique-species counts and between-imputation variance against the five saved fits. The unrounded bias-corrected coefficients and original p/q values are in results.tsv.

A fixed-size replay of these seeds uses `python simulate.py --alt 1 86 --null 1 87 --workers 8` in a fresh study directory.
<!-- end-generated-results -->
## Simulation design

Each independently seeded dataset contains five invented genera with eight species per genus and twelve independent samples, six control and six treated. Each species receives one dataset-specific baseline log abundance from N(0, 0.85²). Each sample adds independent N(0, 0.65²) biological variation per species. The group coefficient on this absolute log-abundance scale is +ln(4) for species indices 0, 10 and 20 and −ln(4) for indices 5, 25 and 35. All 34 other coefficients are zero. All 40 coefficients are zero in global-null datasets. TaxIDs are 100 plus the zero-based index. Exponentiating and closing abundances to proportions, then sampling 200,000 reads multinomially, gives a fixed library depth. Each method receives the same reads and metadata within a dataset.

The first four species in each genus form a quartet. Reads from its four members have probabilities 0.30, 0.40, 0.60 and 0.70, respectively, of also hitting a uniformly chosen subset of sisters in that quartet. The number of additional hits is uniform on 1, 2 and 3. The next two species can hit each other with probability 0.50. The final two have unique reads only. Ambiguity probabilities are constant across groups. Every candidate list contains the true origin, with equal scores for all hits. All six changed species share reads with an unchanged sister. Four belong to quartets with overlapping candidate sets and two to disjoint pairs. Generation uses conditional binomial/multinomial draws and writes one physical hit line per candidate, keeping each read's lines consecutive. Grouping read IDs by candidate set saves generation time and preserves the information seen by the classifier-output parser.

The unequal ambiguity probabilities within quartets are a deliberate challenge to the allocation model. Given a shared hit set, its members can have different probabilities of producing that set. Allocation by abundance alone omits that difference. The equal-probability pairs provide a comparison under the same library depth and biological variation. Constant species-specific loss in LCA can affect abundance levels while canceling from within-species group ratios in expectation.

## Estimands and denominators

For each method and dataset, the simulation records the unmodified coefficient, standard error, p and q for treated versus control. For bias and RMSE, it subtracts the median estimated coefficient of the true-null species before comparing with the generating absolute log fold change. This uses knowledge of the simulation truth. Both centered and original-coefficient coverage are reported. The nominal 95% interval uses 1.96 × SE for each single fit and the 0.975 Student quantile at the reported Barnard–Rubin degrees of freedom times the pooled SE for MI. Centering shifts the interval intact and makes no allowance for the estimated offset's uncertainty. Coverage describes estimable intervals. Missing or untested species count as nondiscoveries for power and false-positive rates.

Raw type I error is the fraction of true-null species with p < 0.05. BH false-positive rate is the fraction of true-null species with q < 0.05. Power is the fraction of changed species called at q < 0.05. Empirical FDR is the mean across datasets of false discoveries divided by all discoveries, assigning zero to a dataset with none. Under the global null this equals the fraction of datasets with any discoveries. Tests use the pipeline's original p and q. The sensitivity-robust call is also recorded. Monte Carlo standard errors treat independently generated datasets as the units, preserving dependence between species in a library and between methods fitted to the same reads.

## Reproduction

For a fresh replay, copy the scripts, env.sh, REPORT.md and pinned/ into an empty directory. The frozen source sets the exact target version even if the shared checkout moves. To re-analyze this completed study in place, start at revision_audit.py after setting the environment. The runner refuses to reuse an existing dataset name during normal operation. The full command and exit status of each production CLI invocation are in the corresponding run.json, with the dataset directory as its working directory. The R wrapper sets its own seed to 11. The generator and taxa-table use seed 202611190 + 2 × replicate number + 1 for alternative datasets and 202611190 + 2 × replicate number for global-null datasets.

```bash
source env.sh
export PYTHONPATH="$PWD/pinned/src"
python simulate.py --alt 1 86 --null 1 87 --workers 8
python revision_audit.py
python analyze.py
python report_tables.py --update
```

Within each dataset the runner executes the following commands. The shell variable SEED denotes that dataset's recorded integer seed.

```bash
python -m unclassified_meta.cli taxa-table --reports reports.tsv --per-read \
  --nodes ../../nodes.dmp --names ../../names.dmp --taxonomy-name sim --rank S --assign lca -o lca
python -m unclassified_meta.cli taxa-table --reports reports.tsv --per-read \
  --nodes ../../nodes.dmp --names ../../names.dmp --taxonomy-name sim --rank S \
  --assign em --imputations 5 --seed "$SEED" -o em
python -m unclassified_meta.cli da --counts lca/counts.tsv --features lca/features.tsv --metadata meta.tsv \
  --group group --formula group --reference-level control --p-adjust BH --alpha 0.05 --threads 1 -o da_lca
python -m unclassified_meta.cli da --counts em/counts.tsv --features em/features.tsv --metadata meta.tsv \
  --group group --formula group --reference-level control --p-adjust BH --alpha 0.05 --threads 1 -o da_em --imputations
```

The final command publishes the ordinary expected-table fit as da_em/ancombc2_results.tsv and the Rubin-pooled fit as da_em/mi_results.tsv. Each dataset also retains the five component fits. The shared checkout began on branch em-imputation at cabe340713d94f9734eeed858b586c63535701de. The tracked source and README matched HEAD at the initial verification. Source provenance for the later runs is described below. The environment used Python 3.11.18, NumPy 1.26.4, pandas 2.2.3, R 4.5.3 and ANCOMBC 2.14.0. env.sh sets the requested Python, module path, Rscript and Matplotlib path and confines temporary files and caches to sim/.

## Execution notes

The first eight invocations of taxa-table exited with status 1 before fitting, with the exact message `ERROR: taxa-table: output folder lca lies inside the input .`. The input-directory argument included the intended output folder. Retrying with a reports.tsv sheet used identical input bytes and seeds. Each affected dataset retains initial_table_lca.log and initial_run.json. The pilot was alt_001, with alt_002–005 and null_001–003 generated concurrently. The continuation fills freed dataset slots, first targeting at least 12 alternative and 8 null datasets and then balancing the two scenarios. The planned stopping rule used measured run durations and a four-hour wall budget, with time reserved for verification. The shared-checkout failures stopped new submissions from that driver. Its active jobs drained before the final source audit. No stopping decision uses estimates, p values or coverage.

The pipeline reports that 40 taxa give limited information for estimating its sample-specific bias, and that two groups deactivate multi-group comparisons. These warnings are retained verbatim in the DA output and logs. Prevalence cutoff was 0.10, library cutoff zero, structural-zero testing off and pseudo-count sensitivity analysis on. The environment probe found that SciPy is absent from the requested Python environment. Interval calculation uses the specified R installation's qt function.

During the first interim report render, pandas interpreted the scenario label `null` as missing and the renderer raised `KeyError: 'null'`. The exported scenario is now `global_null`; aggregation has used the simulation truth throughout. An earlier syntax check on the renderer caught `f-string expression part cannot include a backslash` before it ran. Both corrections concern report generation.

## Scope and interpretation

Library depth and biological replication play different roles here. Each sample has 200,000 reads, but each coefficient compares six independent samples with six others. Most shared reads have substantial singleton evidence nearby. Even a high percentage of ambiguous reads can therefore contribute little between-imputation variance to the group coefficient. The estimates still depend on how each species generates ambiguous hits. In the quartets, the pipeline weights an eligible species by fitted abundance while the generator also gives candidate species unequal probabilities of producing a shared hit set. Imputation propagates uncertainty around that assignment model and inherits its systematic errors.

The CI and testing conventions also matter. The first component fit of alt_001 gives complete-data t degrees of freedom of 10 for all 40 species, recovered from each reported statistic and p value. With a small positive FMI, the implementation's Barnard–Rubin observed-data term approaches 10 × 11/13 = 8.46. At exactly zero FMI the implementation returns the complete-data df. Consequently, small changes to the pooled SE can accompany appreciably wider intervals and larger p values. The single-fit 1.96 intervals requested for this study also use a smaller critical value than a t interval with 10 df. Power uses the pipeline's delivered tests in every comparison.

The results apply to high-depth libraries, 40 synthetic species, stable within-group mapping behavior, independent samples, one continuous biological-variation model and five imputations. All ambiguous hit sets include the true origin. The simulation does not exercise erroneous candidate lists, contaminants, nearly unidentified taxa, pervasive zeros or group-dependent classification quality. The symmetric pairs and heterogeneous quartets test different parts of the allocator. Their comparison also changes the number and overlap of candidates, so it cannot isolate the effect of any one assumption. The per-species Monte Carlo errors describe repeated simulated datasets under these choices.

## Revision drift and verification

Another process advanced the shared branch during the study. The first later commit, aea3ae7, changed documentation at 18:57:02 UTC. Commit 537f358 subsequently changed the MI treatment of p = 1 with a nonzero test statistic. The working tree also received edits to structural-zero reporting while two live R processes were reading the wrapper. The shared-checkout scheduler stopped admitting new jobs when it saw a failure. All existing outputs and errors were retained. I exported `src` from the immutable Git object into `pinned/` inside sim/ and verified every blob against `git ls-tree -r cabe340 src`. This copy supplied the pipeline for retries and the final fresh datasets. The original checkout was only read by this study.

The primary estimates use the conservative cutoff above and the frozen runs. Datasets completed after the initial cutoff enter the extended comparison, even where only documentation had changed. Two separate failures occurred on the shared path. For alt_084 the pipeline reported `ERROR: ANCOM-BC2 failed on imputation 3 (exit 1); see the messages above`. For alt_085 it reported `ERROR: ANCOM-BC2 failed (exit 1); see the messages above` on LCA. Each R log ends with `Error: unexpected ')' in "r(0))"` and `Execution halted`. The frozen retries used the same seed, reads and model settings and reran all four CLI commands. Original attempts remain in their dataset directories and retries are in pinned_fit/. Two additional frozen reruns (alt_040 and null_040) checked later-source results directly. Null_040 completed. Alt_040 had byte-identical count tables, matching LCA and EM expected-table fits, and four matching imputation fits. Its optional duplicate was interrupted during the fifth imputation (driver exit 130) at the four-hour limit. Its original, complete result remains in the extended cohort only.

The recorded source diff is in observed_source.diff. The changes visible at audit time covered MI's handling of untestable features, bookkeeping and structural-zero outputs. Structural-zero testing was disabled for all runs. The R fitting arguments and the count allocation implementation were identical in the recorded diff. For every accepted dataset, revision_audit.py re-runs the cabe340 pooling code on the five saved fit tables and checks coefficients, SE, degrees of freedom, FMI, p, q, fit counts and decisions against the published output. Integer and logical decisions must match exactly. The comparisons use absolute tolerance 1e−10 and relative tolerance 2e−8 for floating-point columns, allowing for the round trip through TSV. Per-dataset maximum differences and the number of p = 1, nonzero-W cases are in repool_audit.tsv. These tests support the extended numerical comparison. Endpoint diffs alone cannot establish the exact source bytes read by every earlier subprocess.


The raw hit listings are retained in runs/ alongside their singleton and origin truth tables. An optional archive attempt was stopped (exit 130) when its throughput would have extended the study past the time budget. Its partial output was discarded. To pack these files later, `bash pack_reads.sh` creates a Zstandard-compressed tar archive and compares every member with its source before removing the redundant listings. `tar -I zstd -xf per_read_inputs.tar.zst` restores the original paths.


At 22:05 UTC the final optional launch batch was interrupted together (driver exit 130) as fit times rose. The original per-dataset logs remain. Alt_087 and null_088 had finished both count tables and LCA DA; alt_088 had finished its LCA count table. These three datasets have status `budget_stopped` in runtimes.tsv. Their partial outputs contribute to no metric. The three were distinct from the same-seed retries. The intended complete-case analysis gives every accepted seed equal weight.

For the frozen reruns the command was `python pinned_rerun.py --workers 4 alt_084 alt_085 alt_040 null_040`. Every subprocess used the environment in env.sh with PYTHONPATH redirected to the verified export. Generation and the original shared-path fits used the user-specified PYTHONPATH. Reruns retained the same Python executable, Rscript, seeds, counts model and thresholds. The successful frozen executions have exact invocations and exit codes in pinned_fit/run.json. The stopped duplicate has a cutoff record there. The original launch history is in pilot.log, wave1.log, their retry logs and continuation.log.

## Files

| File | Contents |
| --- | --- |
| [results.tsv](results.tsv) | One row per accepted dataset, method and species, including truth, ambiguity probability, source cohort, original and centered LFC, SE, p/q, interval and call |
| [strict_metrics.tsv](strict_metrics.tsv) | Primary revision-specific estimates and dataset-level Monte Carlo errors |
| [metrics.tsv](metrics.tsv) and [replicate_metrics.tsv](replicate_metrics.tsv) | Extended summary and its dataset-level denominators |
| [runtimes.tsv](runtimes.tsv) | Effective status, original status, stage times and source cohort for each seed |
| [repool_audit.tsv](repool_audit.tsv) | Original-version re-pooling compared with saved results for each successful execution, including duplicate checks |
| [simulate.py](simulate.py), [pinned_rerun.py](pinned_rerun.py), [revision_audit.py](revision_audit.py), [analyze.py](analyze.py), [report_tables.py](report_tables.py) | Reproducible generation, pipeline execution, version audit, statistical analysis and table rendering |
