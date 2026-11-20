# Design: viral haplotype evidence, EM uncertainty and differential-abundance robustness

_2026-11-19. Extends the panel workflow (`panel-build` → `panel-quant` → `da`) and the EM work of
`src/unclassified_meta/em.py`. Nothing here changes a default of an existing command._

## Principles

1. **Sequence clusters are not haplotypes.** A 97/98/99 % identity cluster of contigs is an *operational
   sequence cluster*. It becomes a *candidate strain-like lineage* only when its members differ from other
   clusters inside homologous, overlapping regions (not merely by covering other parts of the genome), a
   *phased viral haplotype* only when reads link its distinguishing variants, and a *reference-supported viral
   strain* only when its consensus matches a reference genome closely over most of that genome. Every label is
   traceable to the rows that justify it.
2. **Observed counts are never overwritten.** Five matrices are kept apart: unique fragments, EM expected
   fragments, bootstrap realizations, normalized abundance, taxonomically aggregated abundance.
3. **Resampling is not replication.** Bootstrap replicates, EM/posterior realizations and technical replicates
   measure read-sampling and assignment instability. They are never counted as biological replicates, and their
   p-values are never pooled (Rubin's rules in `da --imputations` remain the one justified pooling, for the
   ANCOM-BC2 test of posterior imputations).
4. **Agreement is not validation.** ANCOM-BC2, ALDEx2 and LinDA are complementary frameworks applied to the same
   reads. Agreement is reported as *"robust across three complementary differential abundance frameworks"*,
   never as *"validated by three independent statistical methods"*.
5. **Dimensions stay separate.** Sequence-assignment confidence, abundance-estimation stability,
   differential-abundance support and cross-method agreement are reported as separate columns. The only
   combination is a rule-based `conclusion` (`high_confidence` / `exploratory`) whose rule is printed with it.

## Data flow

```
run ──► strain-evidence ──► viral_strain_clusters.tsv ──► panel-build --strain-clusters F --strain-resolution 98
            │                (+ overlap, haplotype tables)          (one feature per operational cluster)
            ▼
     reads (optional) → local phasing                       panel-quant (now also writes
                                                             equivalence_classes.tsv, unique_counts.tsv)
                                                                   │
                                                             em-uncertainty ──► bootstrap/ realizations/ matrices
                                                                   │            em_bootstrap_abundance.tsv
                                                                   │            em_assignment_uncertainty.tsv
                                                                   │            viral_unique_read_support.tsv
                                                                   ▼
                                                             da-robustness (ANCOM-BC2, ALDEx2, LinDA on the
                                                               primary, unique-only, bootstrap, realization,
                                                               leave-one-out and prevalence-filter matrices)
                                                                   │
                                                             da-resolutions (97 / 98 / 99 % panels)
```

## 1. `strain-evidence` (module `strains.py`)

Input: a finished `run` folder (`contigs/contigs.fa`, `tables/contig_category.tsv`, `references.used.json`),
the viral reference FASTA (the run's, or `--reference`), optionally a sample sheet / `--input-dir` for phasing.

1. Viral contigs with a BLASTN reference hit are grouped by that reference. Each group is aligned to its
   reference genome with `minimap2 -c --cs -x asm20` (all alignments, incl. supplementary).
2. **Overlap analysis.** For each pair of contigs on one reference: reference interval of each, overlap (bp),
   SNVs and indels between the two *inside the overlap only*, identity over the overlap, relation
   (`overlapping_identical`, `overlapping_divergent`, `disjoint` = a coverage difference, not divergence).
3. **Operational clusters** at each `--identity` (default 97,98,99): single linkage on pairs with overlap
   ≥ `--min-overlap` (200 bp) and overlap identity ≥ threshold. Contigs disjoint from all others stay singletons
   flagged `unresolved_disjoint` (they may belong to any cluster).
4. **Chimera / artifact flags**: split or non-collinear alignments, segments best matching different
   references, a jump in windowed identity along the contig, low aligned fraction.
5. **Variant phasing** (with reads): fragments aligned to the reference(s) with Bowtie 2 (pairs kept together);
   SNVs/indels called from the pileup (depth, alt count, allele frequency, base quality, both strands, binomial
   test against `--error-rate`); variant pairs linked by fragments covering both sites (2×2 tables, r², D');
   local blocks = consecutive variant sites spanned jointly by ≥ `--min-link` fragments; local haplotypes =
   complete allele patterns over a block, frequencies with bootstrap intervals; a pattern explainable as
   sequencing error of a more frequent neighbour is flagged `possible_error`. Read length, insert-size
   distribution, depth and the longest phasable span are reported per sample and reference.
   **Blocks are local**: no genome-wide haplotype is ever produced from disconnected blocks.
6. **Classification** per cluster and resolution (`viral_strain_clusters.tsv`):
   `operational_sequence_cluster` → `candidate_strain_lineage` (divergent from another cluster in an overlap
   ≥ min-overlap, not chimeric) → `phased_viral_haplotype` (≥ 2 of its distinguishing variants carried by one
   supported local haplotype; block ids listed) → `reference_supported_viral_strain` (consensus identity ≥
   `--reference-identity` 99 % over ≥ `--reference-cover` 80 % of the genome). The evidence columns name the
   overlap pairs and haplotype blocks used.

Outputs: `viral_cluster_overlap.tsv`, `viral_haplotype_evidence.tsv`, `viral_strain_clusters.tsv`,
`phasing_limits.tsv`, `figures/overlap_tracks_<ref>.png`, `figures/linkage_<ref>_<block>.png`, sealed
`MANIFEST.json`.

`panel-build --strain-clusters viral_strain_clusters.tsv --strain-resolution 98` splits each reference bin into
one feature per operational cluster (`feature_id` gets `|c98_NNN`, column `parent_feature` keeps the bin), so
that the panel can be quantified at each resolution.

## 2. `em-uncertainty` (module `uncertainty.py`)

Input: a sealed `panel-quant` folder (any `--ambiguous` mode), which now also writes
`equivalence_classes.tsv` (sample, candidates, fragments; every fragment once) and `unique_counts.tsv`.

- **Matrices**: `unique_counts.tsv` (copied), `em_expected_counts.tsv`, `bootstrap/counts_NNN.tsv`
  (`--bootstraps`, default 100), `realizations/counts_NNN.tsv` (`--realizations`, default 20; integer draws from
  the posterior assignment, `em.imputations`), `normalized_abundance.tsv` (EM fragments per million input
  fragments), `aggregated_counts.tsv` (by `--aggregate-by` column of features.tsv, `parent_feature` when present,
  or a `panel-taxonomy` table).
- **Read bootstrap**: within each sample, the fragments are resampled with replacement (multinomial over the
  equivalence classes, the unmapped class included, total = input fragments). A fragment is a read pair, so pairs
  stay together. The EM is refitted on every replicate. Seeds: `numpy.random.default_rng([seed, 3, sample_index])`.
- **Alternative mapping criteria**: unique only, even split of shared fragments, EM, and "any candidate"
  (every shared fragment counted for each of its candidates: an upper bound).
- `em_bootstrap_abundance.tsv` (feature × sample): EM point, bootstrap mean, SD, 2.5/97.5 percentiles, CV (when
  the mean ≥ 1 fragment).
- `em_assignment_uncertainty.tsv` (feature × sample): unique, shared (fragments in classes with this feature and
  others), EM, unique fraction of EM, realization SD and interval, even-split and any-candidate counts.
- `viral_unique_read_support.tsv` (feature): totals of the above, `em_driven` (unique fraction < 0.5),
  stability ratio across criteria, bootstrap CV, sequence evidence (identity class, mixed identity), and
  `assignment_confidence` (`high`, `moderate`, `low`) with the rule written in the file header.
- Figures: unique vs shared support per feature, bootstrap abundance distributions.

## 3. `da-robustness` (module `robustness.py`; R scripts `r/aldex2.R`, `r/linda.R`)

Methods: ANCOM-BC2 (existing `r/ancombc2.R`), ALDEx2 (`aldex.clr` + `aldex.glm`, Monte-Carlo Dirichlet
instances, `--mc-samples` 128), LinDA (`MicrobiomeStat::linda`, fixed or mixed formula). All three receive the
same samples, formula, reference level, and the same prevalence-filtered feature set (the filter is applied in
Python before any method runs, then the methods' own filters are switched off). Every method writes one
standardized table: `feature_id, term, effect, se, ci_low, ci_high, stat, p, q, effect_definition`.

- ANCOM-BC2: bias-corrected natural-log fold change; ALDEx2: CLR-scale glm coefficient (mean over Monte-Carlo
  instances); LinDA: bias-corrected CLR-scale coefficient (log2). Effects are compared by **direction and rank**,
  never as one quantity.
- ALDEx2's Dirichlet Monte Carlo models the sampling of counts *given* an integer table; it does not know which
  feature an ambiguous read belongs to. EM expected counts are rounded for ALDEx2 (documented); assignment
  uncertainty enters only through the realization and bootstrap tables.

Analyses (each method; `--bootstrap-fits`, `--realization-fits` choose how many tables are refitted):
primary (EM expected), unique-only, bootstrap tables, posterior realizations, leave-one-sample-out (or
leave-one-unit-out with `--unit`), alternative prevalence cuts (`--prv-cuts`).

Outputs:
- `differential_abundance_stability.tsv` (method × feature × term): primary effect/q, bootstrap effect median
  and 2.5/97.5 %, **directional stability** = share of resampled fits with the primary's sign,
  **significance stability** = share with q < alpha, the same for realizations, unique-only effect/q, filter
  sensitivity.
- `cross_method_concordance.tsv` (feature × term): methods tested, methods with q < alpha, sign per method,
  direction agreement, effects with their definitions, rank correlation summary in the header. No p-value pooling.
- `replicate_sensitivity.tsv`: per left-out sample/unit and summary rows: direction and significance stability,
  `single_sample_driven`, group sizes, largest single-sample share of a group's abundance, minimal detectable
  effect at 80 % power from the primary SE.
- `robustness_summary.tsv`: one row per feature × term with the four dimensions as separate columns
  (`assignment_confidence`, `abundance_stability`, `da_support`, `cross_method_agreement`), replicate
  consistency, `robustness_basis` (which of method agreement / effect-direction stability / read-assignment
  stability / clustering-resolution stability / biological replicate consistency hold) and `conclusion`.
- `robustness_report.md`: design, n per group and power limits, wording rules, high-confidence vs exploratory
  findings, software versions, parameters, seeds.
- Figures: bootstrap effect distributions, direction-stability heatmap, concordance plot, replicate-level
  abundance plot.

`da-resolutions --run 97=DIR --run 98=DIR --run 99=DIR` → `multi_resolution_robustness.tsv`: per parent feature,
significance and direction at each resolution, consistent or not.

## Reproducibility

Every new folder is sealed with `MANIFEST.json` (parameters, seeds, input checksums, Python package versions,
R `sessionInfo()` per method, external tool versions); consumers verify before reading.
