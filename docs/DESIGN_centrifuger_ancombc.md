# Design: differential abundance of viruses / selected taxa (panel workflow)

Status: implemented (2026-11-18). Commands: `panel-build`, `panel-quant`, `da`, `taxa-table`, `control-reads`,
`control-summary`, `da --rand-formula`, `--min-score` / `--min-hitlen` for `run` and `extract`. Ideas for the
last five were contributed by Eli Holschbach.

## Problem

Centrifuger leaves a large unclassified fraction and its index cannot be enlarged (the prebuilt index already
needs about 217 GiB, close to the memory of a typical large compute node). `unclassified-meta run` explains the unclassified reads. The goal here is
a differential-abundance (DA) test of the viruses, or other selected taxa, found that way.

## What changed after the review

The first draft merged `centrifuger-quant` rows (classified reads) with read-back counts (unclassified reads)
into one table. A Codex review (Centrifuger v1.2.0 and ANCOMBC sources checked) found two blockers:

1. `centrifuger-quant` counts are cumulative over the taxonomy and `-k 1` can assign a read to an ancestor, so
   no set of report rows counts every read exactly once.
2. "Closest reference" labels and a pooled "novel" class are not validated biological features.

For a focused question (viruses, selected taxa) the review preferred one small joint reference and one counting
rule for all reads. That is what is implemented.

The classifier's reports came back in a second step, in a form that avoids blocker 1: `taxa-table` takes the
clade reads of the rows at ONE rank, never rows of two ranks, and reports the reads assigned above that rank in
a ledger instead of adding them (section "Classifier taxa as features").

## Workflow

```
unclassified-meta run          pooled assembly of the unclassified reads, contig classes   (existing)
unclassified-meta panel-build  features = selected contigs + the reference genomes they hit (+ extra targets)
                               decoys   = host genome(s), rRNA reference, the run's other contigs
unclassified-meta panel-quant  ALL reads of every sample (original FASTQ, full depth) -> Bowtie 2 -> one
                               outcome per fragment -> counts.tsv (feature x sample) + ledger.tsv
unclassified-meta da           ANCOM-BC2 on counts.tsv + sample metadata
```

Subsampling is decoupled from quantification: `run` may assemble from a subsample, `panel-quant` always counts
every read against the frozen panel.

## Classifier taxa as features (`taxa.py`)

Input: one Kraken-style report per sample (`kraken2 --report`, `centrifuger-kreport`,
`centrifuger-quant --output-format 3`). Column 2 (clade reads) of the rows whose rank code equals `--rank` is
the count; column 1 is never read (in quant reports it is the EM abundance). Rows at one rank are disjoint
clades, so no read is counted twice. `-k 1` assignments to an ancestor, and reads the classifier placed above
the rank for any other reason, are `unresolved_above_rank`: reported, not in the table. The per-sample ledger
`input = unclassified + dropped_lineages + kept_at_rank + unresolved_above_rank + rounding_residual` is built
from two columns independently (clade reads for kept and dropped, own reads for unresolved), so it is a real
check of the tree; a residual beyond the rounding of `centrifuger-quant` (0.5 read per row) is an error.
Feature ids are `namespace|taxID|rank`. Duplicate samples and a file suffix that matches two reports of one
sample are errors.

Second input path, `taxa-table --per-read`: Centrifuge(r) per-read output plus the index taxonomy (`nodes.dmp`,
`names.dmp`; for a prebuilt index `centrifuger-inspect --taxonomy-tree` / `--name-table`). No report is needed.
`extract.centrifuge_reads()` is the single function that decides, per read, between unclassified, below the
thresholds and classified (with the accepted taxIDs); `extract` and `taxa-table --per-read` both call it, so
`unclassified + below_threshold` of the taxa ledger equals the reads `extract` writes (unit-tested over several
thresholds). A classified read with one accepted taxID is projected to its ancestor at the rank; with several
(`-k` > 1) it goes to their lowest common ancestor first. This differs from `centrifuger-quant`, which splits
the read evenly over the taxa: LCA never credits a taxon the read only partly supports, at the price of more
reads in `unresolved_above_rank`. Ledger: `input = unclassified + below_threshold + dropped_lineages +
kept_at_rank + unresolved_above_rank` (no rounding). TaxIDs absent from the taxonomy are counted as unresolved
and reported (`not_in_taxonomy`). `--taxonomy-name` is mandatory and becomes the namespace of the feature IDs.

Two ways to put these taxa and a panel in one DA table (`taxa-table --merge-panel`); `panel-quant` records which
one applies as `reads_scope` in its provenance, and the merge follows the record:

| | scope `unclassified` (disjoint reads) | scope `all` |
|---|---|---|
| panel counted on | `unclassified_reads/` of `extract` / `run` (no-hit reads + low-score reads) | the original FASTQ |
| overlap of the two tables | none: a read is either in the report's classified part or in the extracted reads | reads of the panel's organisms that the classifier also assigned |
| rule | no lineage has to be dropped | `--panel-lineage` required; that lineage is dropped from the reports (`none` only as an explicit statement) |
| check | per sample: classified kept + dropped + unresolved (+ rounding) + panel features + panel rest = total input reads; the merge fails otherwise | not possible |

The reconciliation uses the total of `tables/input_summary.tsv` of the extract run when `panel-quant` recorded it,
otherwise the report's own unclassified row; without either the row is marked `not_checked`. It fails when the
report and the extracted reads were made with different thresholds, which is the loss described in finding 5.

Thresholds: `run` / `extract` accept `--min-score` and `--min-hitlen` for Centrifuge(r) per-read output with the
rule of `centrifuger-quant` v1.2.0 (`LoadReadAssignments`: a line is skipped when hitLength < min, score < min or
taxID = 0); a read is classified when any of its lines passes. Reads below the thresholds are extracted with the
unclassified reads and counted as `low_score_reads`.

## Counting policy (one rule, `fragments.py`)

Unit: fragment = read pair, or single read. Bowtie 2 runs with `-k 10`; for each mate the best-scoring
alignments (all alignments tied at that mate's best score) are kept and looked up in the sequence -> feature map.

| Outcome | Rule |
|---|---|
| a feature | all kept alignments of both mates fall in that one feature (several sequences of the same feature are fine; R1-only and R2-only fragments count) |
| `ambiguous` | kept alignments fall in two or more features |
| `ambiguous_decoy` | kept alignments fall in a feature and in a decoy |
| `decoy_host:<name>`, `decoy_rrna`, `decoy_other_contigs` | kept alignments fall only in that decoy class |
| `decoy_mixed` | only decoys, of more than one class |
| `unmapped` | no alignment |

The outcomes are mutually exclusive, their sum per sample equals the fragments Bowtie 2 read, and the run fails
when it does not. Only feature counts enter the DA table; the other outcomes stay in `ledger.tsv` as QC.

Deviation from the first wording of the policy ("every reported alignment"): lower-scoring alignments are
ignored. Without this, a read that matches its own contig exactly and the divergent reference genome poorly
would be called ambiguous, and most reads of conserved regions would be lost.

## Feature definition (`panel.py`)

A feature is a reference-anchored bin; its id is `namespace|accession|class`, with the reference release in the
namespace (for example `refseq_viral_r238|NC_040842.1|like`). Names are labels only.

| Class | Contigs in the bin | Label |
|---|---|---|
| `same` | BLASTN to one reference genome, read-weighted identity >= 95 % (`--same-identity`); the reference genome is added to the feature | X |
| `like` | BLASTN (or TBLASTX) to one reference genome below that | X-like |
| `reference` | the reference genome of a `like` bin, kept as its own feature; also each `--extra-targets` sequence | X (reference genome) |
| `protein_like` | DIAMOND blastx only; binned by the organism of the protein, so several proteins of one virus give one feature | X-like (protein only) |
| `kraken2` | Kraken2 taxid of the contig | name (Kraken2 contig assignment) |
| `unmatched` | no reference match; one feature per sequence cluster (single linkage, BLASTN >= 90 % over >= 50 % of the shorter contig) | unmatched sequence cluster N |

`features.tsv` also reports `identity_min`, `identity_max` and `mixed_identity` (10-90 % of the read weight on
each side of the threshold), which flags a bin that may hold two entities.

## ANCOM-BC2 settings (`r/ancombc2.R`, ANCOMBC version recorded in the output)

- `tax_level = NULL` (features are already defined), raw integer counts, no rounding step.
- `struc_zero = FALSE` unless `--struc-zero`: zero fragments is not proof of biological absence.
- `pseudo_sens = TRUE`; the call to use is `diff_robust_*`.
- `prv_cut` (default 0.10) and `lib_cut` (default 0) are options.
- Error below 2 samples per group. Warnings, also written into the result header: fewer than 5 samples in a
  group; fewer than 50 features after the prevalence filter.
- `--rand-formula "(1|unit)"` is passed to `rand_formula` (lme4 through ANCOM-BC2). The grouping column is
  checked first: error when it is missing or has one sample per unit; warnings for single-sample units and for
  fewer than 3 units in a group; the unit structure is written into the result header. lme4's default checks stay
  on; `--lme-relax-checks` switches two of them off (`check.nobs.vs.nlev` and `check.nobs.vs.nRE` ignored)
  with a warning.

## Review findings: status

| # | Finding | Status |
|---|---|---|
| 1 | cumulative `centrifuger-quant` rows | addressed: the panel does not use the report; `taxa-table` uses the rows of one rank only (disjoint clades), or with `--per-read` projects each read once from the per-read output, and reports reads above that rank as `unresolved_above_rank`. Open: a read above the rank cannot be projected down, so a species table omits it (a genus table omits fewer) |
| 2 | nearest-reference labels, "novel" pool | addressed: identity class in the feature id, reference genomes separate from `like` bins, protein hits binned by organism, unmatched contigs as sequence clusters, class totals and unmapped reads only in the ledger. Open: bins are not validated species (see limits) |
| 3 | subsampling coupled to quantification | addressed: `panel-quant` counts all reads against the frozen panel |
| 4 | `centrifuger-quant` format assumptions | addressed in `taxa-table`: Kraken-style layout only (formats 0-2 are refused as non-reports), counts from column 2, never the abundance column; no extra rounding (the report's own row rounding is carried as `rounding_residual`). The format-3 report has no unclassified row; the total then comes from the extract summary |
| 5 | reads in = reads counted | addressed: panel ledger with exclusive outcomes, asserted; `--min-score` / `--min-hitlen` in `extract` so that quant-rejected reads land on the unclassified side; `classified + low_score + no_hit = total` asserted; reconciliation of both sides in the merge (scope `unclassified`) |
| 6 | first-mate counting, no ambiguity policy | addressed: one fragment policy; the read-back of `run` now also counts R2-only fragments |
| 7 | taxonomy joins by name | addressed inside the panel and in `taxa-table` (namespaced ids, taxID keys, names as labels; `--drop-lineage` also takes taxIDs). Open: protein-only bins group by organism name within one RefSeq release; there is no crosswalk between panel features and Centrifuger/GTDB taxa, so in scope `unclassified` one organism can be two rows (classifier taxon and X-like bin) |
| 8 | ANCOM-BC2 requirements | addressed as checks and warnings; random effects with checks of the unit column; small panels can be given background features with `taxa-table --merge-panel` |
| 9 | host-like classes are not proof of host | partly addressed: decoys remove host and rRNA reads from the features and report them in the ledger; a decoy hit is "aligns to the host/rRNA reference", not proof of origin |
| 10 | RNA vs DNA interpretation | documented in the README (RNA: taxon-associated RNA signal; DNA: genome abundance) |
| 11 | stale counts, lossy tables | addressed: step markers and count caches are keyed on parameters and inputs; `run` writes raw `contig_counts.tsv`; provenance JSON with input hashes and tool versions |
| 12 | integer / container requirements | addressed: plain matrix input, `tax_level = NULL` |

## Limits and open issues

- **Small panels.** ANCOM-BC2 estimates each sample's bias (sampling fraction) from the features and assumes most
  of them do not differ. A virus panel often has fewer than 50 features; then the log fold changes are relative
  to the other panel features and a shift shared by all of them is not detectable. The script warns. Remedy,
  implemented: `taxa-table --merge-panel` appends the classifier's taxa at one rank as background features
  (synthetic test: 8 panel + 60 background features, no warning, the one differential virus found, no false
  positive, in both merge modes). Residual assumptions: both tables come from the same read files; in scope
  `all`, organisms outside the dropped lineage that are in the panel and in the classifier database would be
  counted twice; features of the two tables are counted by different methods, which is harmless only while the
  efficiency of a feature does not vary between samples.
- **Reads above the rank.** `taxa-table` does not project ancestor-level assignments down to the chosen rank;
  they are reported in the ledger only. Their share depends on the database and on `-k`.
- **Chance hits.** Classifier assignments at permissive settings can be mostly chance hits; `control-reads` and
  `control-summary` measure this with reversed reads, but the classifier run in between is the user's.
- **Mixed models.** With few units per group the unit variance is poorly estimated and the test has little
  power (synthetic check: 4 units x 3 replicates per group found none of three 8-fold features after Holm
  correction over 60 features; 6 units x 2 replicates found both of two over 24 features). ANCOM-BC2 with
  `rand_formula` and `pseudo_sens` took 4 min for 24 features x 24 samples on 4 cores.
- **Bins are closest-reference bins.** Two viruses with the same closest reference and the same identity class
  share a feature; one virus whose contigs hit two references is split. `mixed_identity` flags only the first
  case in part. Check the contigs of a significant feature before naming a species.
- **`-k 10`.** A read with more than 10 equally good alignments can hide a hit in another feature.
- **Recovery differs between samples.** Reads of a strain that is far from the pooled contigs align less well;
  a sample-specific strain therefore lowers its own counts. Not corrected.
- **Host decoy for RNA-seq.** Bowtie 2 does not align spliced reads; `decoy_host` is a lower bound of host reads.
- **Step cache keys** use the absolute path, size, modification time (ns) and a digest of the first and last
  megabyte of each input (whole file below 2 MB, or with UNCLASSIFIED_META_FULL_HASH=1), plus a cache version.
  A change in the middle of a large FASTQ that keeps its size and time is not detected.

## Code review 2 (2026-11-18): fixes

| # | Finding | Fix |
|---|---|---|
| 1 | a read ID reopened in a per-read file could be extracted and counted | the shared iterator refuses any reopened ID (Kraken2 too) |
| 2 | reconciliation accepted 101 vs 100 | exact equality; a report's rounding is its explicit `rounding_residual`, allowed only for estimated-abundance reports (centrifuger-quant format 3); read-share reports must add up exactly |
| 3 | an interrupted panel rebuild could validate stale counts | build in `.staging/`, move into place, write `provenance.json` (with sha256 of every artifact) last; `panel-quant` verifies the artifacts |
| 4 | markers without keys were trusted | they are recomputed; every key includes a cache version |
| 5 | classification and fallback BLAST/DIAMOND databases keyed on paths | keyed on reference identity (see file_id); an adjacent BLAST database older than its FASTA is not used |
| 6 | inputs identified by name, size and second | absolute path, size, mtime in ns, sampled content digest |
| 7 | TBLASTX identity mixed into nucleotide identity | identity of a nucleotide bin comes from BLASTN hits only |
| 8 | structural zeros vanished | `structural_zeros.tsv` and a header line |
| 9 | collaborator results in the README | removed; acknowledgement by name only |
| 10 | a panel lineage below the rank could double count | refused, for reports and for per-read input |
| 11 | missing references silently skipped | configured but missing references are errors; leave one out with --no-host-decoy / --no-rrna-decoy / --no-reference-genomes |
| 12 | extra-virus reference genomes never in the panel | pulled from the run's extra virus FASTA under its own namespace |
| 13 | design checked before ANCOM-BC2's sample filters | checked and reported again on the retained samples |
| 14 | numeric-looking sample IDs | read as text in Python and R; matched by name |
| 15 | a sample named like an output column | reserved names refused |

## Code review 3 (2026-11-18): fixes

| # | Finding | Fix |
|---|---|---|
| 1 | classifier none symlinked inputs into the output; a later run wrote through the link | real copies; every output written to a temporary file and renamed; a link at the output path is removed first; an output that is an input file is refused |
| 2 | missing or duplicate selected reads, unsynchronised mates | every selected read found exactly once, FASTQ record count = classifier total, mates compared record by record by normalised ID |
| 3 | generated names could collide | output read file names and reference-table columns checked across the sample set |
| 4 | detection missed STAR / `_clean` names and chose silently | STAR `.Unmapped.out.mate1/2` and `<s>_clean_R1/R2` recognised; classifier outputs matched by exact sample name; orphan or ambiguous outputs are errors; matching table printed |
| 5 | markers outlived failed reruns; quant tables published before provenance | markers removed before a step changes outputs and written last; panel-quant stages its tables and writes a checksummed provenance last; the merge verifies it |
| 6 | databases not bound to their sources | sidecar with the source sha256 (whole file) for every BLAST/DIAMOND database; builds go to a temporary name; a database without a matching sidecar is rebuilt |
| 7 | a panel lineage also in --drop-lineage escaped the rank check | always checked |
| 8 | missing reference genomes only warned | error, unless --no-reference-genomes |
| 9 | formula variables and missing values | variables must exist, the group must be in the formula, samples with missing values are left out and named, then replication is checked |
| 10 | --struc-zero turned on neg_lb | separate --neg-lb (off); "absent" means zero reads in every retained sample of the group |
| 11 | report format guessed | --report-format required; clade = own + children checked per row (exact, or within the rounding of centrifuger-quant) |
| 12 | versions of reused downloads and panel namespaces | reused databases keep their recorded version; the run records the references it used (references.used.json) and panel-build takes its namespaces from there and refuses a different --refs |
| 13 | identities of different methods averaged together | blastn_identity, tblastx_identity, diamond_identity, minimap2_identity |
| 14 | pooled Sankey weights used subsampled counts | each sample's share is scaled to its unclassified reads before pooling |
| 15 | hygiene | cluster-specific wording removed; tool caches in .gitignore |

Memory of the read-ID checks: about 130-180 bytes per read (2.7 GB for 15 million reads in `extract`).

## Code review 4 (2026-11-18): general mechanisms

- **Output safety** (`io_utils.refuse_overlap`, `io_utils.staged`): every command first collects all its outputs
  (and their temporary names) and all inputs of all samples, mates and references, and refuses when any output is
  an input (same path, symbolic link, or hard link: same device and inode). Every file is written to a temporary
  name in its destination folder and renamed into place. Used by extract/run, control-reads, panel-build,
  panel-quant, taxa-table and download-refs.
- **Generation manifests** (`manifest.py`): a producing step removes its manifest before it changes anything and
  writes it last with the sha256 of every artifact; consumers verify it first. `run` (MANIFEST.json, read by
  panel-build), panel-build and panel-quant (provenance.json, read by panel-quant, taxa-table --merge-panel and
  da), taxa-table (MANIFEST.json, read by the merge and da), download-refs (MANIFEST.json per reference folder,
  read by run and by download-refs before reusing a download). `da` refuses a folder that looks like a
  panel-quant or taxa-table output (ledger.tsv) without a valid manifest.
- **Extraction**: every FASTQ read must be one of the classifier's reads, exactly once, and none may be missing;
  paired files are checked record by record in every mode, also without a classifier.
- **Detection**: classifier files are recognised by their content (Kraken2 C/U lines, Centrifuge(r) numeric
  columns, Kraken-style report rows), plain or gzipped; two candidates for one sample are an error.
- **Decoys**: a configured decoy must hold a sequence; only categories whose reference is loaded lose their
  contig decoys.
- **Provenance**: versions are read per reference key; a reused reference without a record is "unknown".
- **DA labels**: coefficients are labelled "X vs reference", "interaction ..." or "covariate ...".
- **Sankey**: per-sample shares are kept at full precision; rounding is for display only.

## Safety model

Guarantees (each enforced by one mechanism used by every command, not by checks at single sites):

| Guarantee | Mechanism |
|---|---|
| No input is ever modified | `cli._guard` → `io_utils.preflight` before the first write of every command (log file and sample sheet included): no input inside an owned output/scratch folder, no owned folder inside an input folder, no file in an owned folder that is an input by path, symbolic link or hard link. All files are written through `io_utils.staged` (unique `mkstemp` name, `os.replace`), which replaces a link instead of following it. Deletion is limited to staging folders made by `io_utils.make_staging_dir` (ownership marker; symbolic links refused). |
| Interrupted or concurrent runs never publish mixed outputs | `io_utils.folder_lock` (exclusive `flock` per output folder); `manifest.begin` writes the producer marker and removes the old manifest before anything changes, `manifest.seal` writes the checksums of a non-empty artifact list last; R writes into a staging folder that `da.run_ancombc2` publishes after the overlap check; `runner.Steps` reuses a cached step only when its recorded outputs still have their recorded sha256, and classification / read-back keys use the contigs' content. |
| Every consumer verifies its inputs | `manifest.verify` / `manifest.verify_file` (the exact consumed file must be listed); a folder with a producer marker but no complete seal is refused; panel-build compares every consumed reference (overrides included) with `references.used.json`; reference folders need a seal (`download-refs`, or `seal-refs` for a folder made by hand); downloaded parts carry a release record and are refetched for another release. Hand-made tables: `da --unsealed-input` only. |

Out of scope: hand edits inside a sealed folder while a command reads it; deliberate symbolic-link or file swaps
between the preflight check and the write (a time-of-check/time-of-use race by a malicious local user).

## Code review 5 (2026-11-18): safety model completed

1. Preflight before the first write of every command (`cli._guard`), including the CLI's own sample sheet and
   log; R output staged and published atomically by the Python wrapper.
2. Unique temporary names (`mkstemp`) for every staged file and manifest; per-folder lock.
3. Cached steps verify their outputs' sha256; classification and read-back keyed by the contigs' content.
4. taxa-table `--merge-panel` and `da` consume only files listed in a verified manifest.
5. Producer marker written first; a marked folder needs a complete, non-empty seal; `da --unsealed-input` for
   hand-made tables.
6. panel-build verifies each consumed reference against the run's recorded hashes.
7. Download cache bound to URL + release (sidecar with sha256); SILVA checked against its published MD5.
8. Control sampling checks seqtk's exit status inside the staging context.
9. `unclassified-meta seal-refs <dir>` adopts a reference folder prepared by hand; the run's error names it.
10. Sankey reference flows split by the largest-remainder method, so they sum exactly to their parent.

## Code review 6 (2026-11-18): scoped review fixes

1. Indirect inputs (reports of a report sheet, references a panel inherits from the run's config.used.yaml,
   sample-sheet FASTQs) are resolved before `cli._guard`, so they are part of the preflight.
2. `io_utils.preflight` refuses any symbolic link and any file with a second hard link inside an owned folder;
   logs are opened by `runner.open_log` (O_NOFOLLOW, single hard link).
3. `manifest.private_copies`: `da` and `taxa-table --merge-panel` read copies in their own staging folder,
   each checked against the manifest checksum; R reads only the copies.
4. `extract` (and `run`) seal `unclassified_reads/` together with `tables/input_summary.tsv`;
   `manifest.verify_consumed` makes panel-quant, control-reads and extract refuse reads from an unsealed or
   changed generation.
5. Every step marker records its outputs (the pool files included); a marker without outputs is not trusted.
6. The README control recipe writes its summary to `control_summary/`; tagged README blocks
   (`<!-- readme-test: ... -->`) are run by the test suite on the synthetic data.
7. SILVA downloads need the published MD5 (or `--allow-unverified-silva`, recorded in REFS_VERSION.tsv and the
   manifest); a cached part is reused only when it was verified against the same checksum.
Also: a `da` rerun removes outputs it did not write again; no site-specific node sizes in the docs.
