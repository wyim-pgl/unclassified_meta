#!/usr/bin/env Rscript
# Differential abundance with ANCOM-BC2 on a feature x sample count table (unclassified-meta panel-quant).
#
# Usage: Rscript ancombc2.R --counts counts.tsv --features features.tsv --metadata metadata.tsv \
#          --group group --outdir da [--formula "group + batch"] [--reference-level control] \
#          [--prv-cut 0.10] [--lib-cut 0] [--struc-zero] [--alpha 0.05] [--p-adjust holm] [--threads 1]
#          [--ledger ledger.tsv]   (panel-quant ledger: the abundance plot then shows shares of all fragments)
#          [--rand-formula "(1|unit)"] [--lme-relax-checks]   (repeated measures: random intercept per unit)
#
# Random-effects support (repeated samples of one unit) is an idea contributed by Eli Holschbach.
#
# Settings follow the design note (docs/DESIGN_centrifuger_ancombc.md): features are already defined, so
# tax_level = NULL; struc_zero is off unless asked for (no reads is not proof of biological absence);
# pseudo_sens = TRUE and the sensitivity-robust calls (diff_robust_*) are reported next to the raw ones.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(formula = NULL, `reference-level` = NULL, `prv-cut` = "0.10", `lib-cut` = "0", `struc-zero` = FALSE,
            alpha = "0.05", `p-adjust` = "holm", threads = "1", outdir = "da", `rand-formula` = NULL,
            `lme-relax-checks` = FALSE, `neg-lb` = FALSE)
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (key %in% c("struc-zero", "lme-relax-checks", "neg-lb")) { opt[[key]] <- TRUE; i <- i + 1; next }
  if (i == length(args)) stop("option --", key, " needs a value")
  opt[[key]] <- args[i + 1]; i <- i + 2
}
for (need in c("counts", "features", "metadata", "group")) {
  if (is.null(opt[[need]])) stop("missing required option --", need)
}
suppressPackageStartupMessages({ library(ANCOMBC); library(ggplot2) })
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
warn <- character(0)
one_line <- function(x) gsub("\\s+", " ", trimws(x))
note <- function(...) { m <- one_line(paste0(...)); warn <<- c(warn, m); message("WARNING: ", m) }

# identifiers are text: "01" stays "01", and "101" is a name, never a column position
counts <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c(feature_id = "character"))
feats <- read.delim(opt$features, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c(feature_id = "character"))
# na.strings = character(0): a sample called "NA" stays a name; empty cells or "NA" in a model variable are
# treated as missing below
meta <- read.delim(opt$metadata, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c(sample = "character"),
                   na.strings = character(0))
if (!"sample" %in% names(meta)) stop("metadata needs a 'sample' column")
if (names(counts)[1] != "feature_id") stop("the first column of the counts table must be 'feature_id'")
group <- opt$group
if (!group %in% names(meta)) stop("metadata has no column '", group, "'")
mat <- as.matrix(counts[, -1, drop = FALSE]); rownames(mat) <- counts[[1]]
if (any(!is.finite(mat)) || any(mat < 0)) stop("counts must be finite and non-negative")
if (any(duplicated(meta$sample))) stop("duplicate sample names in the metadata")
miss <- setdiff(meta$sample, colnames(mat))
if (length(miss)) stop("samples in the metadata without counts: ", paste(miss, collapse = ", "))
extra <- setdiff(colnames(mat), meta$sample)
if (length(extra)) note("samples with counts but not in the metadata were left out: ", paste(extra, collapse = ", "))
mat <- mat[, match(meta$sample, colnames(mat)), drop = FALSE]
rownames(meta) <- meta$sample

# model variables: every formula variable must be a metadata column, and the group must be in the formula
formula <- if (is.null(opt$formula)) group else opt$formula
fvars <- tryCatch(all.vars(as.formula(paste("~", formula))), error = function(e) stop("--formula '", formula,
                  "' is not a valid model formula: ", conditionMessage(e)))
miss_f <- setdiff(fvars, names(meta))
if (length(miss_f)) stop("--formula uses column(s) not in the metadata: ", paste(miss_f, collapse = ", "))
if (!group %in% fvars) stop("the group variable '", group, "' is not in --formula '", formula,
                            "'; the tested group must be a term of the model")
rvars <- character(0)
if (!is.null(opt$`rand-formula`)) {
  rvars <- tryCatch(all.vars(as.formula(paste("~", gsub("\\|", "+", opt$`rand-formula`)))), error = function(e) character(0))
  miss_r <- setdiff(rvars, names(meta))
  if (length(miss_r)) stop("--rand-formula uses column(s) not in the metadata: ", paste(miss_r, collapse = ", "))
}
model_vars <- unique(c(fvars, rvars))
is_missing <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "NA")
na_rows <- Reduce(`|`, lapply(model_vars, function(v) is_missing(meta[[v]])))
if (any(na_rows)) {
  note(sum(na_rows), " sample(s) with missing values in ", paste(model_vars, collapse = ", "), " were left out: ",
       paste(meta$sample[na_rows], collapse = ", "))
  meta <- meta[!na_rows, , drop = FALSE]
  mat <- mat[, meta$sample, drop = FALSE]
}
for (v in setdiff(model_vars, group)) {      # numeric covariates stay numeric, others become factors
  num <- suppressWarnings(as.numeric(meta[[v]]))
  meta[[v]] <- if (!anyNA(num)) num else factor(meta[[v]])
}
meta[[group]] <- factor(meta[[group]])
if (!is.null(opt$`reference-level`)) {
  if (!opt$`reference-level` %in% levels(meta[[group]])) stop("reference level not in column ", group)
  meta[[group]] <- relevel(meta[[group]], ref = opt$`reference-level`)
}
n_group <- table(meta[[group]])
if (nlevels(meta[[group]]) < 2) stop("column ", group, " needs at least two groups")
if (any(n_group < 2)) stop("every group needs at least 2 samples; found ",
                           paste(names(n_group), n_group, sep = "=", collapse = ", "))
if (any(n_group < 5)) note("fewer than 5 samples in a group (", paste(names(n_group), n_group, sep = "=", collapse = ", "),
                           "): ANCOM-BC2 results have low power and unstable variance estimates")
prv_cut <- as.numeric(opt$`prv-cut`); lib_cut <- as.numeric(opt$`lib-cut`)
# ANCOM-BC2 filters features by prevalence, then samples by library size (lib_cut): check the groups that remain
keep_f <- rowMeans(mat > 0) >= prv_cut
lib_after <- colSums(mat[keep_f, , drop = FALSE])
low_lib <- names(lib_after)[lib_after < lib_cut]
if (length(low_lib)) {
  n_after <- table(droplevels(meta[setdiff(colnames(mat), low_lib), group]))
  if (length(n_after) < 2 || any(n_after < 2) || length(n_after) < nlevels(meta[[group]]))
    stop("after --lib-cut ", lib_cut, " (samples ", paste(low_lib, collapse = ", "), " fall below it) every group needs ",
         "at least 2 retained samples; found ", paste(names(n_after), n_after, sep = "=", collapse = ", "))
}
kept <- sum(keep_f)
if (kept < 50) note(kept, " features pass the prevalence filter (fewer than 50): ANCOM-BC2 estimates the sample-specific ",
                    "bias from the features, so with a small focused panel the bias correction is unreliable; ",
                    "read the log fold changes as relative to the other panel features")
if (kept < 2) stop("fewer than 2 features pass the prevalence filter; nothing to test")

# random effects (repeated measures): check the grouping columns and whether the model can be estimated
rand <- opt$`rand-formula`
rand_note <- "none"
if (!is.null(rand)) {
  rand <- gsub("\\s+", "", rand)
  terms_r <- regmatches(rand, gregexpr("\\([^()]*\\|[^()]*\\)", rand))[[1]]
  if (!length(terms_r) || nchar(gsub("[+]", "", paste(terms_r, collapse = ""))) != nchar(gsub("[+]", "", rand)))
    stop("--rand-formula must look like '(1|unit)' or '(1|unit)+(1|batch)', not '", opt$`rand-formula`, "'")
  unit_cols <- unique(unlist(strsplit(sub("^.*\\|", "", gsub("[()]", "", terms_r)), "[:/]")))
  miss_u <- setdiff(unit_cols, names(meta))
  if (length(miss_u)) stop("--rand-formula uses column(s) not in the metadata: ", paste(miss_u, collapse = ", "))
  if (anyNA(meta[unit_cols]) || any(meta[unit_cols] == "")) stop("empty values in the random-effect column(s) ",
                                                                 paste(unit_cols, collapse = ", "))
  parts <- character(0)
  for (u in unit_cols) {
    meta[[u]] <- factor(meta[[u]])
    per_unit <- table(meta[[u]])
    units_per_group <- tapply(as.character(meta[[u]]), meta[[group]], function(x) length(unique(x)))
    parts <- c(parts, paste0(u, ": ", nlevels(meta[[u]]), " units, ", min(per_unit), "-", max(per_unit),
                             " samples per unit, units per group ",
                             paste(names(units_per_group), units_per_group, sep = "=", collapse = ", ")))
    if (nlevels(meta[[u]]) >= nrow(meta))
      stop("random effect '", u, "' has one sample per unit (", nlevels(meta[[u]]), " units, ", nrow(meta),
           " samples): a unit effect cannot be separated from the residual; remove --rand-formula or give the ",
           "column that groups the repeated samples")
    if (nlevels(meta[[u]]) < 2) stop("random effect '", u, "' has a single unit")
    if (any(per_unit < 2)) note("random effect '", u, "': ", sum(per_unit < 2), " of ", length(per_unit),
                                " units have one sample; they do not help to estimate the unit variance")
    if (any(units_per_group < 3)) note("random effect '", u, "': fewer than 3 units in a group (",
                                       paste(names(units_per_group), units_per_group, sep = "=", collapse = ", "),
                                       "); the unit variance is poorly estimated and the group test rests on very few independent units")
  }
  rand_note <- paste0(rand, " [", paste(parts, collapse = "; "), "]")
  message("random effects: ", rand_note)
  if (isTRUE(opt$`lme-relax-checks`))
    note("--lme-relax-checks: lme4 checks on the number of observations per random-effect level are switched off; ",
         "features observed in few samples get fits whose unit variance is not identifiable (raise --prv-cut instead where possible)")
} else if (isTRUE(opt$`lme-relax-checks`)) stop("--lme-relax-checks needs --rand-formula")

fit_args <- list(data = mat, meta_data = meta, tax_level = NULL, fix_formula = formula, rand_formula = rand,
                 p_adj_method = opt$`p-adjust`, pseudo_sens = TRUE, prv_cut = prv_cut, lib_cut = lib_cut,
                 s0_perc = 0.05, group = group, struc_zero = isTRUE(opt$`struc-zero`), neg_lb = isTRUE(opt$`neg-lb`),
                 alpha = as.numeric(opt$alpha), n_cl = as.integer(opt$threads), verbose = FALSE,
                 global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE)
if (!is.null(rand) && isTRUE(opt$`lme-relax-checks`))
  fit_args$lme_control <- lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
fit_args <- fit_args[names(fit_args) %in% names(formals(ancombc2))]
set.seed(11)
fit <- withCallingHandlers(do.call(ancombc2, fit_args),
                           warning = function(w) { warn <<- c(warn, one_line(paste0("ANCOM-BC2: ", conditionMessage(w))))
                                                   invokeRestart("muffleWarning") })
# validate the design on the samples ANCOM-BC2 actually kept (lib_cut and the prevalence filter can drop some)
kept_samples <- if (!is.null(fit$feature_table)) colnames(fit$feature_table) else colnames(mat)
dropped_samples <- setdiff(colnames(mat), kept_samples)
meta_k <- meta[kept_samples, , drop = FALSE]
meta_k[[group]] <- droplevels(meta_k[[group]])
n_group_k <- table(meta_k[[group]])
if (length(dropped_samples)) note(length(dropped_samples), " sample(s) removed by ANCOM-BC2's filters (lib_cut / prevalence): ",
                                  paste(dropped_samples, collapse = ", "))
if (nlevels(meta_k[[group]]) < 2 || any(n_group_k < 2))
  stop("after ANCOM-BC2's sample filters the groups are ", paste(names(n_group_k), n_group_k, sep = "=", collapse = ", "),
       "; every group needs at least 2 retained samples (lower --lib-cut)")
if (any(n_group_k < 5) && !all(n_group_k == n_group[names(n_group_k)]))
  note("after ANCOM-BC2's sample filters fewer than 5 samples remain in a group (",
       paste(names(n_group_k), n_group_k, sep = "=", collapse = ", "), ")")
if (!is.null(rand)) {
  for (u in unit_cols) {
    upg <- tapply(as.character(meta_k[[u]]), meta_k[[group]], function(x) length(unique(x)))
    if (length(unique(meta_k[[u]])) >= nrow(meta_k))
      stop("after ANCOM-BC2's sample filters random effect '", u, "' has one sample per unit")
    if (any(upg < 3)) note("after ANCOM-BC2's sample filters, random effect '", u, "' has fewer than 3 units in a group (",
                           paste(names(upg), upg, sep = "=", collapse = ", "), ")")
  }
}
res <- fit$res
names(res)[1] <- "feature_id"
for (col in c("label", "identity_class")) if (!col %in% names(feats)) feats[[col]] <- NA_character_
res$label <- feats$label[match(res$feature_id, feats$feature_id)]
res$identity_class <- feats$identity_class[match(res$feature_id, feats$feature_id)]
# sensitivity-robust call for releases that only report diff_* and passed_ss_*
for (d in grep("^diff_", names(res), value = TRUE)) {
  term <- sub("^diff_", "", d)
  if (startsWith(d, "diff_robust_")) next
  ss <- paste0("passed_ss_", term)
  if (!paste0("diff_robust_", term) %in% names(res) && ss %in% names(res)) {
    res[[paste0("diff_robust_", term)]] <- res[[d]] & res[[ss]]
  }
}
res <- res[, c("feature_id", "label", "identity_class", setdiff(names(res), c("feature_id", "label", "identity_class")))]
# what each coefficient is: a main-effect group level against the reference level, an interaction, or a covariate
ref_level <- levels(meta[[group]])[1]
coef_label <- function(t) {
  if (grepl(":", t, fixed = TRUE)) return(paste0("interaction ", t))
  lv <- sub(paste0("^", group), "", t)
  if (startsWith(t, group) && lv %in% levels(meta[[group]])) return(paste0(lv, " vs ", ref_level))
  paste0("covariate ", t)
}
all_terms <- setdiff(sub("^lfc_", "", grep("^lfc_", names(res), value = TRUE)), "(Intercept)")
coef_line <- paste0("# coefficients: ", paste(all_terms, sapply(all_terms, coef_label), sep = " = ", collapse = "; "))
version <- as.character(packageVersion("ANCOMBC"))
header <- c(paste0("# ANCOM-BC2, ANCOMBC ", version, ", ", R.version.string),
            paste0("# rand_formula: ", rand_note, if (isTRUE(opt$`lme-relax-checks`)) "; lme checks relaxed" else ""),
            paste0("# formula: ", formula, "; group: ", group, " (reference level: ", levels(meta[[group]])[1], "); samples per group: ",
                   paste(names(n_group), n_group, sep = "=", collapse = ", "), "; retained after ANCOM-BC2's filters: ",
                   paste(names(n_group_k), n_group_k, sep = "=", collapse = ", ")),
            paste0("# features in: ", nrow(mat), "; tested: ", nrow(res), "; prv_cut=", prv_cut, "; lib_cut=", lib_cut,
                   "; struc_zero=", isTRUE(opt$`struc-zero`), "; neg_lb=", isTRUE(opt$`neg-lb`),
                   "; pseudo_sens=TRUE; p_adj_method=", opt$`p-adjust`,
                   "; alpha=", opt$alpha),
            coef_line,
            "# lfc_* are natural-log fold changes; use diff_robust_* (significant and stable to the pseudo-count) as the call",
            if (length(warn)) paste0("# WARNING: ", warn) else "# warnings: none")
sz_note <- "# structural zeros: not tested (--struc-zero off)"
if (isTRUE(opt$`struc-zero`)) {
  zi <- fit$zero_ind
  sz_file <- file.path(opt$outdir, "structural_zeros.tsv")
  if (!is.null(zi) && nrow(zi)) {
    names(zi)[1] <- "feature_id"
    flag_cols <- setdiff(names(zi), "feature_id")
    zi$label <- feats$label[match(zi$feature_id, feats$feature_id)]
    zi$flagged_by_ancombc2 <- apply(zi[, flag_cols, drop = FALSE], 1, any)
    # 'absent' is a statement about the data: zero reads in every retained sample of the group
    grp_k <- meta_k[[group]]
    absent <- sapply(zi$feature_id, function(f) {
      if (!f %in% rownames(mat)) return("")
      g0 <- tapply(mat[f, kept_samples], grp_k, function(x) all(x == 0))
      paste(names(g0)[g0], collapse = ",")
    })
    zi$absent_in_groups <- absent
    zi$absent_in_a_group <- absent != ""
    zi <- zi[, c("feature_id", "label", flag_cols, "flagged_by_ancombc2", "absent_in_groups", "absent_in_a_group")]
    # the evidence behind the decision, per group level and on the same retained samples: presence is count != 0,
    # the rule ANCOM-BC2 applies; reads keep their fractions (EM expected counts), nothing is rounded or cut
    rows_f <- match(zi$feature_id, rownames(mat))
    for (lv in levels(grp_k)) {
      s_lv <- kept_samples[grp_k == lv]
      sub <- mat[rows_f, s_lv, drop = FALSE]
      tag <- paste0(" (", group, " = ", lv, ")")
      zi[[paste0("samples", tag)]] <- length(s_lv)
      zi[[paste0("present_samples", tag)]] <- rowSums(sub != 0)
      zi[[paste0("reads", tag)]] <- rowSums(sub)
      zi[[paste0("max_count", tag)]] <- apply(sub, 1, max)
    }
    write.table(zi, sz_file, sep = "\t", quote = FALSE, row.names = FALSE)
    n_sz <- sum(zi$absent_in_a_group)
    by_grp <- table(unlist(strsplit(absent[absent != ""], ",")))
    sz_note <- paste0("# structural zeros: ", n_sz, " feature(s) absent from every sample of a group",
                      if (length(by_grp)) paste0(" (", paste(names(by_grp), by_grp, sep = ": ", collapse = ", "), ")") else "",
                      "; ", sum(zi$flagged_by_ancombc2), " flagged by ANCOM-BC2",
                      if (isTRUE(opt$`neg-lb`)) " (neg_lb on: also features with a negative lower bound)" else "",
                      "; flagged features are not tested here, see structural_zeros.tsv; presence follows ",
                      "ANCOM-BC2's rule (count != 0 in a sample), so with EM expected counts a fraction of a shared ",
                      "read counts as present (structural_zeros.tsv gives present samples, reads and the largest ",
                      "count per group; `unclassified-meta detection` sets unique-read, LCA and EM evidence side by ",
                      "side); --neg-lb is a sensitivity analysis")
  } else {
    write.table(data.frame(feature_id = character(0)), sz_file, sep = "\t", quote = FALSE, row.names = FALSE)
    sz_note <- "# structural zeros: none (presence follows ANCOM-BC2's rule, count != 0 in a sample)"
  }
}
header <- append(header, sz_note, after = length(header) - 1)
out <- file.path(opt$outdir, "ancombc2_results.tsv")
writeLines(header, out)
suppressWarnings(write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
writeLines(c(header, "", capture.output(sessionInfo())), file.path(opt$outdir, "session_info.txt"))
writeLines(if (length(warn)) warn else "none", file.path(opt$outdir, "warnings.txt"))

terms <- sub("^lfc_", "", grep("^lfc_", names(res), value = TRUE))
terms <- terms[terms != "(Intercept)" & startsWith(terms, group)]   # group main effects and their interactions
save_plot <- function(p, name, w, h) {
  ggsave(file.path(opt$outdir, paste0(name, ".png")), p, width = w, height = h, dpi = 200)
  ggsave(file.path(opt$outdir, paste0(name, ".pdf")), p, width = w, height = h)
}
if (length(terms)) {
  long <- do.call(rbind, lapply(terms, function(t) {
    robust <- res[[paste0("diff_robust_", t)]]
    data.frame(feature = ifelse(is.na(res$label), res$feature_id, res$label),
               contrast = coef_label(t),
               lfc = res[[paste0("lfc_", t)]] / log(2), q = res[[paste0("q_", t)]],
               call = ifelse(!is.na(robust) & robust, "significant, robust",
                             ifelse(res[[paste0("diff_", t)]], "significant, pseudo-count sensitive", "not significant")))
  }))
  long$mlq <- -log10(pmax(long$q, 1e-300))
  p <- ggplot(long, aes(lfc, mlq, colour = call)) + geom_point(size = 2.2) +
    geom_hline(yintercept = -log10(as.numeric(opt$alpha)), linetype = 2, colour = "grey50") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_text(data = long[long$call != "not significant", ], aes(label = feature), size = 2.6, vjust = -0.8,
              hjust = "inward", show.legend = FALSE, check_overlap = TRUE) +
    scale_x_continuous(expand = expansion(mult = 0.12)) + scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
    scale_colour_manual(values = c("significant, robust" = "#E64B35", "significant, pseudo-count sensitive" = "#F39B7F",
                                   "not significant" = "#8C8C8C")) +
    facet_wrap(~contrast) + theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(x = "log2 coefficient (ANCOM-BC2, bias-corrected; see the panel title for what it compares)",
         y = "-log10 adjusted P")
  save_plot(p, "volcano", 6.5, 5)
}
ord <- order(apply(res[, paste0("q_", terms), drop = FALSE], 1, min, na.rm = TRUE))
top <- res$feature_id[head(ord, 12)]
lib <- colSums(mat)
ylab <- "share of the fragments assigned to panel features (log scale)"
if (!is.null(opt$ledger) && file.exists(opt$ledger)) {
  led <- read.delim(opt$ledger, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c(sample = "character"),
                    na.strings = character(0))
  if ("input_fragments" %in% names(led) && all(colnames(mat) %in% led$sample)) {
    lib <- setNames(led$input_fragments[match(colnames(mat), led$sample)], colnames(mat))
    ylab <- "share of all fragments of the sample (log scale)"
  }
}
ab <- do.call(rbind, lapply(top, function(f) data.frame(
  feature = ifelse(is.na(feats$label[match(f, feats$feature_id)]), f, feats$label[match(f, feats$feature_id)]),
  sample = colnames(mat), group = meta[[group]], share = (mat[f, ] + 0.5) / (lib + 0.5))))
ab$feature <- factor(ab$feature, levels = unique(ab$feature))
p2 <- ggplot(ab, aes(group, share, colour = group)) + geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, size = 1.4) + scale_y_log10() + facet_wrap(~feature, scales = "free_y", ncol = 4,
                                                                   labeller = label_wrap_gen(28)) +
  theme_bw(base_size = 9) + theme(legend.position = "none", strip.text = element_text(size = 7)) +
  labs(x = NULL, y = ylab)
save_plot(p2, "feature_abundance", 9, 2.4 * ceiling(length(top) / 4) + 0.6)
message("results: ", out)
