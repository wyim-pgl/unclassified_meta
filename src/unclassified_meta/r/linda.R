#!/usr/bin/env Rscript
# Differential abundance with LinDA (MicrobiomeStat::linda) for `unclassified-meta da-robustness`.
#
# Usage: Rscript linda.R --counts counts.tsv --metadata metadata.tsv --group group --outdir out \
#          [--formula "group + batch"] [--reference-level control] [--rand-formula "(1|unit)"]
#          [--alpha 0.05] [--p-adjust holm] [--threads 1] [--seed 1]
#
# Writes outdir/results.tsv (feature_id, term, effect, se, ci_low, ci_high, stat, p, q, effect_definition),
# session_info.txt and warnings.txt.
#
# What the numbers are: LinDA fits a linear (or, with --rand-formula, linear mixed) model to the CLR of the counts
# (zeros handled with LinDA's pseudo-count, adaptive), then removes the compositional bias from the coefficients
# (mode of the coefficient distribution). effect = log2FoldChange (the bias-corrected CLR-scale coefficient, log2),
# se = lfcSE, stat, p = pvalue and q = padj with --p-adjust, all as LinDA reports them; ci = effect +- t(0.975, df) se.
# LinDA's own filters are switched off (prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0):
# da-robustness filters the features once, in Python, for all methods. Fractional (EM expected) counts are
# accepted as they are (feature.dat.type = "count"); LinDA's winsorization of outliers (is.winsor, default) stays on.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(formula = NULL, `reference-level` = NULL, alpha = "0.05", `p-adjust` = "holm", threads = "1",
            outdir = "linda", `rand-formula` = NULL, seed = "1")
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (i == length(args)) stop("option --", key, " needs a value")
  opt[[key]] <- args[i + 1]; i <- i + 2
}
for (need in c("counts", "metadata", "group")) {
  if (is.null(opt[[need]])) stop("missing required option --", need)
}
suppressPackageStartupMessages(library(MicrobiomeStat))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
warn <- character(0)
one_line <- function(x) gsub("\\s+", " ", trimws(x))
note <- function(...) { m <- one_line(paste0(...)); warn <<- c(warn, m); message("WARNING: ", m) }

# identifiers are text ("01" stays "01"); a sample called "NA" stays a name (as in ancombc2.R)
counts <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c(feature_id = "character"))
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

formula <- if (is.null(opt$formula)) group else opt$formula
fvars <- tryCatch(all.vars(as.formula(paste("~", formula))), error = function(e) stop("--formula '", formula,
                  "' is not a valid model formula: ", conditionMessage(e)))
miss_f <- setdiff(fvars, names(meta))
if (length(miss_f)) stop("--formula uses column(s) not in the metadata: ", paste(miss_f, collapse = ", "))
if (!group %in% fvars) stop("the group variable '", group, "' is not in --formula '", formula, "'")
rand <- opt$`rand-formula`
rvars <- character(0)
if (!is.null(rand)) {
  rand <- gsub("\\s+", "", rand)
  if (!grepl("^\\([^()]*\\|[^()]*\\)(\\+\\([^()]*\\|[^()]*\\))*$", rand))
    stop("--rand-formula must look like '(1|unit)' or '(1|unit)+(1|batch)', not '", opt$`rand-formula`, "'")
  rvars <- all.vars(as.formula(paste("~", gsub("\\|", "+", rand))))
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
# the same variable coding as ancombc2.R, so the model-matrix column (term) names agree across the methods
for (v in setdiff(fvars, group)) {
  num <- suppressWarnings(as.numeric(meta[[v]]))
  meta[[v]] <- if (!anyNA(num)) num else factor(meta[[v]])
}
for (u in setdiff(rvars, fvars)) meta[[u]] <- factor(meta[[u]])
meta[[group]] <- factor(meta[[group]])
if (!is.null(opt$`reference-level`)) {
  if (!opt$`reference-level` %in% levels(meta[[group]])) stop("reference level not in column ", group)
  meta[[group]] <- relevel(meta[[group]], ref = opt$`reference-level`)
}
n_group <- table(meta[[group]])
if (nlevels(meta[[group]]) < 2) stop("column ", group, " needs at least two groups")
if (any(n_group < 2)) stop("every group needs at least 2 samples; found ",
                           paste(names(n_group), n_group, sep = "=", collapse = ", "))
for (u in setdiff(rvars, fvars)) {
  if (nlevels(meta[[u]]) >= nrow(meta))
    stop("random effect '", u, "' has one sample per unit: a unit effect cannot be separated from the residual")
  if (nlevels(meta[[u]]) < 2) stop("random effect '", u, "' has a single unit")
}
if (nrow(mat) < 2) stop("fewer than 2 features; nothing to test")

fml <- paste("~", formula, if (!is.null(rand)) paste("+", gsub("\\)\\+\\(", ") + (", rand)) else "")
lin_args <- list(feature.dat = mat, meta.dat = meta, formula = fml, feature.dat.type = "count",
                 prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0,
                 alpha = as.numeric(opt$alpha), p.adj.method = opt$`p-adjust`, n.cores = as.integer(opt$threads),
                 verbose = FALSE)
lin_args <- lin_args[names(lin_args) %in% names(formals(MicrobiomeStat::linda))]
set.seed(as.integer(opt$seed))
fit <- withCallingHandlers(suppressMessages(do.call(MicrobiomeStat::linda, lin_args)),
                           warning = function(w) { warn <<- c(warn, one_line(paste0("LinDA: ", conditionMessage(w))))
                                                   invokeRestart("muffleWarning") })
mm_terms <- setdiff(colnames(model.matrix(as.formula(paste("~", formula)), meta)), "(Intercept)")
def <- "LinDA bias-corrected CLR-scale coefficient (log2)"
rows <- list()
for (t in names(fit$output)) {
  d <- fit$output[[t]]
  if (!t %in% mm_terms) next
  df_t <- if ("df" %in% names(d)) d$df else rep(Inf, nrow(d))
  tq <- qt(0.975, df_t)
  rows[[t]] <- data.frame(feature_id = rownames(d), term = t, effect = d$log2FoldChange, se = d$lfcSE,
                          ci_low = d$log2FoldChange - tq * d$lfcSE, ci_high = d$log2FoldChange + tq * d$lfcSE,
                          stat = d$stat, p = d$pvalue, q = d$padj, effect_definition = def,
                          check.names = FALSE, stringsAsFactors = FALSE)
}
missing_terms <- setdiff(mm_terms, names(rows))
if (length(missing_terms)) note("LinDA reported no coefficient for term(s) ", paste(missing_terms, collapse = ", "))
if (!length(rows)) stop("LinDA returned no coefficients for the model terms ", paste(mm_terms, collapse = ", "))
res <- do.call(rbind, rows)
write.table(res, file.path(opt$outdir, "results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
header <- c(paste0("# LinDA, MicrobiomeStat ", as.character(packageVersion("MicrobiomeStat")), ", ", R.version.string),
            paste0("# formula: ", fml, "; group: ", group, " (reference level: ", levels(meta[[group]])[1],
                   "); samples per group: ", paste(names(n_group), n_group, sep = "=", collapse = ", ")),
            paste0("# features: ", nrow(mat), "; prev.filter=0; mean.abund.filter=0; p.adj.method=", opt$`p-adjust`,
                   "; alpha=", opt$alpha))
writeLines(c(header, "", capture.output(sessionInfo())), file.path(opt$outdir, "session_info.txt"))
writeLines(if (length(warn)) warn else "none", file.path(opt$outdir, "warnings.txt"))
message("results: ", file.path(opt$outdir, "results.tsv"))
