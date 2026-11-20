#!/usr/bin/env Rscript
# Differential abundance with ALDEx2 (aldex.clr + aldex.glm) for `unclassified-meta da-robustness`.
#
# Usage: Rscript aldex2.R --counts counts.tsv --metadata metadata.tsv --group group --outdir out \
#          [--formula "group + batch"] [--reference-level control] [--alpha 0.05] [--p-adjust BH]
#          [--mc-samples 128] [--seed 1] [--threads 1] [--rand-formula "(1|unit)"]
#
# Writes outdir/results.tsv (feature_id, term, effect, se, ci_low, ci_high, stat, p, q, effect_definition),
# session_info.txt and warnings.txt.
#
# What the numbers are: aldex.clr draws --mc-samples Dirichlet Monte-Carlo instances of each sample's composition
# (the sampling of reads given an INTEGER table), takes the centred log-ratio (log2, denom = "all") of each, and
# aldex.glm fits the model-matrix glm to every instance. effect / se / stat / p are ALDEx2's means over the
# instances (p is the mean of the per-instance glm t-test p-values, as ALDEx2 reports it; it is not a pooled test).
# q = p.adjust(p, --p-adjust) within each term.
# Limits, written to warnings.txt when they apply:
#  - counts are rounded to integers (EM expected counts are fractional; the Dirichlet model needs integers);
#    read-assignment uncertainty therefore does not enter ALDEx2 itself, only through the bootstrap / realization
#    tables that da-robustness refits;
#  - aldex.glm has no random effects: --rand-formula is ignored and the fixed-effects formula alone is fitted
#    (repeated samples of one unit are then treated as independent; the report says so);
#  - no prevalence filter here: da-robustness filters the features in Python, once, for all methods.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(formula = NULL, `reference-level` = NULL, alpha = "0.05", `p-adjust` = "BH", threads = "1",
            outdir = "aldex2", `rand-formula` = NULL, seed = "1", `mc-samples` = "128")
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (i == length(args)) stop("option --", key, " needs a value")
  opt[[key]] <- args[i + 1]; i <- i + 2
}
for (need in c("counts", "metadata", "group")) {
  if (is.null(opt[[need]])) stop("missing required option --", need)
}
suppressPackageStartupMessages(library(ALDEx2))
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
if (!is.null(opt$`rand-formula`))
  note("ALDEx2 (aldex.glm) has no random effects: --rand-formula '", opt$`rand-formula`, "' was ignored and only the ",
       "fixed effects (", formula, ") were fitted; repeated samples of one unit are treated as independent")
is_missing <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "NA")
na_rows <- Reduce(`|`, lapply(fvars, function(v) is_missing(meta[[v]])))
if (any(na_rows)) {
  note(sum(na_rows), " sample(s) with missing values in ", paste(fvars, collapse = ", "), " were left out: ",
       paste(meta$sample[na_rows], collapse = ", "))
  meta <- meta[!na_rows, , drop = FALSE]
  mat <- mat[, meta$sample, drop = FALSE]
}
# the same variable coding as ancombc2.R, so the model-matrix column (term) names agree across the methods
for (v in setdiff(fvars, group)) {
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

if (any(abs(mat - round(mat)) > 1e-9))
  note("counts were rounded to integers for ALDEx2 (EM expected counts are fractional; the Dirichlet Monte Carlo ",
       "models read sampling of an integer table and not the read-assignment uncertainty)")
reads <- round(mat)
keep <- rowSums(reads) > 0
if (any(!keep)) note(sum(!keep), " feature(s) with no reads after rounding were not tested: ",
                     paste(head(rownames(reads)[!keep], 10), collapse = ", "), if (sum(!keep) > 10) ", ..." else "")
reads <- reads[keep, , drop = FALSE]
if (nrow(reads) < 2) stop("fewer than 2 features with reads; nothing to test")
storage.mode(reads) <- "integer"
mm <- model.matrix(as.formula(paste("~", formula)), meta)
mc <- as.integer(opt$`mc-samples`)
set.seed(as.integer(opt$seed))
catch <- function(expr) withCallingHandlers(expr, warning = function(w) {
  warn <<- c(warn, one_line(paste0("ALDEx2: ", conditionMessage(w)))); invokeRestart("muffleWarning") })
clr <- catch(aldex.clr(as.data.frame(reads), mm, mc.samples = mc, denom = "all", verbose = FALSE,
                       useMC = as.integer(opt$threads) > 1))
# a data frame that keeps aldex.glm's (partly duplicated) column names; [[ ]] takes the first of a name
glm_res <- as.data.frame(catch(aldex.glm(clr, mm)), check.names = FALSE)

# aldex.glm column names differ between releases: "term:Est", "term:SE", "term:t.val", "term:pval" in 1.42 (that
# "term:pval" is ALDEx2's two-sided p from the Monte-Carlo-averaged one-sided p-values; "term:pval.padj" is its own
# adjustment, not used), "term::Est" in some releases, "model.term Estimate" in older ones
pick <- function(term, new, old) {
  for (nm in c(paste0(term, ":", new), paste0(term, "::", new), paste0("model.", term, " ", old)))
    if (nm %in% names(glm_res)) return(glm_res[[nm]])
  NULL
}
terms <- setdiff(colnames(mm), "(Intercept)")
alpha <- as.numeric(opt$alpha)
def <- "ALDEx2 aldex.glm coefficient on the CLR scale (log2; mean over Monte-Carlo Dirichlet instances)"
rows <- list()
for (t in terms) {
  est <- pick(t, "Est", "Estimate"); se <- pick(t, "SE", "Std. Error")
  st <- pick(t, "t.val", "t value"); p <- pick(t, "pval", "Pr(>|t|)")
  if (is.null(est) || is.null(p)) {
    note("aldex.glm reported no coefficient for term '", t, "' (columns: ", paste(head(names(glm_res), 12), collapse = ", "), ")")
    next
  }
  if (is.null(se)) se <- rep(NA_real_, length(est))
  if (is.null(st)) st <- rep(NA_real_, length(est))
  rows[[t]] <- data.frame(feature_id = rownames(glm_res), term = t, effect = est, se = se,
                          ci_low = est - qnorm(0.975) * se, ci_high = est + qnorm(0.975) * se,
                          stat = st, p = p, q = p.adjust(p, method = opt$`p-adjust`), effect_definition = def,
                          check.names = FALSE, stringsAsFactors = FALSE)
}
if (!length(rows)) stop("aldex.glm returned no coefficients for the model terms ", paste(terms, collapse = ", "))
res <- do.call(rbind, rows)
write.table(res, file.path(opt$outdir, "results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
header <- c(paste0("# ALDEx2 ", as.character(packageVersion("ALDEx2")), ", ", R.version.string),
            paste0("# formula: ~ ", formula, "; group: ", group, " (reference level: ", levels(meta[[group]])[1],
                   "); samples per group: ", paste(names(n_group), n_group, sep = "=", collapse = ", ")),
            paste0("# features tested: ", nrow(reads), "; mc.samples=", mc, "; denom=all; seed=", opt$seed,
                   "; p.adjust=", opt$`p-adjust`, "; alpha=", alpha))
writeLines(c(header, "", capture.output(sessionInfo())), file.path(opt$outdir, "session_info.txt"))
writeLines(if (length(warn)) warn else "none", file.path(opt$outdir, "warnings.txt"))
message("results: ", file.path(opt$outdir, "results.tsv"))
