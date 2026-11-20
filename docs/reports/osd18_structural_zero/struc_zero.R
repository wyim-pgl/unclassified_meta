# ANCOM-BC2's own structural-zero rule (ANCOMBC:::.get_struc_zero, neg_lb FALSE and TRUE) on the saved count tables
suppressMessages(library(ANCOMBC))
args <- commandArgs(TRUE); O <- args[1]; R <- args[2]; outp <- args[3]
meta <- read.delim(file.path(dirname(O), "meta_RNA.tsv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
rd <- function(f) { x <- read.delim(f, comment.char = "#", check.names = FALSE, quote = ""); m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m[, meta$sample] }
src <- list(lca = file.path(O, sprintf("taxa_%s_lca/merged/counts.tsv", R)), em_expected = file.path(O, sprintf("taxa_%s_em/merged/counts.tsv", R)))
imp <- sort(Sys.glob(file.path(O, sprintf("taxa_%s_em/merged/imputations/counts_*.tsv", R))))
names(imp) <- sprintf("imp_%03d", seq_along(imp)); src <- c(src, as.list(imp))
cat("ANCOMBC", as.character(packageVersion("ANCOMBC")), "\n")
res <- NULL
for (nm in names(src)) for (nl in c(FALSE, TRUE)) {
  m <- rd(src[[nm]])
  z <- ANCOMBC:::.get_struc_zero(data = m, meta_data = meta, group = "condition", neg_lb = nl)
  z$any <- apply(z[, -1], 1, any)
  colnames(z) <- c("feature_id", "sz_healthy", "sz_macho", "sz_any")
  z$source <- nm; z$neg_lb <- nl
  res <- rbind(res, z)
  if (!nl) cat(nm, ": features", nrow(z), " struc-zero healthy", sum(z$sz_healthy), " macho", sum(z$sz_macho), "\n")
}
write.table(res, outp, sep = "\t", quote = FALSE, row.names = FALSE)
