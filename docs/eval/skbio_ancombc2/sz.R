# R structural zeros (ANCOMBC:::.get_struc_zero, the function ancombc2(struc_zero = TRUE) calls) on a count table
# usage: Rscript sz.R counts.tsv meta.tsv out_prefix
a <- commandArgs(trailingOnly = TRUE)
counts <- read.delim(a[1], check.names = FALSE, colClasses = c(feature_id = "character"))
meta <- read.delim(a[2], check.names = FALSE, colClasses = "character")
mat <- as.matrix(counts[, meta$sample]); rownames(mat) <- counts$feature_id
rownames(meta) <- meta$sample
for (nl in c(FALSE, TRUE)) {
  z <- ANCOMBC:::.get_struc_zero(data = mat, meta_data = meta, group = "condition", neg_lb = nl)
  names(z) <- c("feature_id", "healthy", "macho")
  write.table(z, paste0(a[3], if (nl) "_neglb.tsv" else ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
