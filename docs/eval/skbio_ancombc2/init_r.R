# R initial ANCOM-BC2 estimates (the .ancombc2_core steps before bias EM) on prv_cut-filtered counts, pseudo 0.
# usage: Rscript init_r.R counts.tsv meta.tsv out_prefix
suppressPackageStartupMessages(library(ANCOMBC))
a <- commandArgs(trailingOnly = TRUE)
counts <- read.delim(a[1], check.names = FALSE, colClasses = c(feature_id = "character"))
meta <- read.delim(a[2], check.names = FALSE, colClasses = "character")
mat <- as.matrix(counts[, meta$sample]); rownames(mat) <- counts$feature_id
mat <- mat[rowMeans(mat != 0) >= 0.10, , drop = FALSE]
rownames(meta) <- meta$sample
meta$condition <- relevel(factor(meta$condition), ref = "healthy"); meta$plant <- factor(meta$plant)
fix_formula <- "condition + plant"
options(na.action = "na.pass"); x <- model.matrix(formula(paste0("~", fix_formula)), data = meta)
options(na.action = "na.omit")
o1 <- log(mat); o1[is.infinite(o1)] <- NA
y1 <- o1 - rowMeans(o1, na.rm = TRUE)
set.seed(11)
p1 <- ANCOMBC:::.iter_mle(x = x, y = y1, meta_data = meta, formula = fix_formula, theta = NULL, tol = 0.01, max_iter = 20)
bias <- t(sapply(seq_len(ncol(p1$beta)), function(i)
  ANCOMBC:::.bias_em(beta = p1$beta[, i], var_hat = p1$var_hat[, i], tol = 1e-5, max_iter = 100)))
rownames(bias) <- colnames(x)
write.table(data.frame(feature_id = rownames(p1$beta), p1$beta, check.names = FALSE), paste0(a[3], "_beta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(feature_id = rownames(p1$var_hat), p1$var_hat, check.names = FALSE), paste0(a[3], "_var.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(term = rownames(bias), bias, check.names = FALSE), paste0(a[3], "_bias.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(as.character(p1$theta), paste0(a[3], "_theta.txt"))
