# 36_target_evidence_matrix.R
# Build an auditable, sample-level target evidence matrix for the fixed 200-gene
# aggressiveness programme. This is prioritisation evidence, not target validation.
suppressMessages({library(Seurat); library(Matrix)})

root <- normalizePath(".")
outdir <- file.path(root, "results", "drug_repurposing")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

program <- read.csv(file.path(root, "docs", "Table_S2_aggressiveness_program_genes.csv"),
                    stringsAsFactors=FALSE, check.names=FALSE)
obj <- readRDS(file.path(root, "results", "scrna", "GSE206647_processed.rds"))
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data
stopifnot(all(c("celltype", "grade", "gsm") %in% colnames(md)))

tum.cells <- rownames(md)[md$celltype %in% c("Meningioma", "Tumor", "Tumour")]
tum <- subset(obj, cells=tum.cells)
counts <- GetAssayData(tum, assay="RNA", layer="counts")
genes <- intersect(program$symbol, rownames(counts))

# Patient-level pseudobulk avoids treating cells from one tumour as replicates.
sample.cells <- split(colnames(counts), tum$gsm)
# Library sizes must be computed from the full transcriptome, not only the 200
# programme genes. The former implementation subset first and inflated gene-level
# grade associations; the independent submission audit detected this error.
lib.size <- vapply(sample.cells, function(cells) sum(counts[, cells, drop=FALSE]), numeric(1))
pb <- sapply(sample.cells, function(cells) Matrix::rowSums(counts[genes, cells, drop=FALSE]))
logcpm <- log1p(t(t(pb) / lib.size) * 1e6)
sample.grade <- vapply(sample.cells, function(cells) {
  values <- as.character(tum@meta.data[cells, "grade"])
  names(sort(table(values), decreasing=TRUE))[1]
}, character(1))
grade.ord <- setNames(c(1,2,3), c("I","II","III"))[sample.grade]

rho <- pval <- rep(NA_real_, nrow(logcpm)); names(rho) <- names(pval) <- rownames(logcpm)
for (gene in rownames(logcpm)) {
  keep <- !is.na(grade.ord) & is.finite(logcpm[gene, ])
  if (sum(keep) >= 6 && length(unique(logcpm[gene, keep])) > 1) {
    test <- suppressWarnings(cor.test(logcpm[gene, keep], grade.ord[keep], method="spearman", exact=FALSE))
    rho[gene] <- unname(test$estimate); pval[gene] <- test$p.value
  }
}
fdr <- p.adjust(pval, method="BH")

# Tumour specificity is descriptive: log2 ratio of mean expression and detection
# frequency in meningioma cells versus all annotated non-tumour cells.
all.counts <- GetAssayData(obj, assay="RNA", layer="counts")[genes, , drop=FALSE]
is.tum <- colnames(all.counts) %in% tum.cells
mean.tum <- Matrix::rowMeans(all.counts[, is.tum, drop=FALSE])
mean.other <- Matrix::rowMeans(all.counts[, !is.tum, drop=FALSE])
pct.tum <- Matrix::rowMeans(all.counts[, is.tum, drop=FALSE] > 0)
pct.other <- Matrix::rowMeans(all.counts[, !is.tum, drop=FALSE] > 0)
specificity <- log2((mean.tum + 0.01) / (mean.other + 0.01))

# Independent-cohort recovery flags reuse only the saved gene sets, not old scores.
v1 <- readRDS(file.path(root, "results", "deg", "GSE16581_program_validation.rds"))
v2 <- readRDS(file.path(root, "results", "deg", "GSE74385_program_recurrence.rds"))
v1.up <- unique(v1$up_sym); v1.dn <- unique(v1$dn_sym)
v2.up <- unique(v2$up); v2.dn <- unique(v2$dn)

out <- program
idx <- match(out$symbol, genes)
out$scrna_grade_rho <- rho[out$symbol]
out$scrna_grade_p <- pval[out$symbol]
out$scrna_grade_fdr <- fdr[out$symbol]
out$scrna_mean_count_tumor <- mean.tum[out$symbol]
out$scrna_mean_count_other <- mean.other[out$symbol]
out$scrna_pct_tumor <- pct.tum[out$symbol]
out$scrna_pct_other <- pct.other[out$symbol]
out$scrna_tumor_specificity_log2 <- specificity[out$symbol]
out$aligned_scrna_grade <- ifelse(out$direction == "up_in_WHO_II_vs_I",
                                  out$scrna_grade_rho > 0, out$scrna_grade_rho < 0)
out$recovered_GSE16581 <- ifelse(out$direction == "up_in_WHO_II_vs_I",
                                 out$symbol %in% v1.up, out$symbol %in% v1.dn)
out$recovered_GSE74385 <- ifelse(out$direction == "up_in_WHO_II_vs_I",
                                 out$symbol %in% v2.up, out$symbol %in% v2.dn)

# Internal evidence score (0-7), deliberately separate from druggability/clinical
# evidence that will be joined later from curated external sources.
out$discovery_strength <- pmax(0, 1 - (out$program_rank - 1) / 199)
out$internal_evidence_score <-
  2 * out$discovery_strength +
  ifelse(!is.na(out$aligned_scrna_grade) & out$aligned_scrna_grade, 1, 0) +
  ifelse(!is.na(out$scrna_grade_fdr) & out$scrna_grade_fdr < 0.10, 1, 0) +
  ifelse(!is.na(out$scrna_tumor_specificity_log2) & out$scrna_tumor_specificity_log2 > 0, 1, 0) +
  ifelse(out$recovered_GSE16581, 1, 0) + ifelse(out$recovered_GSE74385, 1, 0)
out <- out[order(-out$internal_evidence_score, out$program_rank), ]

write.csv(out, file.path(outdir, "target_evidence_matrix.csv"), row.names=FALSE, na="")
write.csv(out[seq_len(min(40, nrow(out))), ],
          file.path(outdir, "target_evidence_top40.csv"), row.names=FALSE, na="")

cat("Fixed programme genes:", nrow(program), "\n")
cat("Genes observed in scRNA-seq:", length(genes), "\n")
cat("Tumour cells:", length(tum.cells), " patient pseudobulks:", ncol(pb), "\n")
cat("Aligned grade association (FDR<0.10):",
    sum(out$aligned_scrna_grade & out$scrna_grade_fdr < 0.10, na.rm=TRUE), "\n")
print(out[1:min(20,nrow(out)), c("symbol","direction","internal_evidence_score",
  "scrna_grade_rho","scrna_grade_fdr","scrna_tumor_specificity_log2",
  "recovered_GSE16581","recovered_GSE74385")], row.names=FALSE)
