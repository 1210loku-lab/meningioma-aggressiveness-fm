# Finalize the target-evidence matrix from the independent raw-count audit.
# This is a deterministic correction step for environments where the historical
# processed Seurat RDS is an unavailable external-drive symlink.
suppressMessages(library(utils))

root <- normalizePath(".")
old.path <- file.path(root, "results", "drug_repurposing", "target_evidence_matrix.csv")
audit.path <- file.path(root, "results", "audit_submission", "independent_GSE206647_target_evidence.csv")
old <- read.csv(old.path, stringsAsFactors=FALSE, check.names=FALSE)
audit <- read.csv(audit.path, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(nrow(old)==200, nrow(audit)==200, !anyDuplicated(old$symbol), !anyDuplicated(audit$symbol))

i <- match(old$symbol, audit$symbol)
stopifnot(!anyNA(i))
old$scrna_grade_rho <- audit$rho[i]
old$scrna_grade_p <- audit$p[i]
old$scrna_grade_fdr <- audit$fdr[i]
old$scrna_tumor_specificity_log2 <- audit$tumour_specificity_log2[i]
old$aligned_scrna_grade <- audit$aligned[i]
old$internal_evidence_score <-
  2 * old$discovery_strength +
  ifelse(!is.na(old$aligned_scrna_grade) & old$aligned_scrna_grade, 1, 0) +
  ifelse(!is.na(old$scrna_grade_fdr) & old$scrna_grade_fdr < 0.10, 1, 0) +
  ifelse(!is.na(old$scrna_tumor_specificity_log2) & old$scrna_tumor_specificity_log2 > 0, 1, 0) +
  ifelse(old$recovered_GSE16581, 1, 0) + ifelse(old$recovered_GSE74385, 1, 0)
old <- old[order(-old$internal_evidence_score, old$program_rank),]

write.csv(old, old.path, row.names=FALSE, na="")
write.csv(head(old,40), file.path(root,"results","drug_repurposing","target_evidence_top40.csv"), row.names=FALSE, na="")
write.csv(old, file.path(root,"docs","Table_S5_target_evidence_matrix.csv"), row.names=FALSE, na="")

aligned <- sum(old$aligned_scrna_grade & old$scrna_grade_fdr<0.10,na.rm=TRUE)
cat("Corrected target matrix: measurable=",sum(!is.na(old$scrna_grade_rho)),
    " aligned_FDR10=",aligned,"\n",sep="")
print(old[old$symbol %in% c("PI3","PITX1","NF2","LTK"),
  c("symbol","scrna_grade_rho","scrna_grade_fdr","internal_evidence_score")],row.names=FALSE)
