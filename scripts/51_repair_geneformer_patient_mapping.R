# Recover stable cell and patient identifiers for the historical Geneformer
# embedding. The original EmbExtractor output retained grade and AggrScore but
# not cell/gsm, and tokenisation reordered rows. AggrScore is unique per cell in
# this frozen 1,800-cell export, so it provides a deterministic recovery key.
# This script never modifies the historical embedding.

root <- normalizePath(".")
emb.path <- file.path(root, "results", "scrna", "gf_out", "tumor_emb.csv")
cells.path <- file.path(root, "results", "scrna", "gf_export", "cells.csv")
out.path <- file.path(root, "results", "scrna", "gf_out",
                      "tumor_emb_patient_mapped.csv")
audit.dir <- file.path(root, "results", "audit_submission")
dir.create(audit.dir, recursive=TRUE, showWarnings=FALSE)

emb <- read.csv(emb.path, check.names=FALSE, stringsAsFactors=FALSE)
cells <- read.csv(cells.path, stringsAsFactors=FALSE)

stopifnot(nrow(emb) == nrow(cells))
stopifnot(all(c("grade", "AggrScore") %in% names(emb)))
stopifnot(all(c("cell", "grade", "gsm", "n_counts", "AggrScore") %in% names(cells)))
stopifnot(!anyNA(emb$AggrScore), !anyNA(cells$AggrScore))
stopifnot(!anyDuplicated(emb$AggrScore), !anyDuplicated(cells$AggrScore))

# Direct row agreement is diagnostic only; it must not be used for mapping.
direct.grade.agreement <- mean(as.character(emb$grade) == as.character(cells$grade))

mapped.row <- integer(nrow(emb))
score.delta <- numeric(nrow(emb))
for (i in seq_len(nrow(emb))) {
  candidate <- which(as.character(cells$grade) == as.character(emb$grade[i]))
  stopifnot(length(candidate) > 0L)
  delta <- abs(cells$AggrScore[candidate] - emb$AggrScore[i])
  j <- candidate[which.min(delta)]
  mapped.row[i] <- j
  score.delta[i] <- abs(cells$AggrScore[j] - emb$AggrScore[i])
}

tolerance <- 1e-12
stopifnot(max(score.delta) < tolerance)
stopifnot(!anyDuplicated(mapped.row))
stopifnot(length(unique(mapped.row)) == nrow(cells))
stopifnot(all(as.character(emb$grade) == as.character(cells$grade[mapped.row])))

mapped <- emb
mapped$cell_id <- cells$cell[mapped.row]
mapped$gsm <- cells$gsm[mapped.row]
mapped$n_counts_export <- cells$n_counts[mapped.row]
mapped$cells_csv_row <- mapped.row

write.csv(mapped, out.path, row.names=FALSE)

mapping.audit <- data.frame(
  embedding_row=seq_len(nrow(emb)),
  cells_csv_row=mapped.row,
  cell_id=cells$cell[mapped.row],
  gsm=cells$gsm[mapped.row],
  grade=emb$grade,
  AggrScore=emb$AggrScore,
  score_delta=score.delta,
  stringsAsFactors=FALSE
)
write.csv(mapping.audit,
          file.path(audit.dir, "geneformer_patient_mapping_audit.csv"),
          row.names=FALSE)

report <- c(
  "Geneformer historical embedding patient-mapping audit",
  sprintf("Embedding rows: %d", nrow(emb)),
  sprintf("Patients recovered: %d", length(unique(mapped$gsm))),
  sprintf("Direct row-wise grade agreement before mapping: %.3f", direct.grade.agreement),
  sprintf("Unique one-to-one mappings recovered: %d", length(unique(mapped.row))),
  sprintf("Maximum absolute AggrScore matching difference: %.3g", max(score.delta)),
  sprintf("Mapping tolerance: %.1e", tolerance),
  "All mapped grade labels agree: TRUE",
  paste("Mapped embedding:", out.path)
)
writeLines(report, file.path(audit.dir, "geneformer_patient_mapping_audit.txt"))
cat(paste(report, collapse="\n"), "\n")
