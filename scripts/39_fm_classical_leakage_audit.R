# Compatibility entry point for the revised patient-level audit.
#
# The historical implementation joined a tokenisation-reordered embedding to
# cells.csv by row position. That code has been retired. Script 52 consumes the
# explicitly audited one-to-one patient mapping created by script 51 and keeps
# patients as the inferential unit throughout.

if (!file.exists("results/scrna/gf_out/tumor_emb_patient_mapped.csv")) {
  source("scripts/51_repair_geneformer_patient_mapping.R")
}
source("scripts/52_geneformer_patient_level_audit.R")
