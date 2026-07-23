# Compatibility entry point for the revised patient-level Geneformer figure.
# The former composite consumed cell-level permutation and row-misaligned
# benchmark outputs; Figure 5 is now rebuilt from patient-level audit tables.

required <- c(
  "results/audit_submission/geneformer_patient_level_axis_scores.csv",
  "results/audit_submission/geneformer_patient_level_axis_summary.csv",
  "results/audit_submission/geneformer_patient_level_lopo_per_repeat.csv"
)
if (!all(file.exists(required))) {
  source("scripts/52_geneformer_patient_level_audit.R")
}
source("scripts/56_rebuild_review_response_figures.R")
