# 75_patient_level_metaprogram_correlation.R
# Formal, reproducible computation of PATIENT-LEVEL metaprogram <-> aggressiveness-
# score correlation for every reproducible metaprogram (used by Fig. S5B and Table
# S10). Runs after scripts/71 (which writes the metaprogram signatures/summary and
# the lineage-filtered clean-cell list). Aggregates each metaprogram score and the
# aggressiveness score to the patient level over the NMF-eligible patients
# (>=200 lineage-filtered cells) and computes Spearman correlations, so the
# 0.568/0.589 values in the manuscript regenerate from a numbered script rather
# than a scratchpad file. Also appends BH-FDR of the metaprogram grade p-values.
#
# Inputs : results/scrna/GSE206647_processed.rds
#          results/scrna/cleansed/clean_cells.rds               (from scripts/71)
#          results/scrna/metaprograms_cleansed/metaprogram_signatures.csv (scripts/71)
#          results/scrna/metaprograms_cleansed/metaprogram_summary.csv    (scripts/71)
# Outputs: results/scrna/metaprograms_cleansed/metaprogram_patient_activity.csv
#          results/scrna/metaprograms_cleansed/metaprogram_summary_patientlevel.csv
# Run: Rscript scripts/75_patient_level_metaprogram_correlation.R

suppressMessages({library(Seurat); library(dplyr)})
set.seed(1)
mpdir <- "results/scrna/metaprograms_cleansed"
MIN_CELLS <- 200L                                   # NMF-eligibility threshold (scripts/71)

obj <- readRDS("results/scrna/GSE206647_processed.rds")
clean_cells <- readRDS("results/scrna/cleansed/clean_cells.rds")$clean_cells
md  <- obj@meta.data
sig <- read.csv(file.path(mpdir, "metaprogram_signatures.csv"))
mps <- split(sig$gene, sig$metaprogram)
mps <- mps[paste0("MP", seq_along(mps))]            # preserve MP1..MPn order

sub <- obj[, clean_cells]
sub <- AddModuleScore(sub, features = mps, name = "MP", seed = 1)
sc_cols <- paste0("MP", seq_along(mps))
m <- sub@meta.data
m <- m[m$grade %in% c("I", "II", "III"), ]

# NMF-eligible patients: >= MIN_CELLS lineage-filtered tumour cells
elig <- names(which(table(as.character(m$gsm)) >= MIN_CELLS))
m <- m[m$gsm %in% elig, ]
pa <- m %>% group_by(gsm) %>%
  summarise(across(all_of(c(sc_cols, "AggrScore")), mean),
            grade = grade[1], n_cells = n(), .groups = "drop")
write.csv(pa, file.path(mpdir, "metaprogram_patient_activity.csv"), row.names = FALSE)
cat(sprintf("Patient-level activity: %d NMF-eligible patients (>= %d cells)\n", nrow(pa), MIN_CELLS))

aggr_patient_rho <- sapply(sc_cols, function(s)
  suppressWarnings(cor(pa[[s]], pa$AggrScore, method = "spearman")))

summ <- read.csv(file.path(mpdir, "metaprogram_summary.csv"))
summ$aggr_patient_rho <- round(aggr_patient_rho[summ$metaprogram], 3)
summ$grade_fdr <- round(p.adjust(summ$grade_p, "BH"), 4)
summ <- summ[order(-summ$grade_rho), ]
write.csv(summ, file.path(mpdir, "metaprogram_summary_patientlevel.csv"), row.names = FALSE)

cat("=== metaprogram summary with PATIENT-LEVEL aggr correlation ===\n")
print(summ[, c("metaprogram", "n_patients", "aggr_spearman", "aggr_patient_rho",
               "grade_rho", "grade_p", "grade_fdr")], row.names = FALSE)
cat(sprintf("\nCHECK MP5=%.3f MP6=%.3f (expect 0.568 / 0.589)\n",
            summ$aggr_patient_rho[summ$metaprogram == "MP5"],
            summ$aggr_patient_rho[summ$metaprogram == "MP6"]))
