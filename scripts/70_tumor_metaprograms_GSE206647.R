# 70_tumor_metaprograms_GSE206647.R
# Discover recurrent malignant cell-state metaprograms (MPs) across patients in
# GSE206647 tumour/tumour-like cells, following the per-patient NMF ->
# cross-patient program clustering paradigm (Gavish/Tirosh 2023). Goal: find a
# RECURRENT aggressive cell state (not patient-private clusters) that (i) tracks
# the fixed bulk aggressiveness program and (ii) expands with WHO grade at the
# patient level. Make-or-break novelty analysis for the resubmission.
#
# Inputs : results/scrna/GSE206647_processed.rds (Seurat; celltype, grade, AggrScore)
# Outputs: results/scrna/metaprograms/
#   - per_patient_programs.csv        (all NMF factors, top genes)
#   - metaprogram_signatures.csv      (recurrent MP gene signatures)
#   - metaprogram_patient_activity.csv(per-patient mean MP activity + grade)
#   - metaprogram_summary.csv         (recurrence, AggrScore corr, grade rho)
#   - cell_metaprogram_scores.rds     (per-cell MP module scores)
# Run from repo root: Rscript scripts/70_tumor_metaprograms_GSE206647.R

suppressMessages({
  library(Seurat); library(RcppML); library(Matrix); library(dplyr)
})
set.seed(1)
outdir <- "results/scrna/metaprograms"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Loading GSE206647 ...")
obj <- readRDS("results/scrna/GSE206647_processed.rds")
tum <- subset(obj, subset = celltype == "Meningioma" & grade %in% c("I","II","III"))
rm(obj); gc()
message(sprintf("Tumour cells: %d across %d patients", ncol(tum), length(unique(tum$gsm))))

patients <- sort(unique(as.character(tum$gsm)))
grade_map <- tapply(as.character(tum$grade), as.character(tum$gsm), function(x) x[1])

## ---- 1. Per-patient NMF -----------------------------------------------------
TOPG <- 50L          # top genes defining each program
KS   <- 4:9          # NMF ranks per patient
prog_list <- list()  # each entry: character vector of top genes
prog_meta <- list()  # patient, k, factor

for (p in patients) {
  cells <- colnames(tum)[tum$gsm == p]
  if (length(cells) < 200) { message(sprintf("  skip %s (%d cells)", p, length(cells))); next }
  sub <- tum[, cells]
  sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  hvg <- VariableFeatures(sub)
  X <- as.matrix(GetAssayData(sub, layer = "data")[hvg, ])   # lognorm, genes x cells
  # relative expression: center each gene, clip negatives (non-negative NMF input)
  X <- X - rowMeans(X); X[X < 0] <- 0
  X <- X[rowSums(X) > 0, , drop = FALSE]
  for (k in KS) {
    fit <- tryCatch(RcppML::nmf(X, k = k, verbose = FALSE, seed = 1),
                    error = function(e) NULL)
    if (is.null(fit)) next
    W <- fit$w; rownames(W) <- rownames(X)
    for (j in seq_len(ncol(W))) {
      top <- names(sort(W[, j], decreasing = TRUE))[1:TOPG]
      id <- sprintf("%s_k%d_f%d", p, k, j)
      prog_list[[id]] <- top
      prog_meta[[id]] <- data.frame(id = id, patient = p, grade = grade_map[p], k = k, factor = j)
    }
  }
  message(sprintf("  %s: %d cells, %d programs", p, length(cells), length(KS) * 1L))
}
prog_meta <- bind_rows(prog_meta)
message(sprintf("Total raw programs: %d", length(prog_list)))

# dump raw programs
prog_df <- do.call(rbind, lapply(names(prog_list), function(id)
  data.frame(id = id, rank = seq_along(prog_list[[id]]), gene = prog_list[[id]])))
prog_df <- merge(prog_meta, prog_df, by = "id")
write.csv(prog_df, file.path(outdir, "per_patient_programs.csv"), row.names = FALSE)

## ---- 2. Cluster programs across patients into metaprograms -------------------
ids <- names(prog_list)
n <- length(ids)
J <- matrix(0, n, n, dimnames = list(ids, ids))
for (i in 1:(n-1)) for (j in (i+1):n) {
  a <- prog_list[[i]]; b <- prog_list[[j]]
  jac <- length(intersect(a, b)) / length(union(a, b))
  J[i, j] <- jac; J[j, i] <- jac
}
d <- as.dist(1 - J)
hc <- hclust(d, method = "average")
# cut so that within-cluster programs share meaningful overlap
cl <- cutree(hc, h = 0.9)   # 1 - Jaccard < 0.9  => Jaccard > 0.1 average linkage
prog_meta$cluster <- cl[prog_meta$id]

# keep metaprograms recurrent across >=3 distinct patients
cl_tab <- prog_meta %>% group_by(cluster) %>%
  summarise(n_prog = n(), n_pat = n_distinct(patient), .groups = "drop") %>%
  filter(n_pat >= 3) %>% arrange(desc(n_pat))
message(sprintf("Recurrent metaprograms (>=3 patients): %d", nrow(cl_tab)))

# signature genes for each recurrent MP: genes appearing in >=25% of member programs
mp_sig <- list()
for (cc in cl_tab$cluster) {
  members <- prog_meta$id[prog_meta$cluster == cc]
  g <- table(unlist(prog_list[members]))
  thr <- max(2, ceiling(0.25 * length(members)))
  sig <- names(sort(g[g >= thr], decreasing = TRUE))
  if (length(sig) < 10) sig <- names(sort(g, decreasing = TRUE))[1:min(30, length(g))]
  mp_sig[[as.character(cc)]] <- head(sig, 50)
}
# name MPs MP1..MPn by patient recurrence
names(mp_sig) <- paste0("MP", seq_along(mp_sig))
mp_cluster_id <- setNames(cl_tab$cluster, names(mp_sig))

sig_df <- do.call(rbind, lapply(names(mp_sig), function(m)
  data.frame(metaprogram = m, cluster = mp_cluster_id[m],
             n_patients = cl_tab$n_pat[cl_tab$cluster == mp_cluster_id[m]],
             rank = seq_along(mp_sig[[m]]), gene = mp_sig[[m]])))
write.csv(sig_df, file.path(outdir, "metaprogram_signatures.csv"), row.names = FALSE)

## ---- 3. Score each MP on all tumour cells -----------------------------------
tum <- AddModuleScore(tum, features = mp_sig, name = "MP", seed = 1)
mp_cols <- paste0("MP", seq_along(mp_sig))
# AddModuleScore appends indices; align names
score_cols <- grep("^MP[0-9]+$", colnames(tum@meta.data), value = TRUE)
score_cols <- paste0("MP", seq_along(mp_sig))
saveRDS(tum@meta.data[, c("gsm","grade","AggrScore", score_cols)],
        file.path(outdir, "cell_metaprogram_scores.rds"))

## ---- 4. Characterise: AggrScore corr + per-patient grade association --------
md <- tum@meta.data
# per-cell correlation of each MP with the fixed bulk AggrScore
aggr_cor <- sapply(score_cols, function(s) cor(md[[s]], md$AggrScore, method = "spearman"))

# per-patient mean MP activity
pat_act <- md %>% group_by(gsm) %>%
  summarise(across(all_of(score_cols), mean), grade = grade[1], .groups = "drop")
pat_act$grade_ord <- as.integer(factor(pat_act$grade, levels = c("I","II","III")))
write.csv(pat_act, file.path(outdir, "metaprogram_patient_activity.csv"), row.names = FALSE)

grade_rho <- sapply(score_cols, function(s)
  suppressWarnings(cor(pat_act[[s]], pat_act$grade_ord, method = "spearman")))
grade_p <- sapply(score_cols, function(s)
  suppressWarnings(cor.test(pat_act[[s]], pat_act$grade_ord, method = "spearman")$p.value))

summ <- data.frame(
  metaprogram = score_cols,
  cluster = mp_cluster_id[score_cols],
  n_patients = cl_tab$n_pat[match(mp_cluster_id[score_cols], cl_tab$cluster)],
  aggr_spearman = round(aggr_cor, 3),
  grade_rho = round(grade_rho, 3),
  grade_p = signif(grade_p, 3)
) %>% arrange(desc(aggr_spearman))
write.csv(summ, file.path(outdir, "metaprogram_summary.csv"), row.names = FALSE)

message("\n==== METAPROGRAM SUMMARY (sorted by AggrScore correlation) ====")
print(summ, row.names = FALSE)
message("\nTop genes of the AggrScore-aligned metaprogram:")
top_mp <- summ$metaprogram[1]
message(top_mp, ": ", paste(head(mp_sig[[top_mp]], 30), collapse = ", "))
message("\nDone. Outputs in ", outdir)
