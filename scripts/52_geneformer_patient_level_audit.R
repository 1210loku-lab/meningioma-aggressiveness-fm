# Patient-level Geneformer audit.
#
# Primary axis analysis: average the 768-dimensional embedding within each
# patient, fit PCA to the 16 patient means, and assess the best absolute
# Spearman correlation among the top 10 PCs using patient-label permutations.
#
# Sensitivity axis analysis: fit the unsupervised PCA in cell space, average PC
# scores within patient, and run the same patient-level test. This is retained
# only to show how the result depends on the order of aggregation and PCA.
#
# Classification: compare patient-mean Geneformer embeddings with patient-mean
# program-gene expression using repeated nested leave-one-patient-out ridge
# logistic regression. No cell-weighted AUC is used.

suppressPackageStartupMessages({
  library(Matrix)
  library(glmnet)
})

set.seed(20260717)
root <- normalizePath(".")
outdir <- file.path(root, "results", "audit_submission")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

emb <- read.csv(file.path(root, "results", "scrna", "gf_out",
                          "tumor_emb_patient_mapped.csv"),
                check.names=FALSE, stringsAsFactors=FALSE)
emb.cols <- grep("^[0-9]+$", names(emb), value=TRUE)
stopifnot(length(emb.cols) == 768L)
stopifnot(all(c("grade", "AggrScore", "gsm", "cell_id") %in% names(emb)))
stopifnot(nrow(emb) == 1800L, length(unique(emb$gsm)) == 16L)

grade.map <- c(I=1, II=2, III=3)
patient.meta <- unique(emb[,c("gsm", "grade")])
stopifnot(!anyDuplicated(patient.meta$gsm))
patient.meta$grade_ord <- unname(grade.map[patient.meta$grade])
patient.meta$binary_grade <- as.integer(patient.meta$grade != "I")
stopifnot(table(patient.meta$grade)["I"] == 8L)
stopifnot(sum(patient.meta$grade %in% c("II", "III")) == 8L)

mean_by_patient <- function(x, gsm, patient.order) {
  z <- rowsum(as.matrix(x), group=gsm, reorder=FALSE)
  n <- as.numeric(table(factor(gsm, levels=rownames(z))))
  z <- z / n
  z[match(patient.order, rownames(z)),,drop=FALSE]
}

patient.order <- patient.meta$gsm
x.cell <- as.matrix(emb[,emb.cols,drop=FALSE])
x.patient <- mean_by_patient(x.cell, emb$gsm, patient.order)
aggr.patient <- as.numeric(mean_by_patient(matrix(emb$AggrScore,ncol=1),
                                           emb$gsm, patient.order))
y.ord <- patient.meta$grade_ord
y.bin <- patient.meta$binary_grade

best_patient_permutation <- function(scores, y, nperm=100000L, seed=1L) {
  scores <- as.matrix(scores)
  ranked <- apply(scores, 2, rank, ties.method="average")
  ranked <- sweep(ranked, 2, colMeans(ranked), "-")
  ranked <- sweep(ranked, 2, sqrt(colSums(ranked^2)), "/")
  yr <- rank(y, ties.method="average")
  yr <- yr - mean(yr)
  ynorm <- sqrt(sum(yr^2))
  obs.rho <- as.numeric(crossprod(yr / ynorm, ranked))
  best <- which.max(abs(obs.rho))

  set.seed(seed)
  null <- numeric(nperm)
  batch.size <- 5000L
  starts <- seq.int(1L, nperm, by=batch.size)
  for (start in starts) {
    b <- min(batch.size, nperm - start + 1L)
    yp <- t(replicate(b, sample(yr, replace=FALSE))) / ynorm
    cor.mat <- abs(yp %*% ranked)
    null[start:(start+b-1L)] <- apply(cor.mat, 1, max)
  }
  list(
    best_pc=best,
    rho=obs.rho[best],
    all_rho=obs.rho,
    null_mean=mean(null),
    null_q95=unname(quantile(null, 0.95)),
    null_q99=unname(quantile(null, 0.99)),
    empirical_p=(sum(null >= abs(obs.rho[best])) + 1) / (nperm + 1),
    nperm=nperm
  )
}

# Primary: patients are the independent observations before PCA.
patient.pca <- prcomp(x.patient, center=TRUE, scale.=FALSE, rank.=10)
patient.pc <- patient.pca$x[,seq_len(10),drop=FALSE]
primary.axis <- best_patient_permutation(patient.pc, y.ord,
                                         nperm=100000L, seed=20260717L)
primary.aggr.rho <- cor(patient.pc[,primary.axis$best_pc], aggr.patient,
                        method="spearman")

# Sensitivity: unsupervised cell-space PCA followed by patient aggregation.
cell.pca <- prcomp(x.cell, center=TRUE, scale.=FALSE, rank.=10)
cell.pc.patient <- mean_by_patient(cell.pca$x[,seq_len(10),drop=FALSE],
                                   emb$gsm, patient.order)
sensitivity.axis <- best_patient_permutation(cell.pc.patient, y.ord,
                                             nperm=100000L, seed=20260718L)
sensitivity.aggr.rho <- cor(cell.pc.patient[,sensitivity.axis$best_pc],
                            aggr.patient, method="spearman")

axis.summary <- data.frame(
  analysis=c("patient_mean_embedding_then_PCA_primary",
             "cell_PCA_then_patient_mean_sensitivity"),
  n_patients=16L,
  best_pc=c(primary.axis$best_pc, sensitivity.axis$best_pc),
  spearman_grade=c(primary.axis$rho, sensitivity.axis$rho),
  spearman_program_score=c(primary.aggr.rho, sensitivity.aggr.rho),
  permutation_unit="patient",
  n_permutations=c(primary.axis$nperm, sensitivity.axis$nperm),
  null_mean_abs_rho_max=c(primary.axis$null_mean, sensitivity.axis$null_mean),
  null_q95_abs_rho_max=c(primary.axis$null_q95, sensitivity.axis$null_q95),
  null_q99_abs_rho_max=c(primary.axis$null_q99, sensitivity.axis$null_q99),
  empirical_p=c(primary.axis$empirical_p, sensitivity.axis$empirical_p)
)
write.csv(axis.summary,
          file.path(outdir, "geneformer_patient_level_axis_summary.csv"),
          row.names=FALSE)

patient.axis.data <- data.frame(
  gsm=patient.order,
  grade=patient.meta$grade,
  grade_ord=y.ord,
  binary_grade=y.bin,
  n_cells=as.integer(table(factor(emb$gsm, levels=patient.order))),
  mean_AggrScore=aggr.patient,
  primary_axis=patient.pc[,primary.axis$best_pc],
  sensitivity_axis=cell.pc.patient[,sensitivity.axis$best_pc]
)
write.csv(patient.axis.data,
          file.path(outdir, "geneformer_patient_level_axis_scores.csv"),
          row.names=FALSE)

# Build patient-mean classical program-gene expression from the same cells.
cells <- read.csv(file.path(root, "results", "scrna", "gf_export", "cells.csv"),
                  stringsAsFactors=FALSE)
counts <- t(readMM(file.path(root, "results", "scrna", "gf_export", "counts.mtx")))
genes <- read.csv(file.path(root, "results", "scrna", "gf_export", "genes.csv"),
                  stringsAsFactors=FALSE)
program <- read.csv(file.path(root, "docs", "Table_S2_aggressiveness_program_genes.csv"),
                    stringsAsFactors=FALSE)
stopifnot(nrow(cells) == nrow(counts))
program.ids <- unique(c(program$symbol, program$ensembl))
keep <- as.character(genes$symbol) %in% program.ids
stopifnot(sum(keep) == 159L)

libsize <- Matrix::rowSums(counts)
libsize[libsize == 0] <- 1
x.classical.cell <- Diagonal(x=1e4/libsize) %*% counts[,keep,drop=FALSE]
x.classical.cell@x <- log1p(x.classical.cell@x)
x.classical.patient <- mean_by_patient(x.classical.cell, cells$gsm, patient.order)

make_stratified_foldid <- function(y, k, seed) {
  set.seed(seed)
  foldid <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- sample(which(y == cls))
    foldid[idx] <- rep(seq_len(k), length.out=length(idx))
  }
  foldid
}

nested_lopo <- function(x, y, feature_set, repeats=10L) {
  x <- as.matrix(x)
  all.pred <- vector("list", repeats)
  summary <- vector("list", repeats)
  for (r in seq_len(repeats)) {
    pred <- numeric(length(y))
    selected.lambda <- numeric(length(y))
    for (i in seq_along(y)) {
      tr <- setdiff(seq_along(y), i)
      te <- i
      foldid <- make_stratified_foldid(y[tr], 5L,
                                       seed=20260717L + 1000L*r + i)
      # glmnet warns whenever a training class has <8 observations. That is an
      # expected property of this 16-patient audit, not a numerical failure;
      # suppress the repeated console warning and disclose it in the report.
      cvfit <- suppressWarnings(cv.glmnet(
        x[tr,,drop=FALSE], y[tr], family="binomial",
        alpha=0, foldid=foldid, type.measure="deviance",
        standardize=TRUE, nlambda=100))
      pred[te] <- as.numeric(predict(cvfit, x[te,,drop=FALSE],
                                     s="lambda.1se", type="response"))
      selected.lambda[te] <- cvfit$lambda.1se
    }
    all.pred[[r]] <- data.frame(
      feature_set=feature_set,
      repeat_id=r,
      gsm=patient.order,
      grade=patient.meta$grade,
      y=y,
      probability=pred,
      selected_lambda=selected.lambda
    )
    pos <- pred[y == 1]; neg <- pred[y == 0]
    auc <- mean(outer(pos, neg, ">")) + 0.5*mean(outer(pos, neg, "=="))
    summary[[r]] <- data.frame(feature_set=feature_set, repeat_id=r,
                               n_patients=length(y), lopo_auc=auc)
  }
  list(predictions=do.call(rbind,all.pred), summary=do.call(rbind,summary))
}

fm.cv <- nested_lopo(x.patient, y.bin, "Geneformer_patient_mean_embedding")
classical.cv <- nested_lopo(x.classical.patient, y.bin,
                            "Program_gene_patient_mean_expression")
cv.pred <- rbind(fm.cv$predictions, classical.cv$predictions)
cv.repeat <- rbind(fm.cv$summary, classical.cv$summary)

# Median prediction across repeated inner-fold assignments is used for the
# primary patient-level AUC and paired patient bootstrap interval.
median.pred <- aggregate(probability ~ feature_set + gsm + grade + y,
                         cv.pred, median)
auc_binary <- function(y, p) {
  pos <- p[y == 1]; neg <- p[y == 0]
  mean(outer(pos, neg, ">")) + 0.5*mean(outer(pos, neg, "=="))
}

features <- unique(median.pred$feature_set)
primary.auc <- setNames(numeric(length(features)), features)
for (f in features) {
  d <- median.pred[median.pred$feature_set == f,]
  primary.auc[f] <- auc_binary(d$y, d$probability)
}

pred.wide <- reshape(median.pred[,c("gsm","grade","y","feature_set","probability")],
                     idvar=c("gsm","grade","y"), timevar="feature_set",
                     direction="wide")
fm.col <- paste0("probability.", "Geneformer_patient_mean_embedding")
cl.col <- paste0("probability.", "Program_gene_patient_mean_expression")
set.seed(20260719)
B <- 10000L
boot <- matrix(NA_real_, nrow=B, ncol=3,
               dimnames=list(NULL,c("fm_auc","classical_auc","delta_classical_minus_fm")))
idx0 <- which(pred.wide$y == 0); idx1 <- which(pred.wide$y == 1)
for (b in seq_len(B)) {
  idx <- c(sample(idx0, length(idx0), replace=TRUE),
           sample(idx1, length(idx1), replace=TRUE))
  yy <- pred.wide$y[idx]
  af <- auc_binary(yy, pred.wide[[fm.col]][idx])
  ac <- auc_binary(yy, pred.wide[[cl.col]][idx])
  boot[b,] <- c(af, ac, ac-af)
}

cv.summary <- data.frame(
  feature_set=features,
  n_patients=16L,
  repeated_nested_LOPO_auc_mean=sapply(features, function(f)
    mean(cv.repeat$lopo_auc[cv.repeat$feature_set == f])),
  repeated_nested_LOPO_auc_sd=sapply(features, function(f)
    sd(cv.repeat$lopo_auc[cv.repeat$feature_set == f])),
  repeated_nested_LOPO_auc_min=sapply(features, function(f)
    min(cv.repeat$lopo_auc[cv.repeat$feature_set == f])),
  repeated_nested_LOPO_auc_max=sapply(features, function(f)
    max(cv.repeat$lopo_auc[cv.repeat$feature_set == f])),
  median_prediction_auc=unname(primary.auc[features]),
  bootstrap_ci_low=c(quantile(boot[,"fm_auc"],0.025),
                     quantile(boot[,"classical_auc"],0.025))[
                       match(features,c("Geneformer_patient_mean_embedding",
                                        "Program_gene_patient_mean_expression"))],
  bootstrap_ci_high=c(quantile(boot[,"fm_auc"],0.975),
                      quantile(boot[,"classical_auc"],0.975))[
                        match(features,c("Geneformer_patient_mean_embedding",
                                         "Program_gene_patient_mean_expression"))]
)

delta.summary <- data.frame(
  contrast="Program_gene_expression_minus_Geneformer_embedding",
  observed_auc_difference=unname(primary.auc["Program_gene_patient_mean_expression"] -
                                 primary.auc["Geneformer_patient_mean_embedding"]),
  paired_patient_bootstrap_ci_low=unname(quantile(boot[,3],0.025)),
  paired_patient_bootstrap_ci_high=unname(quantile(boot[,3],0.975)),
  n_bootstrap=B
)

write.csv(cv.pred,
          file.path(outdir, "geneformer_patient_level_lopo_predictions.csv"),
          row.names=FALSE)
write.csv(cv.repeat,
          file.path(outdir, "geneformer_patient_level_lopo_per_repeat.csv"),
          row.names=FALSE)
write.csv(cv.summary,
          file.path(outdir, "geneformer_patient_level_lopo_summary.csv"),
          row.names=FALSE)
write.csv(delta.summary,
          file.path(outdir, "geneformer_patient_level_lopo_auc_difference.csv"),
          row.names=FALSE)

report <- c(
  "Patient-level Geneformer audit",
  sprintf("Patients: %d (WHO I/II/III = %s)", nrow(patient.meta),
          paste(as.integer(table(factor(patient.meta$grade,
                                        levels=c("I","II","III")))),collapse="/")),
  sprintf(paste0("Primary patient-mean-embedding PCA: PC%d rho(grade)=%.3f, ",
                 "rho(program)=%.3f, patient-permutation p=%.4g ",
                 "(null mean %.3f; 95th %.3f; %d permutations)"),
          primary.axis$best_pc, primary.axis$rho, primary.aggr.rho,
          primary.axis$empirical_p, primary.axis$null_mean,
          primary.axis$null_q95, primary.axis$nperm),
  sprintf(paste0("Sensitivity cell-PCA-then-patient-mean: PC%d rho(grade)=%.3f, ",
                 "rho(program)=%.3f, patient-permutation p=%.4g ",
                 "(null mean %.3f; 95th %.3f; %d permutations)"),
          sensitivity.axis$best_pc, sensitivity.axis$rho,
          sensitivity.aggr.rho, sensitivity.axis$empirical_p,
          sensitivity.axis$null_mean, sensitivity.axis$null_q95,
          sensitivity.axis$nperm),
  sprintf("Geneformer repeated nested LOPO mean AUC=%.3f; median-prediction AUC=%.3f (patient bootstrap 95%% CI %.3f-%.3f)",
          cv.summary$repeated_nested_LOPO_auc_mean[cv.summary$feature_set=="Geneformer_patient_mean_embedding"],
          cv.summary$median_prediction_auc[cv.summary$feature_set=="Geneformer_patient_mean_embedding"],
          cv.summary$bootstrap_ci_low[cv.summary$feature_set=="Geneformer_patient_mean_embedding"],
          cv.summary$bootstrap_ci_high[cv.summary$feature_set=="Geneformer_patient_mean_embedding"]),
  sprintf("Classical repeated nested LOPO mean AUC=%.3f; median-prediction AUC=%.3f (patient bootstrap 95%% CI %.3f-%.3f)",
          cv.summary$repeated_nested_LOPO_auc_mean[cv.summary$feature_set=="Program_gene_patient_mean_expression"],
          cv.summary$median_prediction_auc[cv.summary$feature_set=="Program_gene_patient_mean_expression"],
          cv.summary$bootstrap_ci_low[cv.summary$feature_set=="Program_gene_patient_mean_expression"],
          cv.summary$bootstrap_ci_high[cv.summary$feature_set=="Program_gene_patient_mean_expression"]),
  sprintf("AUC difference classical-minus-Geneformer=%.3f (paired patient bootstrap 95%% CI %.3f-%.3f)",
          delta.summary$observed_auc_difference,
          delta.summary$paired_patient_bootstrap_ci_low,
          delta.summary$paired_patient_bootstrap_ci_high),
  "Interpretation: all inferential units are patients; cell-weighted AUC is not reported.",
  "Small-sample warning: every outer training set has fewer than 8 patients in one class; AUCs and conditional bootstrap intervals are descriptive and cannot establish comparative superiority."
)
writeLines(report, file.path(outdir, "geneformer_patient_level_audit.txt"))
cat(paste(report, collapse="\n"), "\n")
