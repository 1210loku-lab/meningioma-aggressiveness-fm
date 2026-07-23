# Independent nested-CV audit of the classical grade classifier.
# Keeps the original analysis intact and writes a separate robustness result.
suppressMessages({
  library(DESeq2)
  library(glmnet)
  library(pROC)
})

set.seed(20260713)

a <- readRDS("data/raw/GSE136661_assembled.rds")
counts <- a$counts
key <- a$key
rownames(key) <- key$gsm
key <- key[colnames(counts), , drop = FALSE]
grade <- gsub("WHO ", "", key[["pathology:ch1"]])
keep_samples <- grade %in% c("I", "II", "III")
counts <- counts[, keep_samples, drop = FALSE]
grade <- grade[keep_samples]
y <- as.integer(grade != "I")

storage.mode(counts) <- "integer"
counts <- counts[rowSums(counts >= 5) >= 10, , drop = FALSE]
dds <- DESeqDataSetFromMatrix(counts, data.frame(y = y), design = ~1)
x <- t(assay(vst(dds, blind = TRUE)))

z_train_apply <- function(x_train, x_test) {
  mu <- colMeans(x_train)
  sdv <- apply(x_train, 2, sd)
  valid <- is.finite(sdv) & sdv > 0
  list(
    train = sweep(sweep(x_train[, valid, drop = FALSE], 2, mu[valid], "-"), 2, sdv[valid], "/"),
    test = sweep(sweep(x_test[, valid, drop = FALSE], 2, mu[valid], "-"), 2, sdv[valid], "/")
  )
}

n_repeats <- 10
n_outer <- 5
repeat_auc <- numeric(n_repeats)
all_predictions <- vector("list", n_repeats)

for (r in seq_len(n_repeats)) {
  set.seed(20260713 + r)
  fold_id <- integer(length(y))
  for (cls in 0:1) {
    idx <- sample(which(y == cls))
    fold_id[idx] <- rep(seq_len(n_outer), length.out = length(idx))
  }
  pred <- rep(NA_real_, length(y))
  for (fold in seq_len(n_outer)) {
    te <- which(fold_id == fold)
    tr <- which(fold_id != fold)
    z <- z_train_apply(x[tr, , drop = FALSE], x[te, , drop = FALSE])
    set.seed(20260713 + r * 100 + fold)
    inner <- cv.glmnet(
      z$train, y[tr], family = "binomial", alpha = 1,
      nfolds = 5, type.measure = "auc", standardize = FALSE
    )
    pred[te] <- as.numeric(predict(inner, z$test, s = "lambda.min", type = "response"))
  }
  repeat_auc[r] <- as.numeric(auc(roc(y, pred, levels = c(0, 1), direction = "<", quiet = TRUE)))
  all_predictions[[r]] <- data.frame(repeat_id = r, sample = colnames(counts), grade = grade, y = y, probability = pred)
}

out <- c(
  "=== Independent nested-CV audit: GSE136661 WHO I vs II/III LASSO ===",
  sprintf("n=%d; WHO I=%d; WHO II/III=%d", length(y), sum(y == 0), sum(y == 1)),
  "Design: 10 repeats of stratified 5-fold outer CV; lambda selected only within each outer-training set by 5-fold inner CV; z-scaling fitted on each outer-training set.",
  sprintf("Mean outer-CV AUC=%.3f", mean(repeat_auc)),
  sprintf("SD across repeats=%.3f", sd(repeat_auc)),
  sprintf("Range across repeats=%.3f-%.3f", min(repeat_auc), max(repeat_auc)),
  paste0("Per-repeat AUC: ", paste(sprintf("%.3f", repeat_auc), collapse = ", ")),
  "Interpretation: use this estimate for internal discrimination; the original cv.glmnet curve estimate is retained only as a historical analysis output."
)

writeLines(out, "results/deg/R5_nested_lasso_validation.txt")
write.csv(do.call(rbind, all_predictions), "results/deg/R5_nested_lasso_predictions.csv", row.names = FALSE)
cat(paste(out, collapse = "\n"), "\n")
