# GSE74385 endpoint-sensitivity analysis.
# Separates recurrence (R) from malignant progression (M), fits grade- and
# batch-adjusted Firth models, and tests a composite endpoint after excluding
# WHO grade III. This analysis uses the existing frozen score; it does not
# redefine or tune the program in the validation cohort.

suppressPackageStartupMessages({
  library(GEOquery)
  library(logistf)
})

root <- normalizePath(".")
outdir <- file.path(root, "results", "audit_submission")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

m <- read.csv(file.path(root, "results", "deg",
                        "GSE74385_program_recurrence.csv"),
              stringsAsFactors=FALSE)
gse <- getGEO(filename=file.path(root, "data", "raw", "GSE74385",
                                "GSE74385_series_matrix.txt.gz"),
              getGPL=FALSE)
ph <- Biobase::pData(gse)
batch.col <- grep("batch", colnames(ph), ignore.case=TRUE, value=TRUE)[1]
batch <- gsub("batch: ", "", ph[[batch.col]])
m$batch <- factor(batch[match(m$title, ph$title)])
m$grade <- factor(m$grade)
stopifnot(!anyNA(m$batch))

auc_binary <- function(y, score) {
  pos <- score[y == 1]; neg <- score[y == 0]
  mean(outer(pos, neg, ">")) + 0.5*mean(outer(pos, neg, "=="))
}

fit_scenario <- function(label, keep.outcomes, event.outcomes,
                         exclude.grade3=FALSE) {
  d <- m[m$outcome %in% keep.outcomes, , drop=FALSE]
  if (exclude.grade3) d <- d[as.character(d$grade) != "3", , drop=FALSE]
  d$event <- as.integer(d$outcome %in% event.outcomes)
  fit <- logistf(event ~ scale(score) + grade + batch, data=d)
  term <- "scale(score)"
  data.frame(
    scenario=label,
    included_outcomes=paste(keep.outcomes, collapse="/"),
    event_definition=paste(event.outcomes, collapse="/"),
    exclude_grade_III=exclude.grade3,
    n=nrow(d),
    events=sum(d$event),
    nonevents=sum(1-d$event),
    score_mean_event=mean(d$score[d$event==1]),
    score_mean_nonevent=mean(d$score[d$event==0]),
    wilcoxon_p=wilcox.test(score ~ event, data=d, exact=FALSE)$p.value,
    auc=auc_binary(d$event, d$score),
    firth_score_OR=unname(exp(fit$coefficients[term])),
    firth_score_CI_low=unname(exp(fit$ci.lower[term])),
    firth_score_CI_high=unname(exp(fit$ci.upper[term])),
    firth_score_p=unname(fit$prob[term]),
    stringsAsFactors=FALSE
  )
}

results <- rbind(
  fit_scenario("recurrence_only_vs_nonrecurrent",
               c("NR","R"), "R"),
  fit_scenario("malignant_progression_only_vs_nonrecurrent",
               c("NR","M"), "M"),
  fit_scenario("recurrence_or_progression_composite",
               c("NR","R","M"), c("R","M")),
  fit_scenario("composite_excluding_grade_III",
               c("NR","R","M"), c("R","M"), exclude.grade3=TRUE)
)

write.csv(results,
          file.path(outdir, "GSE74385_endpoint_sensitivity.csv"),
          row.names=FALSE)

lines <- c(
  "GSE74385 endpoint sensitivity (all Firth models adjust grade and batch)",
  apply(results, 1, function(x) {
    sprintf(paste0("%s: n=%s, events=%s, AUC=%.3f, Wilcoxon p=%.3g; ",
                   "Firth OR/SD=%.3f (95%% CI %.3f-%.3f), p=%.3g"),
            x[["scenario"]], x[["n"]], x[["events"]],
            as.numeric(x[["auc"]]), as.numeric(x[["wilcoxon_p"]]),
            as.numeric(x[["firth_score_OR"]]),
            as.numeric(x[["firth_score_CI_low"]]),
            as.numeric(x[["firth_score_CI_high"]]),
            as.numeric(x[["firth_score_p"]]))
  }),
  "Interpretation: recurrence and malignant progression are reported separately; the composite is retained only as a sensitivity endpoint."
)
writeLines(lines, file.path(outdir, "GSE74385_endpoint_sensitivity.txt"))
cat(paste(lines, collapse="\n"), "\n")
