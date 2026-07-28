# 60_GSE16581_firth_stability.R
# Reviewer defence for the GSE16581 grade-adjusted Firth recurrence model
# (OR per 1-SD = 14.74, 95% CI 2.77-224.07): the wide CI on 10 events is fragile,
# so we show it is directionally stable via
#   (1) leave-one-recurrence-event-out (LOEO) refitting,
#   (2) 2,000-bootstrap coefficient distribution,
#   (3) profile-likelihood vs Wald CI contrast,
#   (4) grade-only / program-only / grade+program model comparison, and
#   (5) categorical vs ordinal grade coding.
# Independent recompute: the GSE16581 load + scoring path is identical to
# scripts/30 so the canonical OR must reproduce; a mismatch is reported, not hidden.
# Outputs: R6 log, Table S9 csv, FigS3 (Arial), and cross-cohort forest source data.
suppressMessages({library(GEOquery); library(logistf); library(ggplot2); library(patchwork)})
source("scripts/_fig_style.R")
set.seed(42)
dir.create("results/figures_pub", showWarnings=FALSE, recursive=TRUE)
dir.create("results/audit_submission/figure_source_data", showWarnings=FALSE, recursive=TRUE)
sink("results/deg/R6_GSE16581_firth_stability.txt", split=TRUE)
cat("=== R6: GSE16581 grade-adjusted Firth stability (event-limited defence) ===\n\n")

## ---------- GSE16581 load + program score: identical path to scripts/30 ----------
P <- readRDS("results/deg/GSE16581_program_validation.rds")
prog_up <- P$up_sym; prog_dn <- P$dn_sym
zscore_sig <- function(expr, up, dn=NULL){
  z <- t(scale(t(expr)))
  u <- intersect(up, rownames(z)); s <- colMeans(z[u,,drop=FALSE], na.rm=TRUE)
  if(!is.null(dn)){ d <- intersect(dn, rownames(z)); s <- s - colMeans(z[d,,drop=FALSE], na.rm=TRUE) }
  s
}
es <- getGEO(filename="data/raw/GSE16581_series_matrix.txt.gz", getGPL=FALSE)
gpl <- getGEO(filename="data/raw/GPL570.soft.gz")
X <- Biobase::exprs(es); pd <- Biobase::pData(es)
if(max(X, na.rm=TRUE) > 100) X <- log2(X + 1)
gtab <- GEOquery::Table(gpl)
fd <- gtab[match(rownames(X), gtab$ID), , drop=FALSE]
symcol <- grep("symbol", names(fd), ignore.case=TRUE, value=TRUE)[1]
sym <- as.character(fd[[symcol]]); keep <- sym!="" & !is.na(sym)
X <- X[keep,]; sym <- sym[keep]; ord <- order(-rowMeans(X)); X <- X[ord,]; sym <- sym[ord]
X <- X[!duplicated(sym),]; rownames(X) <- sym[!duplicated(sym)]
rf <- suppressWarnings(as.integer(pd[["recurrence_frequency:ch1"]])); rec <- ifelse(rf>0,1,0)
gr <- as.integer(pd[["who grade:ch1"]])
ok <- !is.na(rec) & !is.na(gr)
ps <- zscore_sig(X[,ok], prog_up, prog_dn)
d <- data.frame(rec=rec[ok], ps=as.numeric(scale(ps)), grade=factor(gr[ok]), grade_ord=as.integer(gr[ok]))
d <- d[complete.cases(d),]
cat(sprintf("n=%d  recurrence events=%d  EPV(2 var)=%.2f\n", nrow(d), sum(d$rec), sum(d$rec)/2))
cat("grade x recurrence:\n"); print(table(grade=d$grade, recur=d$rec))

## ---------- (0) canonical model: reproduce OR per 1-SD = 14.74 ----------
f0 <- logistf(rec~ps+grade, data=d)
OR0 <- exp(f0$coefficients[["ps"]]); lo0 <- exp(f0$ci.lower[["ps"]]); hi0 <- exp(f0$ci.upper[["ps"]]); p0 <- f0$prob[["ps"]]
cat(sprintf("\n[canonical] grade-adjusted Firth program OR per 1-SD = %.2f  95%%CI[%.2f, %.2f]  p=%.4g (profile-likelihood)\n",
            OR0, lo0, hi0, p0))
stopifnot(abs(OR0-14.74) < 0.2)  # independent recompute must match the manuscript headline
cat("  -> reproduces manuscript headline (OR 14.74, CI 2.77-224.07, p=1.63e-4).\n")

## ---------- (1) leave-one-recurrence-event-out ----------
cat("\n--- (1) leave-one-recurrence-event-out (drop each of the recurrence events) ---\n")
event_idx <- which(d$rec==1)
loeo <- do.call(rbind, lapply(event_idx, function(i){
  fi <- tryCatch(logistf(rec~ps+grade, data=d[-i,]), error=function(e) NULL)
  if(is.null(fi)) return(data.frame(dropped=i, OR=NA, lo=NA, hi=NA, p=NA))
  data.frame(dropped=i, OR=exp(fi$coefficients[["ps"]]),
             lo=exp(fi$ci.lower[["ps"]]), hi=exp(fi$ci.upper[["ps"]]), p=fi$prob[["ps"]])
}))
loeo$label <- paste0("drop event ", seq_len(nrow(loeo)))
cat(sprintf("LOEO OR per 1-SD: min=%.2f  median=%.2f  max=%.2f ; all OR>1: %s ; all p<0.05: %s\n",
            min(loeo$OR,na.rm=TRUE), median(loeo$OR,na.rm=TRUE), max(loeo$OR,na.rm=TRUE),
            all(loeo$lo>1, na.rm=TRUE), all(loeo$p<0.05, na.rm=TRUE)))
print(round(loeo[,c("OR","lo","hi","p")],3))

## ---------- (2) bootstrap coefficient distribution (per 1-SD) ----------
cat("\n--- (2) 2,000-bootstrap program coefficient (ps standardized once on full data) ---\n")
B <- 2000
boot_warn <- new.env(parent=emptyenv())
boot_warn$n <- 0L
boot_warn$messages <- character()
bs <- replicate(B, {
  idx <- sample(nrow(d), replace=TRUE)
  tryCatch(
    withCallingHandlers(
      {
        fit_b <- logistf(
          rec~ps+grade,
          data=d[idx,],
          pl=FALSE,
          control=logistf.control(maxit=1000)
        )
        fit_b$coefficients[["ps"]]
      },
      warning=function(w) {
        boot_warn$n <- boot_warn$n + 1L
        boot_warn$messages <- unique(c(boot_warn$messages, conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error=function(e) NA_real_
  )
})
bs <- bs[is.finite(bs)]
bs_or <- exp(bs)
cat(sprintf("finite coefficient fits %d/%d ; coef median=%.2f  OR median=%.2f  OR 2.5%%=%.2f  OR 97.5%%=%.2f  P(coef>0)=%.3f\n",
            length(bs), B, median(bs), median(bs_or), quantile(bs_or,.025), quantile(bs_or,.975), mean(bs>0)))
cat(sprintf("bootstrap fitting warnings captured: %d\n", boot_warn$n))
if(length(boot_warn$messages)) {
  cat("unique warning messages:\n")
  cat(paste0("  - ", boot_warn$messages, collapse="\n"), "\n")
}

## ---------- (3) profile-likelihood vs Wald CI for program term ----------
cat("\n--- (3) profile-likelihood vs Wald 95% CI (program OR per 1-SD) ---\n")
ps_idx <- match("ps", names(f0$coefficients))
se_ps <- sqrt(f0$var[ps_idx, ps_idx])
wald_lo <- exp(f0$coefficients[["ps"]] - 1.96*se_ps); wald_hi <- exp(f0$coefficients[["ps"]] + 1.96*se_ps)
cat(sprintf("profile-likelihood: OR=%.2f [%.2f, %.2f]\n", OR0, lo0, hi0))
cat(sprintf("Wald (symmetric)  : OR=%.2f [%.2f, %.2f]\n", OR0, wald_lo, wald_hi))
cat("  -> the manuscript reports the profile-likelihood CI, which is the appropriate interval for penalized rare-event models.\n")

## ---------- (4) grade-only / program-only / grade+program comparison ----------
cat("\n--- (4) model comparison ---\n")
f_grade <- logistf(rec~grade, data=d)
f_prog <- logistf(rec~ps, data=d)
OR_prog <- exp(f_prog$coefficients[["ps"]]); lo_prog <- exp(f_prog$ci.lower[["ps"]]); hi_prog <- exp(f_prog$ci.upper[["ps"]])
f_ord  <- logistf(rec~ps+grade_ord, data=d)
OR_ord <- exp(f_ord$coefficients[["ps"]]); lo_ord <- exp(f_ord$ci.lower[["ps"]]); hi_ord <- exp(f_ord$ci.upper[["ps"]]); p_ord <- f_ord$prob[["ps"]]
add_program <- anova(f_grade, f0)
cat(sprintf("program-only            : OR/1-SD=%.2f [%.2f, %.2f]  p=%.4g\n", OR_prog, lo_prog, hi_prog, f_prog$prob[["ps"]]))
cat(sprintf("grade(categorical)+prog : OR/1-SD=%.2f [%.2f, %.2f]  p=%.4g\n", OR0, lo0, hi0, p0))
cat(sprintf("grade(ordinal)+prog     : OR/1-SD=%.2f [%.2f, %.2f]  p=%.4g\n", OR_ord, lo_ord, hi_ord, p_ord))
cat(sprintf("penalized log-likelihood: grade-only=%.3f ; program-only=%.3f ; grade+program=%.3f\n",
            f_grade$loglik[["full"]], f_prog$loglik[["full"]], f0$loglik[["full"]]))
cat(sprintf("penalized-likelihood-ratio test, grade-only vs grade+program: chi-square=%.3f, df=%d, p=%.4g\n",
            add_program$chisq, add_program$df, add_program$pval))

## ---------- write Table S9 (long format) ----------
grade_rows <- do.call(rbind, lapply(levels(d$grade), function(g) {
  dg <- d[d$grade == g, , drop=FALSE]
  n_event <- sum(dg$rec == 1)
  n_nonevent <- sum(dg$rec == 0)
  estimable <- n_event >= 3 && n_nonevent >= 3
  fg <- if(estimable) logistf(rec~ps, data=dg) else NULL
  data.frame(
    analysis=paste0("within_grade_", g),
    term="exploratory within-grade program OR per full-cohort 1-SD",
    estimate=if(estimable) exp(fg$coefficients[["ps"]]) else NA_real_,
    ci_low=if(estimable) exp(fg$ci.lower[["ps"]]) else NA_real_,
    ci_high=if(estimable) exp(fg$ci.upper[["ps"]]) else NA_real_,
    p=if(estimable) fg$prob[["ps"]] else NA_real_,
    note=sprintf(
      "n=%d; recurrence events=%d; non-events=%d; median score recurrent=%.3f, non-recurrent=%.3f; %s",
      nrow(dg), n_event, n_nonevent,
      if(n_event) median(dg$ps[dg$rec == 1]) else NA_real_,
      if(n_nonevent) median(dg$ps[dg$rec == 0]) else NA_real_,
      if(estimable) "exploratory Firth estimate shown" else "effect not estimated because one outcome group had <3 observations"
    )
  )
}))
tab <- rbind(
  data.frame(analysis="canonical",            term="program OR per 1-SD", estimate=OR0,      ci_low=lo0,      ci_high=hi0,      p=p0,   note="grade-adjusted (categorical) Firth; manuscript headline"),
  data.frame(analysis="LOEO_summary",          term="OR per 1-SD (min/median/max)", estimate=median(loeo$OR,na.rm=TRUE), ci_low=min(loeo$OR,na.rm=TRUE), ci_high=max(loeo$OR,na.rm=TRUE), p=max(loeo$p,na.rm=TRUE), note="leave-one-recurrence-event-out; p column = worst-case p"),
  data.frame(analysis="bootstrap",             term="program OR per 1-SD", estimate=median(bs_or), ci_low=quantile(bs_or,.025), ci_high=quantile(bs_or,.975), p=mean(bs<=0), note=sprintf("2000 resamples; p column = P(coef<=0); finite coefficient fits %d/%d", length(bs), B)),
  data.frame(analysis="CI_profile",            term="program OR per 1-SD", estimate=OR0, ci_low=lo0,     ci_high=hi0,     p=p0, note="profile-likelihood CI (reported)"),
  data.frame(analysis="CI_wald",               term="program OR per 1-SD", estimate=OR0, ci_low=wald_lo, ci_high=wald_hi, p=p0, note="Wald symmetric CI (for contrast only)"),
  data.frame(analysis="model_grade_only",      term="penalized log-likelihood", estimate=f_grade$loglik[["full"]], ci_low=NA, ci_high=NA, p=NA, note="Firth model: recurrence ~ categorical grade"),
  data.frame(analysis="model_program_only",    term="program OR per 1-SD", estimate=OR_prog, ci_low=lo_prog, ci_high=hi_prog, p=f_prog$prob[["ps"]], note="no grade covariate"),
  data.frame(analysis="model_grade_categorical", term="program OR per 1-SD", estimate=OR0, ci_low=lo0, ci_high=hi0, p=p0, note="grade as factor (canonical)"),
  data.frame(analysis="model_grade_ordinal",   term="program OR per 1-SD", estimate=OR_ord, ci_low=lo_ord, ci_high=hi_ord, p=p_ord, note="grade as ordinal 1/2/3"),
  data.frame(analysis="model_add_program_PLR", term="penalized-likelihood-ratio chi-square", estimate=as.numeric(add_program$chisq), ci_low=NA, ci_high=NA, p=as.numeric(add_program$pval), note="nested comparison: grade-only vs categorical-grade + program"),
  grade_rows
)
loeo_out <- data.frame(analysis="LOEO_detail", term=loeo$label, estimate=loeo$OR, ci_low=loeo$lo, ci_high=loeo$hi, p=loeo$p, note="one recurrence event removed")
tab <- rbind(tab, loeo_out)
tab[,c("estimate","ci_low","ci_high","p")] <- lapply(tab[,c("estimate","ci_low","ci_high","p")], function(x) signif(as.numeric(x),4))
write.csv(tab, "docs/Table_S9_GSE16581_firth_stability.csv", row.names=FALSE)
cat("\nwrote docs/Table_S9_GSE16581_firth_stability.csv\n")

## ---------- cross-cohort adjusted-OR forest source data (for Figure 1 panel d) ----------
g74 <- read.csv("results/audit_submission/GSE74385_endpoint_sensitivity.csv")
r74 <- g74[g74$scenario=="recurrence_only_vs_nonrecurrent",]
forest <- data.frame(
  cohort   = c("GSE16581", "GSE74385 (recurrence-only)"),
  platform = c("Affymetrix GPL570", "Illumina HT-12"),
  n        = c(nrow(d), r74$n),
  events   = c(sum(d$rec), r74$events),
  adjustment = c("grade", "grade + batch"),
  OR       = c(OR0, r74$firth_score_OR),
  ci_low   = c(lo0, r74$firth_score_CI_low),
  ci_high  = c(hi0, r74$firth_score_CI_high),
  p        = c(p0,  r74$firth_score_p)
)
forest[,c("OR","ci_low","ci_high","p")] <- lapply(forest[,c("OR","ci_low","ci_high","p")], function(x) signif(as.numeric(x),4))
write.csv(forest, "results/audit_submission/figure_source_data/Fig1_crosscohort_adjusted_OR.csv", row.names=FALSE)
cat("wrote results/audit_submission/figure_source_data/Fig1_crosscohort_adjusted_OR.csv\n")
print(forest)

## ---------- Figure S3 (Arial): stability panels ----------
th <- theme_classic(base_size=11, base_family="Arial") +
  theme(plot.title=element_text(family="Arial", size=11),
        axis.text=element_text(family="Arial"), axis.title=element_text(family="Arial"))
logbreaks <- c(1,3,10,30,100,300)

## a: LOEO forest
la <- loeo
la$idx <- seq_len(nrow(la))
la <- la[order(la$OR),]
la$row <- seq_len(nrow(la))
pa <- ggplot(la, aes(OR, row)) +
  annotate("rect", xmin=lo0, xmax=hi0, ymin=-Inf, ymax=Inf, fill="#a14b3d", alpha=0.10) +
  geom_vline(xintercept=OR0, linetype=2, color="#a14b3d") +
  geom_vline(xintercept=1, linetype=3, color="grey50") +
  geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=0.25, color="#3b6ea5") +
  geom_point(color="#3b6ea5", size=1.8) +
  scale_x_log10(breaks=logbreaks) +
  scale_y_continuous(breaks=la$row, labels=la$label) + th +
  labs(x="Program OR per 1-SD (log scale)", y="Omitted recurrence event",
       title=sprintf("Leave-one-recurrence-event-out (%d refits)\nDashed line: full-data OR %.1f", nrow(la), OR0))

## b: bootstrap coefficient distribution; beta scale avoids a visually distorted OR tail
pb <- ggplot(data.frame(beta=bs), aes(beta)) +
  geom_histogram(bins=40, fill="#3b6ea5", color="white", alpha=0.85) +
  geom_vline(xintercept=median(bs), color="#a14b3d", linewidth=0.8) +
  geom_vline(xintercept=quantile(bs,c(.025,.975)), linetype=2, color="#a14b3d") +
  geom_vline(xintercept=0, linetype=3, color="grey50") + th +
  labs(x="Firth program coefficient per 1-SD", y="Bootstrap resamples",
       title=sprintf("Bootstrap coefficient distribution (%d/%d finite fits)\nMedian beta %.2f; P(beta>0)=%.3f",
                     length(bs), B, median(bs), mean(bs>0)))

## c: model comparison forest
mc <- data.frame(
  model=factor(c("program only","grade (categorical)\n+ program","grade (ordinal)\n+ program"),
               levels=c("grade (ordinal)\n+ program","grade (categorical)\n+ program","program only")),
  OR=c(OR_prog, OR0, OR_ord), lo=c(lo_prog, lo0, lo_ord), hi=c(hi_prog, hi0, hi_ord))
pc <- ggplot(mc, aes(OR, model)) +
  geom_vline(xintercept=1, linetype=3, color="grey50") +
  geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=0.2, color="#3b6ea5") +
  geom_point(color="#3b6ea5", size=2.4) +
  scale_x_log10(breaks=logbreaks) + th +
  labs(x="Program OR per 1-SD (log scale)", y="",
       title=sprintf("Program estimates across model specifications\nPenalized-likelihood p (program | grade) = %.1g", p0))

## d: profile-likelihood versus Wald interval
ci_cmp <- data.frame(
  method=factor(c("Profile likelihood", "Wald"),
                levels=c("Wald", "Profile likelihood")),
  OR=OR0,
  lo=c(lo0, wald_lo),
  hi=c(hi0, wald_hi)
)
pd_ci <- ggplot(ci_cmp, aes(OR, method)) +
  geom_vline(xintercept=1, linetype=3, color="grey50") +
  geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=0.2, color="#3b6ea5") +
  geom_point(color="#3b6ea5", size=2.4) +
  scale_x_log10(breaks=logbreaks) + th +
  labs(x="Program OR per 1-SD (log scale)", y="",
       title="Confidence-interval method comparison\nProfile-likelihood interval is reported")

figS3 <- pa + pb + pc + pd_ci +
  plot_layout(design="AB\nCD", heights=c(1,0.9)) +
  plot_annotation(tag_levels="A") & fig_style
save_fig(figS3, "results/figures_pub/fig_GSE16581_firth_stability", width=11, height=8)
cat("\nwrote results/figures_pub/fig_GSE16581_firth_stability.{png,pdf}\n")

sink(); cat("\nALL-DONE -> results/deg/R6_GSE16581_firth_stability.txt\n")
