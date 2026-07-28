# 13_bulk_figures.R — bulk 臂发表级图（program/GSEA/ROC/survival），不依赖 scRNA
suppressMessages({library(ggplot2); library(patchwork); library(pROC); library(survival); library(GEOquery); library(DESeq2)})
source("scripts/_fig_style.R")
set.seed(42); options(timeout=300); dir.create("results/figures_pub", showWarnings=FALSE)
th <- theme_classic(base_size=12, base_family="Arial")

## A. 侵袭性评分 by WHO grade（GSE136661 训练 + GSE16581 外验）
P136 <- readRDS("results/deg/GSE136661_aggressiveness_program.rds")
d136 <- data.frame(score=P136$score, grade=P136$key$grade, cohort="GSE136661 (RNA-seq)")
V16 <- readRDS("results/deg/GSE16581_program_validation.rds")
d16 <- data.frame(score=V16$score, grade=as.character(V16$grade), cohort="GSE16581 (microarray)")
d16$grade <- c("1"="I","2"="II","3"="III")[d16$grade]
dd <- rbind(d136, d16); dd$grade<-factor(dd$grade,levels=c("I","II","III"))
pA <- ggplot(dd,aes(grade,score,fill=grade))+geom_boxplot(outlier.size=.5)+facet_wrap(~cohort,scales="free_y")+
  th+labs(x="WHO grade",y="Program score",title="Program score across WHO grade")+
  theme(legend.position="none")

## B. Independent recurrence validation (GSE74385)
m74385 <- read.csv("results/deg/GSE74385_program_recurrence.csv")
m74385$outcome <- factor(m74385$outcome, levels=c("NR","M","R"))
m74385_plot <- m74385[!is.na(m74385$outcome) & m74385$outcome %in% c("NR","R"),]
m74385_plot$outcome <- droplevels(m74385_plot$outcome)
p74385 <- wilcox.test(score~outcome, m74385_plot)$p.value
endpoint <- read.csv("results/audit_submission/GSE74385_endpoint_sensitivity.csv")
recurrence_row <- endpoint[endpoint$scenario=="recurrence_only_vs_nonrecurrent",]
pB <- ggplot(m74385_plot, aes(outcome, score, fill=outcome))+
  geom_boxplot(width=0.6, outlier.shape=NA)+
  geom_jitter(width=0.12, size=1.6, alpha=0.7)+
  scale_fill_manual(values=c(NR="#3b6ea5", R="#a14b3d"))+
  th+theme(legend.position="none")+
  labs(x="Recurrence outcome", y="Program score",
       title=sprintf("GSE74385 recurrence association (n=%d)\nFirth OR %.2f (grade+batch), p=%.3f; AUC %.2f",
                     nrow(m74385_plot),recurrence_row$firth_score_OR,
                     recurrence_row$firth_score_p,recurrence_row$auc))

## C. Grade-I recurrence check in GSE74385 (batch-confounded)
gg74385 <- getGEO(filename="data/raw/GSE74385/GSE74385_series_matrix.txt.gz", getGPL=FALSE)
ph74385 <- pData(gg74385)
bt <- gsub("batch: ", "", ph74385[, grep("batch", colnames(ph74385))[1]])
names(bt) <- ph74385$title
g1 <- m74385[m74385$grade==1,]
g1$recur <- factor(g1$recur, levels=c("NonRecur","Recur/Malig"))
g1$batch <- factor(bt[g1$title])
pC <- ggplot(g1, aes(recur, score))+
  geom_boxplot(width=0.55, outlier.shape=NA, fill="grey92")+
  geom_jitter(aes(color=batch), width=0.12, size=2.3, alpha=0.9)+
  scale_color_manual(values=c("1"="#3b6ea5","2"="#d98a00"))+
  th+
  labs(x="WHO grade I only", y="Aggressiveness score", color="Batch",
       title="Grade-I recurrence signal is batch-confounded")

## D. GSEA top pathways
gs <- read.csv("results/deg/GSEA_GOBP_aggressiveness.csv")
gs <- gs[order(gs$NES),]; topu<-tail(gs[gs$NES>0,],8); topd<-head(gs[gs$NES<0,],8)
gg <- rbind(topu,topd); gg$Description<-factor(gg$Description,levels=gg$Description)
pD <- ggplot(gg,aes(NES,Description,fill=NES>0))+geom_col()+th+
  scale_fill_manual(values=c("#3b6ea5","#c0392b"),labels=c("down in high-grade","up in high-grade"))+
  labs(y="",title="GSEA: OXPHOS/ribosome up, synaptic/Wnt down",fill="")+theme(legend.position="bottom")

## E. ROC 曲线（经典 LASSO + program score）
cb <- readRDS("results/deg/classical_baseline.rds")
nested_lines <- readLines("results/deg/R5_nested_lasso_validation.txt")
nested_auc <- as.numeric(sub(".*=", "", grep("^Mean outer-CV AUC=", nested_lines, value=TRUE)))
nested_sd <- as.numeric(sub(".*=", "", grep("^SD across repeats=", nested_lines, value=TRUE)))
roc_df <- function(r,lab) data.frame(spec=rev(r$specificities),sens=rev(r$sensitivities),model=lab)
rc <- rbind(roc_df(cb$roc_g16,sprintf("LASSO grade (AUC %.2f)",auc(cb$roc_g16))),
            roc_df(cb$roc_r16,sprintf("LASSO recurrence (AUC %.2f)",auc(cb$roc_r16))))
pE <- ggplot(rc,aes(1-spec,sens,color=model))+geom_line(linewidth=1)+geom_abline(lty=2,color="grey")+th+
  labs(x="1 - Specificity",y="Sensitivity",title=sprintf("Internal nested CV AUC=%.2f±%.2f\nCurves: external validation",nested_auc,nested_sd),color="")+
  theme(legend.position=c(.6,.2))

## F. 生存 KM（GSE16581 program-score 高 vs 低）
g16<-getGEO(filename="data/raw/GSE16581_series_matrix.txt.gz",getGPL=FALSE); pd<-Biobase::pData(g16)
tts<-suppressWarnings(as.numeric(pd[["tts:ch1"]])); ev<-ifelse(pd[["vital status:ch1"]]=="Deceased",1,ifelse(pd[["vital status:ch1"]]=="Alive",0,NA))
grp<-ifelse(V16$score>median(V16$score),"High","Low")
sd2<-data.frame(tts,ev,grp,score=V16$score)[!is.na(tts)&!is.na(ev),]
fit<-survfit(Surv(tts,ev)~grp,data=sd2); cox_cont<-summary(coxph(Surv(tts,ev)~score,data=sd2))
# 手动构建 KM step dataframe（不依赖 survminer）
km<-do.call(rbind,lapply(seq_along(fit$strata),function(i){
  idx<-(if(i==1)1 else sum(fit$strata[1:(i-1)])+1):sum(fit$strata[1:i])
  data.frame(time=fit$time[idx],surv=fit$surv[idx],grp=sub("grp=","",names(fit$strata)[i]))}))
km<-rbind(data.frame(time=0,surv=1,grp=unique(km$grp)),km)
pF<-ggplot(km,aes(time,surv,color=grp))+geom_step(linewidth=1)+th+ylim(0,1)+
  scale_color_manual(values=c(High="#c0392b",Low="#3b6ea5"))+
  labs(x="Time (days)",y="Overall survival",color="Program",
       title=sprintf("Overall survival (continuous Cox HR=%.2f, p=%.2f)",cox_cont$coefficients[2],cox_cont$coefficients[5]))+
  theme(legend.position=c(.20,.22))

## Forest: cross-cohort grade-adjusted recurrence OR per 1-SD (replaces ROC panel in submission Fig 1)
fst <- read.csv("results/audit_submission/figure_source_data/Fig1_crosscohort_adjusted_OR.csv")
fst$cohort <- factor(fst$cohort, levels=rev(fst$cohort))
fst$lab <- sprintf("OR %.2f (%.2f-%.0f)\nn=%d, %d events; adj. %s",
                   fst$OR, fst$ci_low, fst$ci_high, fst$n, fst$events, fst$adjustment)
pForest <- ggplot(fst, aes(OR, cohort))+
  geom_vline(xintercept=1, linetype=3, color="grey50")+
  geom_errorbar(aes(xmin=ci_low, xmax=ci_high), orientation="y", width=0.18, color="#a14b3d")+
  geom_point(color="#a14b3d", size=3)+
  geom_text(aes(label=lab), hjust=0, nudge_y=0.24, size=3, family="Arial")+
  scale_x_log10(breaks=c(1,3,10,30,100,300), limits=c(0.9,600))+
  th+labs(x="Firth recurrence OR per 1-SD (log scale)", y="",
          title="Grade-adjusted recurrence association")

fig <- (pA|pD)/(pE|pF)+plot_annotation(title="Grade-associated meningioma transcriptomic program - bulk analyses", tag_levels="A") & fig_style
save_fig(fig, "results/figures_pub/Figure_bulk_program", width=15, height=11)
fig_complete <- (pA|pB)/(pD|pForest)+
  plot_annotation(tag_levels="A") & fig_style
save_fig(fig_complete, "results/figures_pub/Figure1_complete", width=15, height=10)
cat("saved results/figures_pub/Figure_bulk_program.{png,pdf} and Figure1_complete.{png,pdf}\n")
