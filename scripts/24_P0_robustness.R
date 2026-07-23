# 24_P0_robustness.R — P0 审稿防守：GSE74385 batch 混杂重析 + GSE16581 Firth 惩罚多因素
suppressMessages({library(GEOquery)})
if(!requireNamespace("logistf",quietly=TRUE)) install.packages("logistf",repos="https://cloud.r-project.org")
suppressMessages(library(logistf))
set.seed(42)
sink("results/deg/P0_robustness_results.txt", split=TRUE)

cat("=========== P0-A: GSE74385 batch 混杂分析 ===========\n")
m<-read.csv("results/deg/GSE74385_program_recurrence.csv")
g<-getGEO(filename="data/raw/GSE74385/GSE74385_series_matrix.txt.gz",getGPL=FALSE); ph<-pData(g)
batch<-gsub("batch: ","",ph[,grep("batch",colnames(ph))[1]])
m$batch<-factor(batch[match(m$title, ph$title)])
m$recbin<-ifelse(m$recur=="NonRecur",0,ifelse(m$recur=="Recur/Malig",1,NA))
m_known<-m[!is.na(m$recbin),]

cat("\n--- 全体 recurrence × batch（检验总体结论是否抗批次）---\n")
print(table(outcome=m_known$recur, batch=m_known$batch))
cat("\n[总体] logistic recbin ~ score (未校正):\n")
print(summary(glm(recbin~score, m_known, family=binomial))$coef["score",])
cat("\n[总体] logistic recbin ~ score + batch (校正批次):\n")
fit_b<-glm(recbin~score+batch, m_known, family=binomial)
print(summary(fit_b)$coef)
cat("→ 解读：若 score 在校正 batch 后仍显著，总体复发结论抗批次。\n")

cat("\n--- grade-I 子集：batch 与 recurrence 完全共线（致命混杂）---\n")
g1<-m[m$grade==1,]
print(table(outcome=g1$recur, batch=g1$batch))
cat("1R(复发) 全部 batch 2，1NR/1M 全部 batch 1 → grade-I 内无法区分生物学与批次。\n")
cat("grade-I logistic 含 batch 不可估计（完全分离）。结论：grade-I-specific 须降级为 hypothesis-generating。\n")

cat("\n=========== P0-B: GSE16581 Firth 惩罚多因素（EPV=3.3 修正）===========\n")
v<-readRDS("results/deg/GSE16581_program_validation.rds")
d<-data.frame(rec=ifelse(v$rec=="Recur",1,0), score=as.numeric(v$score), grade=factor(v$grade))
d<-d[complete.cases(d),]
cat("n=",nrow(d)," events=",sum(d$rec)," EPV(3 var)=",round(sum(d$rec)/3,2),"\n")

cat("\n--- 标准 logistic（未惩罚对照，存在过拟合风险）---\n")
print(summary(glm(rec~score+grade, d, family=binomial))$coef)

cat("\n--- Firth 惩罚 logistic（小样本/罕见事件稳健）---\n")
ff<-logistf(rec~score+grade, data=d)
res<-data.frame(coef=ff$coefficients, lower95=ff$ci.lower, upper95=ff$ci.upper, p=ff$prob)
print(round(res,4))
cat("→ score 行 p=",signif(res["score","p"],3)," (Firth profile-likelihood)\n")

cat("\n--- bootstrap 1000 次 score 系数稳定性 ---\n")
bs<-replicate(1000,{i<-sample(nrow(d),replace=TRUE)
  cc<-tryCatch(coef(logistf(rec~score+grade, data=d[i,]))["score"],error=function(e)NA); cc})
bs<-bs[!is.na(bs)]
cat("score coef bootstrap: median=",round(median(bs),2)," 2.5%=",round(quantile(bs,.025),2),
    " 97.5%=",round(quantile(bs,.975),2)," P(coef>0)=",round(mean(bs>0),3),"\n")
cat("→ 若 95%CI 不跨 0 且 P(>0)≈1，则'独立于分级'防守在惩罚回归下仍成立。\n")

sink()
cat("\nALL-DONE → results/deg/P0_robustness_results.txt\n")
