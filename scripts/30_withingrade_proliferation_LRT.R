# 30_withingrade_proliferation_LRT.R — R3-A
# 回应独立审稿 R2 两条核心方法学缺陷：
#   #5  程序由 WHO I/II 定义→复发验证 grade 混杂：补 within-grade & grade-residualized 复发分析
#   #10 "beyond proliferation" 需正式增量检验：rec~prolif vs rec~prolif+program 的 LRT + 残差化
# 独立重算（不复用旧 AUC），新旧不一致主动报告。
suppressMessages({library(GEOquery); library(pROC)})
HAS_FIRTH <- requireNamespace("logistf", quietly=TRUE)
set.seed(42); sink("results/deg/R3A_withingrade_proliferation_LRT.txt", split=TRUE)
cat("=== R3-A within-grade recurrence + proliferation incremental LRT ===\n")
cat("Firth(logistf) available:", HAS_FIRTH, "\n\n")

## 程序基因（与 validation 一致）
P <- readRDS("results/deg/GSE16581_program_validation.rds")
prog_up <- P$up_sym; prog_dn <- P$dn_sym
## 增殖签名（与 script 25 一致，offline curated cell-cycle/proliferation）
prolif <- c("MKI67","TOP2A","CCNB1","CCNB2","CCNA2","CDK1","CDC20","BUB1","BUB1B","AURKA",
  "AURKB","PLK1","FOXM1","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CENPA","CENPE","CENPF",
  "KIF11","KIF23","TYMS","RRM2","BIRC5","UBE2C","NUSAP1","ASPM","TPX2","CCNE2","E2F1","CDC45","GINS2")

zscore_sig <- function(expr, up, dn=NULL){
  z <- t(scale(t(expr)))
  u <- intersect(up, rownames(z)); s <- colMeans(z[u,,drop=FALSE], na.rm=TRUE)
  if(!is.null(dn)){ d <- intersect(dn, rownames(z)); s <- s - colMeans(z[d,,drop=FALSE], na.rm=TRUE) }
  s
}
auc_of <- function(score, rec) as.numeric(pROC::auc(rec, score, quiet=TRUE, direction="<"))

## 通用分析：给定 program score(ps), proliferation score(pr), rec(0/1), grade(factor)
analyze <- function(ps, pr, rec, grade, label){
  cat("\n========== ", label, " ==========\n")
  d <- data.frame(rec=rec, ps=as.numeric(scale(ps)), pr=as.numeric(scale(pr)), grade=factor(grade))
  d <- d[complete.cases(d),]
  cat(sprintf("n=%d events=%d ; program-prolif cor r=%.3f\n", nrow(d), sum(d$rec), cor(d$ps,d$pr)))
  cat(sprintf("AUC program=%.3f  proliferation=%.3f\n", auc_of(d$ps,d$rec), auc_of(d$pr,d$rec)))

  ## --- (a) proliferation 增量 LRT：nested glm ---
  m_pr   <- glm(rec~pr,      data=d, family=binomial)
  m_prps <- glm(rec~pr+ps,   data=d, family=binomial)
  lrt <- anova(m_pr, m_prps, test="LRT")
  cat("\n[LRT] rec~prolif  vs  rec~prolif+program :\n")
  cat(sprintf("  added program: deltaDeviance=%.3f  df=%d  p=%.4g\n",
      lrt$Deviance[2], lrt$Df[2], lrt$`Pr(>Chi)`[2]))
  ## 反向：program 之上加 prolif 是否还有增量（对称报告）
  m_ps   <- glm(rec~ps,    data=d, family=binomial)
  m_pspr <- glm(rec~ps+pr, data=d, family=binomial)
  lrt2 <- anova(m_ps, m_pspr, test="LRT")
  cat(sprintf("  reverse (add prolif on program): p=%.4g\n", lrt2$`Pr(>Chi)`[2]))

  ## --- (b) program 对 proliferation 残差化，再测复发 & grade ---
  resid_ps <- residuals(lm(ps~pr, data=d))
  cat("\n[Residualized] program ⟂ proliferation:\n")
  cat(sprintf("  resid-program recurrence AUC=%.3f ; cor(resid, grade-ord)=%.3f\n",
      auc_of(resid_ps,d$rec), suppressWarnings(cor(resid_ps, as.integer(d$grade)))))
  rg <- glm(rec~resid_ps, data=d, family=binomial)
  cat(sprintf("  resid-program ~ recurrence: beta=%.3f p=%.4g\n",
      coef(summary(rg))[2,1], coef(summary(rg))[2,4]))

  ## --- (c) within-grade-I 复发（最临床关键、最易被 grade 混杂）---
  g1 <- d[d$grade==levels(d$grade)[1],]
  cat(sprintf("\n[within grade-I] n=%d events=%d\n", nrow(g1), sum(g1$rec)))
  if(sum(g1$rec)>=3 && nrow(g1)-sum(g1$rec)>=3){
    if(HAS_FIRTH){
      fz <- logistf::logistf(rec~ps, data=g1)
      cat(sprintf("  Firth program-in-gradeI: coef=%.3f CI[%.3f,%.3f] p=%.4g\n",
          fz$coefficients[2], fz$ci.lower[2], fz$ci.upper[2], fz$prob[2]))
    } else {
      gz <- glm(rec~ps, data=g1, family=binomial)
      cat(sprintf("  glm program-in-gradeI: beta=%.3f p=%.4g\n",
          coef(summary(gz))[2,1], coef(summary(gz))[2,4]))
    }
    cat(sprintf("  within-gradeI program recurrence AUC=%.3f\n", auc_of(g1$ps,g1$rec)))
  } else cat("  events 不足，within-grade-I underpowered → 标注不可结论。\n")

  ## --- (d) grade-adjusted Firth：program 在控制 grade 后是否仍有信号（OR per 1-SD）---
  if(HAS_FIRTH){
    f <- logistf::logistf(rec~ps+grade, data=d)
    or <- exp(f$coefficients[2]); orl<-exp(f$ci.lower[2]); oru<-exp(f$ci.upper[2])
    cat(sprintf("\n[grade-adjusted Firth] program OR per 1-SD = %.2f  95%%CI[%.2f, %.2f]  p=%.4g\n",
        or, orl, oru, f$prob[2]))
  }
}

## ---------- GSE16581 (GPL570) ----------
es <- getGEO(filename="data/raw/GSE16581_series_matrix.txt.gz", getGPL=FALSE)
gpl <- getGEO(filename="data/raw/GPL570.soft.gz")
X <- Biobase::exprs(es); pd<-Biobase::pData(es)
if(max(X, na.rm=TRUE) > 100) X <- log2(X + 1)
gtab <- GEOquery::Table(gpl)
fd <- gtab[match(rownames(X), gtab$ID), , drop=FALSE]
symcol <- grep("symbol", names(fd), ignore.case=TRUE, value=TRUE)[1]
sym <- as.character(fd[[symcol]]); keep <- sym!="" & !is.na(sym)
X<-X[keep,]; sym<-sym[keep]; ord<-order(-rowMeans(X)); X<-X[ord,]; sym<-sym[ord]
X<-X[!duplicated(sym),]; rownames(X)<-sym[!duplicated(sym)]
rf <- suppressWarnings(as.integer(pd[["recurrence_frequency:ch1"]])); rec<-ifelse(rf>0,1,0)
gr <- as.integer(pd[["who grade:ch1"]]); if(all(is.na(gr))) gr<-as.integer(gsub("\\D","",pd[["who grade:ch1"]]))
ok <- !is.na(rec) & !is.na(gr)
ps <- zscore_sig(X[,ok], prog_up, prog_dn); pr <- zscore_sig(X[,ok], prolif)
analyze(ps, pr, rec[ok], gr[ok], "GSE16581 (GPL570, n=68)")

## ---------- GSE74385 (Illumina) ----------
suppressMessages(library(illuminaHumanv4.db))
mx<-read.delim(gzfile("data/raw/GSE74385/GSE74385_normalized.txt.gz"),check.names=FALSE,stringsAsFactors=FALSE)
rownames(mx)<-mx[,1]; mx<-mx[,-1]; mx<-mx[,!grepl("Detection",colnames(mx)),drop=FALSE]
sy<-AnnotationDbi::mapIds(illuminaHumanv4.db,keys=rownames(mx),column="SYMBOL",keytype="PROBEID",multiVals="first")
mp<-data.frame(probe=names(sy),sym=as.character(sy)); mp<-mp[mp$sym!="" & !is.na(mp$sym),]
mx$probe<-rownames(mx); mm<-merge(mp,mx,by="probe")
meta<-read.csv("results/deg/GSE74385_program_recurrence.csv")
ev<-as.matrix(sapply(mm[,meta$title],as.numeric)); rownames(ev)<-mm$sym
ev<-ev[order(-rowMeans(ev)),]; ev<-ev[!duplicated(rownames(ev)),]
rec2<-ifelse(meta$recur=="NonRecur",0,ifelse(meta$recur=="Recur/Malig",1,NA))
ps2<-zscore_sig(ev, prog_up, prog_dn); pr2<-zscore_sig(ev, prolif)
ok2<-!is.na(rec2) & !is.na(meta$grade)
analyze(ps2[ok2], pr2[ok2], rec2[ok2], meta$grade[ok2], "GSE74385 (Illumina, known outcomes)")

cat("\n解读：\n - 若 [LRT] 加 program 显著(p<0.05) → program 携带超出 proliferation 的复发信息(可写 'not reducible to proliferation')。\n - 若 within-grade-I events 不足 → 诚实标注 underpowered，不作 grade-I 复发主张。\n - grade-adjusted Firth OR per 1-SD = 手稿应报告的诚实效应量与 CI。\n")
sink(); cat("\nALL-DONE -> results/deg/R3A_withingrade_proliferation_LRT.txt\n")
