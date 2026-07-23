# 12_classical_baseline.R — 经典 ML 基线 + 跨队列预测器（防过拟合，诚实基准）
suppressMessages({library(DESeq2); library(glmnet); library(pROC); library(org.Hs.eg.db); library(AnnotationDbi); library(GEOquery)})
set.seed(42); options(timeout=400)
zrow <- function(M){ t(scale(t(M))) }   # gene-wise z-score（跨平台对齐）

## ---- 训练集 GSE136661 (RNA-seq) ----
A <- readRDS("data/raw/GSE136661_assembled.rds"); M<-A$counts; key<-A$key; rownames(key)<-key$gsm; key<-key[colnames(M),]
key$grade <- gsub("WHO ","",key[["pathology:ch1"]])
storage.mode(M)<-"integer"; M<-M[rowSums(M>=5)>=10,]
vsd <- assay(vst(DESeqDataSetFromMatrix(M, key, ~1), blind=TRUE))
sym <- mapIds(org.Hs.eg.db, rownames(vsd), "SYMBOL","ENSEMBL"); keep<-!is.na(sym)&!duplicated(sym)
vsd<-vsd[keep,]; rownames(vsd)<-sym[keep]
y <- factor(ifelse(key$grade=="I","LowGrade","HighGrade"), levels=c("LowGrade","HighGrade"))

## ---- 外验集 GSE16581 (microarray) ----
g16 <- getGEO("GSE16581",GSEMatrix=TRUE,getGPL=TRUE,destdir="data/raw"); e16<-g16[[1]]
fd<-Biobase::fData(e16); X16<-Biobase::exprs(e16); if(max(X16,na.rm=TRUE)>100)X16<-log2(X16+1)
s16<-as.character(fd[[grep("symbol",names(fd),ignore.case=TRUE)[1]]])
ok<-!is.na(s16)&s16!=""; X16<-X16[ok,]; s16<-s16[ok]
X16<-rowsum(X16,s16)/as.vector(table(s16)[rownames(rowsum(X16,s16))])  # collapse 探针->基因均值
pd16<-Biobase::pData(e16); g16grade<-pd16[["who grade:ch1"]]; y16<-factor(ifelse(g16grade=="1","LowGrade","HighGrade"),levels=c("LowGrade","HighGrade"))
rf16<-suppressWarnings(as.integer(pd16[["recurrence_frequency:ch1"]])); rec16<-factor(ifelse(rf16>0,"Recur","NoRecur"))

## ---- 共有基因 + z-score ----
gg <- intersect(rownames(vsd), rownames(X16))
cat("common genes train/val:", length(gg), "\n")
Xtr<-zrow(vsd[gg,]); X16z<-zrow(X16[gg,])

## ---- LASSO 训练 + 重复CV AUC（内部）----
cvfit <- cv.glmnet(t(Xtr), y, family="binomial", alpha=1, nfolds=10, type.measure="auc")
rep_auc <- sapply(1:5, function(r){ set.seed(r); max(cv.glmnet(t(Xtr),y,family="binomial",alpha=1,nfolds=10,type.measure="auc")$cvm) })
cat(sprintf("[INTERNAL] GSE136661 repeated-CV AUC = %.3f (sd %.3f)\n", mean(rep_auc), sd(rep_auc)))

## ---- 外部验证 AUC ----
pr16 <- as.numeric(predict(cvfit, t(X16z), s="lambda.min", type="response"))
roc_g16 <- roc(y16, pr16, quiet=TRUE, levels=c("LowGrade","HighGrade"), direction="<")
cat(sprintf("[EXTERNAL] GSE16581 grade AUC = %.3f (95%%CI %.2f-%.2f)\n", auc(roc_g16), ci.auc(roc_g16)[1], ci.auc(roc_g16)[3]))
roc_r16 <- roc(rec16, pr16, quiet=TRUE, levels=c("NoRecur","Recur"), direction="<")
cat(sprintf("[EXTERNAL] GSE16581 recurrence AUC = %.3f (95%%CI %.2f-%.2f)\n", auc(roc_r16), ci.auc(roc_r16)[1], ci.auc(roc_r16)[3]))

## ---- 对比：侵袭性程序评分作为简单预测器 ----
P<-readRDS("results/deg/GSE16581_program_validation.rds")
score16 <- colMeans(X16[intersect(P$up_sym,rownames(X16)),,drop=FALSE]) - colMeans(X16[intersect(P$dn_sym,rownames(X16)),,drop=FALSE])
cat(sprintf("[PROGRAM-score] GSE16581 grade AUC=%.3f recurrence AUC=%.3f\n",
    auc(roc(y16,score16,quiet=TRUE,direction="<")), auc(roc(rec16,score16,quiet=TRUE,direction="<"))))

nsel <- sum(coef(cvfit,s="lambda.min")!=0)-1
saveRDS(list(cvfit=cvfit, rep_auc=rep_auc, roc_g16=roc_g16, roc_r16=roc_r16, nsel=nsel, genes=gg), "results/deg/classical_baseline.rds")
cat(sprintf("\nLASSO selected %d genes. saved results/deg/classical_baseline.rds\n", nsel))
