suppressMessages({library(GEOquery);library(org.Hs.eg.db);library(AnnotationDbi);library(survival)}); options(timeout=400)
P <- readRDS("results/deg/GSE136661_aggressiveness_program.rds")
up_sym <- na.omit(mapIds(org.Hs.eg.db, P$up, "SYMBOL","ENSEMBL"))
dn_sym <- na.omit(mapIds(org.Hs.eg.db, P$dn, "SYMBOL","ENSEMBL"))
cat("program mapped to symbols: up=",length(up_sym)," dn=",length(dn_sym),"\n")
es <- getGEO(filename="data/raw/GSE16581_series_matrix.txt.gz", getGPL=FALSE)
gpl <- getGEO(filename="data/raw/GPL570.soft.gz")
X <- Biobase::exprs(es); pd<-Biobase::pData(es)
if(max(X,na.rm=TRUE)>100) X<-log2(X+1)
gtab <- GEOquery::Table(gpl)
fd <- gtab[match(rownames(X), gtab$ID),,drop=FALSE]
symcol <- grep("symbol", names(fd), ignore.case=TRUE, value=TRUE)[1]; cat("symbol col:",symcol,"\n")
sym <- as.character(fd[[symcol]])
# Canonical cross-platform score: keep the maximum-mean probe per symbol,
# z-score each gene within this cohort, then mean(up)-mean(down).
keep <- sym!="" & !is.na(sym); X <- X[keep,,drop=FALSE]; sym <- sym[keep]
ord <- order(-rowMeans(X)); X <- X[ord,,drop=FALSE]; sym <- sym[ord]
dedup <- !duplicated(sym); X <- X[dedup,,drop=FALSE]; rownames(X) <- sym[dedup]
z <- t(scale(t(X)))
score <- colMeans(z[intersect(up_sym,rownames(z)),,drop=FALSE],na.rm=TRUE) -
         colMeans(z[intersect(dn_sym,rownames(z)),,drop=FALSE],na.rm=TRUE)
grade <- factor(pd[["who grade:ch1"]]); rf <- suppressWarnings(as.integer(pd[["recurrence_frequency:ch1"]]))
rec <- factor(ifelse(rf>0,"Recur","NoRecur")); 
tts <- suppressWarnings(as.numeric(pd[["tts:ch1"]])); vital<-pd[["vital status:ch1"]]; ev<-ifelse(vital=="Deceased",1,ifelse(vital=="Alive",0,NA))
cat("\n[VAL1] aggressiveness score by WHO grade (GSE16581, independent platform):\n"); print(round(tapply(score,grade,mean),3))
cat("ANOVA score~grade p =", signif(anova(lm(score~grade))$`Pr(>F)`[1],3),"\n")
cat("\n[VAL2] score by recurrence:\n"); print(round(tapply(score,rec,mean,na.rm=TRUE),3))
cat("t-test Recur vs NoRecur p =", signif(t.test(score~rec)$p.value,3),"\n")
ok<-!is.na(tts)&!is.na(ev)
cat("\n[VAL3] survival (OS) Cox ~ aggressiveness score (n=",sum(ok),", events=",sum(ev[ok]),"):\n")
cx<-coxph(Surv(tts,ev)~score, subset=ok); print(summary(cx)$coefficients)
saveRDS(list(score=score,grade=grade,rec=rec,cox=cx,up_sym=up_sym,dn_sym=dn_sym),"results/deg/GSE16581_program_validation.rds")
cat("\nsaved validation -> results/deg/GSE16581_program_validation.rds\n")
