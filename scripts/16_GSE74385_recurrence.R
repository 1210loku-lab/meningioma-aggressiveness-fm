# 16_GSE74385_recurrence.R — 独立 Illumina 队列复发外验：侵袭性程序在 NR/R/M 及 grade-I 特异复发上的分离
suppressMessages({library(GEOquery); library(dplyr)})
set.seed(42)
raw<-"data/raw/GSE74385"; outdir<-"results/deg"; dir.create(outdir,showWarnings=FALSE,recursive=TRUE)

## 1) pheno（subtype）从本地 series matrix 解析
gse<-getGEO(filename=file.path(raw,"GSE74385_series_matrix.txt.gz"), getGPL=FALSE)
ph<-pData(gse)
subt<-ph[,grep("subtype",colnames(ph))[1]]
subt<-gsub("subtype: ","",subt)
ttl<-as.character(ph$title)
cat("samples:",length(ttl)," subtypes:\n"); print(table(subt))
meta<-data.frame(title=ttl, subtype=subt,
                 grade=substr(subt,1,1),
                 outcome=gsub("[0-9]","",subt), stringsAsFactors=FALSE)

## 2) 表达矩阵（去 Detection Pval 列）
mx<-read.delim(gzfile(file.path(raw,"GSE74385_normalized.txt.gz")),check.names=FALSE,stringsAsFactors=FALSE)
rownames(mx)<-mx[,1]; mx<-mx[,-1]
mx<-mx[, !grepl("Detection",colnames(mx)), drop=FALSE]   # 仅样本表达列
cat("expr matrix:",nrow(mx),"probes x",ncol(mx),"samples\n")
stopifnot(all(meta$title %in% colnames(mx)))
mx<-mx[, meta$title, drop=FALSE]

## 3) ILMN 探针 -> symbol（优先 illuminaHumanv4.db，回退 GPL soft）
if(requireNamespace("illuminaHumanv4.db",quietly=TRUE)){
  suppressMessages(library(illuminaHumanv4.db))
  sy<-AnnotationDbi::mapIds(illuminaHumanv4.db, keys=rownames(mx), column="SYMBOL", keytype="PROBEID", multiVals="first")
  map<-data.frame(probe=names(sy), sym=as.character(sy), stringsAsFactors=FALSE)
} else {
  gpl<-getGEO("GPL10558"); tab<-Table(gpl)
  symcol<-intersect(c("Symbol","ILMN_Gene","Gene_Symbol"),colnames(tab))[1]
  map<-data.frame(probe=as.character(tab$ID), sym=as.character(tab[[symcol]]), stringsAsFactors=FALSE)
}
map<-map[map$sym!="" & !is.na(map$sym),]
mx$probe<-rownames(mx)
mxm<-merge(map,mx,by="probe")
# 每基因取最大平均表达探针
mxm$mean<-rowMeans(sapply(mxm[,meta$title],as.numeric))
mxm<-mxm[order(-mxm$mean),]; mxm<-mxm[!duplicated(mxm$sym),]
expr<-as.matrix(sapply(mxm[,meta$title],as.numeric)); rownames(expr)<-mxm$sym
cat("genes after collapse:",nrow(expr),"\n")

## 4) 投影侵袭性程序（z-score within cohort，up - dn）
P<-readRDS("results/deg/GSE16581_program_validation.rds")
up<-intersect(P$up_sym,rownames(expr)); dn<-intersect(P$dn_sym,rownames(expr))
cat("program genes mapped: up",length(up),"dn",length(dn),"\n")
z<-t(scale(t(expr)))
score<-colMeans(z[up,,drop=FALSE],na.rm=TRUE)-colMeans(z[dn,,drop=FALSE],na.rm=TRUE)
meta$score<-score[meta$title]

## 5) 统计
meta$recur<-ifelse(meta$outcome=="NR","NonRecur",
                   ifelse(meta$outcome %in% c("R","M"),"Recur/Malig",NA_character_))
meta$recur<-factor(meta$recur,levels=c("NonRecur","Recur/Malig"))
cat("\n=== AggrScore by outcome (all grades) ===\n"); print(round(tapply(meta$score,meta$outcome,mean),3))
cat("known-outcome recur(R/M) vs NR Wilcoxon p=",signif(wilcox.test(score~recur,meta)$p.value,3),"\n")
av<-summary(aov(score~outcome,meta)); cat("3-group(NR/R/M) ANOVA p=",signif(av[[1]]$`Pr(>F)`[1],3),"\n")

## ⭐ grade-I 特异复发（临床关键，欠功效）
g1<-meta[meta$grade=="1",];
cat("\n=== ⭐ grade-I only: 1NR(n=",sum(g1$recur=="NonRecur"),") vs 1R/1M(n=",sum(g1$recur!="NonRecur"),") ===\n",sep="")
print(round(tapply(g1$score,g1$recur,mean),3))
w1<-wilcox.test(score~recur,g1); cat("grade-I recur Wilcoxon p=",signif(w1$p.value,3)," (underpowered, suggestive)\n")
# grade trend
cat("\nAggrScore by grade:\n"); print(round(tapply(meta$score,meta$grade,mean),3))
cat("grade Spearman rho=",signif(cor(as.numeric(meta$grade),meta$score,method="spearman"),3),"\n")

write.csv(meta[,c("title","subtype","grade","outcome","recur","score")],
          file.path(outdir,"GSE74385_program_recurrence.csv"),row.names=FALSE)
saveRDS(list(meta=meta,up=up,dn=dn),file.path(outdir,"GSE74385_program_recurrence.rds"))
cat("\nALL-DONE saved GSE74385_program_recurrence.csv/.rds\n")
