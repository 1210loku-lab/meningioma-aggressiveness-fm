# 25_signature_benchmark.R — P1 新颖性：侵袭性程序 vs 已知生物学签名(增殖/OXPHOS) vs 随机零模型
# 检验数据驱动的程序在复发预测上是否优于"已知高级别生物学"(增殖)及通用 OXPHOS，并超过随机同尺寸基因集
suppressMessages({library(GEOquery); library(pROC)})
set.seed(42); sink("results/deg/P1_signature_benchmark.txt", split=TRUE)

## 已知生物学签名（offline curated canonical sets，避免 msigdbr 联网依赖；可引用 KEGG/cell-cycle 文献）
sig_list <- list(
  Proliferation = c("MKI67","TOP2A","CCNB1","CCNB2","CCNA2","CDK1","CDC20","BUB1","BUB1B","AURKA",
    "AURKB","PLK1","FOXM1","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CENPA","CENPE","CENPF",
    "KIF11","KIF23","TYMS","RRM2","BIRC5","UBE2C","NUSAP1","ASPM","TPX2","CCNE2","E2F1","CDC45","GINS2"),
  OXPHOS = c("NDUFA1","NDUFA2","NDUFA4","NDUFB1","NDUFB2","NDUFS1","NDUFS2","NDUFV1","SDHA","SDHB",
    "SDHC","SDHD","UQCRC1","UQCRC2","UQCRB","UQCRH","CYC1","COX4I1","COX5A","COX5B","COX6A1","COX6B1",
    "COX6C","COX7A2","COX7B","COX8A","ATP5F1A","ATP5F1B","ATP5F1C","ATP5MC1","ATP5PB","ATP5PO")
)
P <- readRDS("results/deg/GSE16581_program_validation.rds")
prog_up <- P$up_sym; prog_dn <- P$dn_sym
cat("program: up",length(prog_up)," dn",length(prog_dn),"\n")
for(n in names(sig_list)) cat(n,":",length(sig_list[[n]]),"genes\n")

score_dir <- function(expr, up, dn=NULL){ # z within cohort, mean(up)-mean(dn)
  z <- t(scale(t(expr)))
  u <- intersect(up, rownames(z)); s <- colMeans(z[u,,drop=FALSE],na.rm=TRUE)
  if(!is.null(dn)){ d <- intersect(dn, rownames(z)); s <- s - colMeans(z[d,,drop=FALSE],na.rm=TRUE) }
  s
}
auc_rec <- function(score, rec){ as.numeric(pROC::auc(rec, score, quiet=TRUE, direction="<")) }

bench_cohort <- function(expr, rec, label){
  cat("\n========== ",label," (n=",ncol(expr)," rec events=",sum(rec),") ==========\n")
  prog_s <- score_dir(expr, prog_up, prog_dn)
  a_prog <- auc_rec(prog_s, rec)
  cat(sprintf("Aggressiveness program recurrence AUC = %.3f\n", a_prog))
  for(n in names(sig_list)){ a <- auc_rec(score_dir(expr, sig_list[[n]]), rec); cat(sprintf("  %-20s AUC = %.3f\n", n, a)) }
  # 随机零模型：同尺寸(=程序up数) 随机基因集 1000 次
  allg <- rownames(expr); k <- length(intersect(prog_up,allg))
  rnd <- replicate(1000, auc_rec(score_dir(expr, sample(allg,k)), rec))
  emp_p <- mean(rnd >= a_prog)
  cat(sprintf("  Random %d-gene null: mean AUC=%.3f, 95%%=%.3f ; program empirical p=%.3f\n",
              k, mean(rnd), quantile(rnd,.95), emp_p))
  invisible(list(prog=a_prog, rnd=rnd))
}

## --- GSE16581 (GPL570) ---
es <- getGEO(filename="data/raw/GSE16581_series_matrix.txt.gz", getGPL=FALSE)
gpl <- getGEO(filename="data/raw/GPL570.soft.gz")
X <- Biobase::exprs(es); pd<-Biobase::pData(es)
if(max(X, na.rm=TRUE) > 100) X <- log2(X + 1)
gtab <- GEOquery::Table(gpl)
fd <- gtab[match(rownames(X), gtab$ID), , drop=FALSE]
symcol <- grep("symbol", names(fd), ignore.case=TRUE, value=TRUE)[1]
sym <- as.character(fd[[symcol]]); keep <- sym!="" & !is.na(sym)
X <- X[keep,]; sym <- sym[keep]
# collapse to symbol (max mean)
ord <- order(-rowMeans(X)); X<-X[ord,]; sym<-sym[ord]; X<-X[!duplicated(sym),]; rownames(X)<-sym[!duplicated(sym)]
rf <- suppressWarnings(as.integer(pd[["recurrence_frequency:ch1"]])); rec <- ifelse(rf>0,1,0)
ok <- !is.na(rec); bench_cohort(X[,ok], rec[ok], "GSE16581")

## --- GSE74385 (Illumina, 复用 16 的 collapse) ---
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
ok2<-!is.na(rec2)
bench_cohort(ev[,ok2,drop=FALSE], rec2[ok2], "GSE74385 known outcomes")

cat("\n解读：若程序 AUC 同时 > 增殖(E2F/G2M)、通用 OXPHOS、且随机零模型 emp p<0.05，\n说明该数据驱动程序携带超出'已知高级别增殖/通用代谢'的复发预测信息(新颖性)。\n")
sink(); cat("\nALL-DONE -> results/deg/P1_signature_benchmark.txt\n")
