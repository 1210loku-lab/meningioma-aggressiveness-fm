suppressMessages({library(GEOquery)}); options(timeout=300)
dir <- "data/raw/GSE136661_ex"
fs <- list.files(dir, pattern="htseq.counts.txt.gz$", full.names=TRUE)
cat("count files:", length(fs), "\n")
# assemble counts matrix
rl <- lapply(fs, function(f){ d<-read.table(gzfile(f),sep="\t",col.names=c("gene","cnt")); setNames(d$cnt,d$gene) })
genes <- rl[[1]]; gn <- names(genes)
M <- sapply(rl, function(x) x[gn])
gsm <- sub("_.*","",basename(fs)); colnames(M) <- gsm
# drop HTSeq summary rows
M <- M[!grepl("^__", rownames(M)),]
cat("counts matrix:", nrow(M),"genes x",ncol(M),"samples\n")
# join pData
g <- getGEO("GSE136661",GSEMatrix=TRUE,getGPL=FALSE,destdir="data/raw"); pd<-Biobase::pData(g[[1]])
pd$gsm <- rownames(pd)
key <- data.frame(gsm=colnames(M))
key <- merge(key, pd[,c("gsm","cohort:ch1","pathology:ch1","title")], by="gsm", all.x=TRUE)
cat("\n== cohort x pathology ==\n"); print(table(key[["cohort:ch1"]], key[["pathology:ch1"]], useNA="always"))
cat("\n== title 'R' suffix (recurrence?) sample ==\n"); print(head(key$title, 12))
# does title encode recurrence? count titles containing R-like recurrence marks
saveRDS(list(counts=M, key=key, pd=pd), "data/raw/GSE136661_assembled.rds")
cat("\nsaved data/raw/GSE136661_assembled.rds\n")
