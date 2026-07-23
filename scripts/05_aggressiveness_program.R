suppressMessages({library(DESeq2)}); set.seed(42)
A <- readRDS("data/raw/GSE136661_assembled.rds"); M <- A$counts; key <- A$key
rownames(key) <- key$gsm; key <- key[colnames(M),]
key$grade <- factor(gsub("WHO ","",key[["pathology:ch1"]]), levels=c("I","II","III"))
key$cohort <- key[["cohort:ch1"]]
cat("gene id sample:", paste(head(rownames(M),3),collapse=" | "), "\n")
# Define aggressiveness program on PRIMARY tumors (Discovery+Validation), WHO I vs II
prim <- key$cohort %in% c("Discovery","Validation") & key$grade %in% c("I","II")
Mp <- M[, prim]; kp <- droplevels(key[prim,])
cat("primary subset:", ncol(Mp), "samples |", paste(names(table(kp$grade)),table(kp$grade),collapse=" "), "\n")
storage.mode(Mp) <- "integer"
Mp <- Mp[rowSums(Mp>=5) >= 10, ]   # filter low counts
cat("genes after filter:", nrow(Mp), "\n")
dds <- DESeqDataSetFromMatrix(Mp, kp, ~ grade)
dds <- DESeq(dds, quiet=TRUE)
res <- as.data.frame(results(dds, contrast=c("grade","II","I")))
res <- res[order(res$padj), ]
nsig <- sum(res$padj < 0.05, na.rm=TRUE)
cat(sprintf("\n[AGGRESSIVENESS PROGRAM] WHO II vs I (primaries): padj<0.05 = %d genes\n", nsig))
cat("top15:\n"); print(head(res[,c("log2FoldChange","padj")], 15))
# program = top up/down genes
up <- rownames(head(res[res$log2FoldChange>1 & !is.na(res$padj) & res$padj<0.05,], 100))
dn <- rownames(head(res[res$log2FoldChange< -1 & !is.na(res$padj) & res$padj<0.05,], 100))
cat(sprintf("\nprogram size: up=%d dn=%d (|log2FC|>1, padj<0.05)\n", length(up), length(dn)))
# Validate: program score (vst) rises with grade incl WHO III (Recurrence cohort)?
vsd <- assay(vst(DESeqDataSetFromMatrix(M[rownames(Mp),][rowSums(M[rownames(Mp),]>=5)>=10,, drop=FALSE], key, ~1), blind=TRUE))
upx <- intersect(up, rownames(vsd)); score <- colMeans(vsd[upx,,drop=FALSE]) - colMeans(vsd[intersect(dn,rownames(vsd)),,drop=FALSE])
cat("\n[VALIDATION] aggressiveness score (up-dn, vst) by WHO grade across ALL 160:\n")
print(round(tapply(score, key$grade, mean), 3))
cat("score by cohort:\n"); print(round(tapply(score, key$cohort, mean), 3))
saveRDS(list(res=res, up=up, dn=dn, score=score, key=key), "results/deg/GSE136661_aggressiveness_program.rds")
write.csv(res, "results/deg/GSE136661_aggressiveness_WHO_II_vs_I.csv")
cat("\nsaved program -> results/deg/\n")
