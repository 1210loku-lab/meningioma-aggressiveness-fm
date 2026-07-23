suppressMessages({library(DESeq2); library(MCPcounter); library(ggplot2); library(org.Hs.eg.db)})
outdir <- "results/descriptive"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

cat("--- Bulk deconvolution using MCPcounter (GSE136661) ---\n")
A <- readRDS("data/raw/GSE136661_assembled.rds")
M <- A$counts; key <- A$key
rownames(key) <- key$gsm; key <- key[colnames(M),]
key$grade <- factor(gsub("WHO ", "", key[["pathology:ch1"]]), levels=c("I", "II", "III"))

M <- M[rowSums(M >= 5) >= 10, ]
dds <- DESeqDataSetFromMatrix(M, key, ~ 1)
vsd <- assay(vst(dds, blind=TRUE))
syms <- mapIds(org.Hs.eg.db, keys=rownames(vsd), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
vsd_sym <- vsd[!is.na(syms), ]
rownames(vsd_sym) <- syms[!is.na(syms)]
vsd_sym <- vsd_sym[!duplicated(rownames(vsd_sym)), ]

mcp <- MCPcounter.estimate(vsd_sym, featuresType="HUGO_symbols")
key$T_cells <- mcp["T cells",]
key$Fibroblasts <- mcp["Fibroblasts",]
key$Endothelial <- mcp["Endothelial cells",]
for (feature in c("T_cells", "Fibroblasts", "Endothelial")) {
  cat(feature, "by grade:\n")
  print(tapply(key[[feature]], key$grade, mean))
  cat("ANOVA p-value:", signif(anova(lm(key[[feature]] ~ key$grade))$`Pr(>F)`[1], 3), "\n")
}

write.csv(t(mcp), file.path(outdir, "19_GSE136661_mcpcounter_scores.csv"))
saveRDS(list(mcp=mcp, key=key), file.path(outdir, "19_GSE136661_mcpcounter.rds"))
pdf(file.path(outdir, "19_fig_mcpcounter.pdf"), width=8, height=4, family="Arial")
par(mfrow=c(1,3))
boxplot(T_cells ~ grade, key, main="T cells (MCPcounter)", ylab="Score", col="lightblue")
boxplot(Fibroblasts ~ grade, key, main="Fibroblasts (MCPcounter)", ylab="Score", col="lightgreen")
boxplot(Endothelial ~ grade, key, main="Endothelial (MCPcounter)", ylab="Score", col="lightpink")
dev.off()
