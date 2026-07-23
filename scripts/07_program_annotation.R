suppressMessages({library(clusterProfiler);library(org.Hs.eg.db);library(AnnotationDbi)})
P <- readRDS("results/deg/GSE136661_aggressiveness_program.rds")
res <- P$res
sym <- function(e) mapIds(org.Hs.eg.db,e,"SYMBOL","ENSEMBL")
# top program genes with symbols + log2FC
top_up <- head(res[res$log2FoldChange>0 & !is.na(res$padj),],20)
top_dn <- head(res[res$log2FoldChange<0 & !is.na(res$padj),],20)
cat("=== TOP UP (high-grade/aggressive) ===\n")
print(data.frame(symbol=sym(rownames(top_up)), log2FC=round(top_up$log2FoldChange,2), padj=signif(top_up$padj,2)))
cat("\n=== TOP DOWN (low-grade) ===\n")
print(data.frame(symbol=sym(rownames(top_dn)), log2FC=round(top_dn$log2FoldChange,2), padj=signif(top_dn$padj,2)))
# GO enrichment (BP) on up & dn programs (ENTREZ)
toentrez <- function(e) na.omit(mapIds(org.Hs.eg.db,e,"ENTREZID","ENSEMBL"))
univ <- toentrez(rownames(res))
for(set in c("up","dn")){
  genes <- toentrez(P[[set]])
  cat(sprintf("\n=== GO:BP enrichment — %s program (%d genes) ===\n", set, length(genes)))
  eg <- tryCatch(enrichGO(genes, OrgDb=org.Hs.eg.db, ont="BP", universe=univ, pAdjustMethod="BH", qvalueCutoff=0.1, readable=TRUE), error=function(e)NULL)
  if(!is.null(eg) && nrow(as.data.frame(eg))>0) print(head(as.data.frame(eg)[,c("Description","GeneRatio","p.adjust")],10)) else cat("(no significant GO terms)\n")
  if(!is.null(eg)) write.csv(as.data.frame(eg), sprintf("results/deg/GO_BP_%s_program.csv",set))
}
