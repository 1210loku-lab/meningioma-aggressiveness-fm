# 17a_export_tumor_for_geneformer.R — 导出 grade 分层肿瘤细胞 subsample 供 Geneformer tokenize
# 输出 MatrixMarket counts + genes(ENSG) + cells(grade,n_counts)，Python 端组 AnnData
suppressMessages({library(Seurat); library(Matrix)})
set.seed(42)
o<-readRDS("results/scrna/GSE206647_processed.rds")
o<-subset(o,celltype=="Meningioma" & grade %in% c("I","II","III"))
cat("tumor cells (graded):",ncol(o)," by grade:\n"); print(table(o$grade))
# 每 grade subsample（CPU 推理可承受）
ncell<-600
cells<-unlist(lapply(c("I","II","III"),function(g){i<-which(o$grade==g); if(length(i)>ncell) sample(i,ncell) else i}))
o<-o[,cells]; cat("subsampled:",ncol(o)," by grade:\n"); print(table(o$grade))

cts<-GetAssayData(o,assay="RNA",layer="counts")
# symbol -> ENSG
suppressMessages(library(org.Hs.eg.db))
ens<-AnnotationDbi::mapIds(org.Hs.eg.db,keys=rownames(cts),column="ENSEMBL",keytype="SYMBOL",multiVals="first")
keep<-!is.na(ens); cts<-cts[keep,]; ens<-ens[keep]
# 去重 ENSG（保留行和最大）
ord<-order(-Matrix::rowSums(cts)); cts<-cts[ord,]; ens<-ens[ord]
dup<-duplicated(ens); cts<-cts[!dup,]; ens<-ens[!dup]
cat("genes with unique ENSG:",nrow(cts),"\n")

outdir<-"results/scrna/gf_export"; dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
Matrix::writeMM(cts,file.path(outdir,"counts.mtx"))   # genes x cells
write.csv(data.frame(symbol=rownames(cts),ensembl_id=ens),file.path(outdir,"genes.csv"),row.names=FALSE)
write.csv(data.frame(cell=colnames(cts),grade=o$grade,gsm=o$gsm,
                     n_counts=Matrix::colSums(cts),AggrScore=o$AggrScore),
          file.path(outdir,"cells.csv"),row.names=FALSE)
cat("ALL-DONE exported to",outdir,"\n")
