# 10_scrna_process.R — GSE183655 脑膜瘤 scRNA 处理 + 侵袭性程序投影
# 前置：data/raw/GSE183655_RAW.tar 完整下载（741MB）。staged，下载完成后运行。
suppressMessages({library(Seurat); library(Matrix); library(harmony); library(dplyr); library(ggplot2)})
set.seed(42); options(future.globals.maxSize=8e9)
raw <- "data/raw"; full <- file.path(raw,"10x_full")

## 1. 读 10x_full/（解压自 RAW.tar，文件名 GSM_tag_{barcodes,features,matrix}），用 ReadMtx
mtxs <- list.files(full, pattern="_matrix.mtx.gz$", full.names=TRUE)
stubs <- sub("_matrix.mtx.gz$","",basename(mtxs))   # e.g. GSM5567093_MSC1
cat("samples found:", length(stubs), "\n"); print(stubs)
read_one <- function(stub){
  pre <- file.path(full, stub)
  m <- ReadMtx(mtx=paste0(pre,"_matrix.mtx.gz"),
               cells=paste0(pre,"_barcodes.tsv.gz"),
               features=paste0(pre,"_features.tsv.gz"), feature.column=2)
  colnames(m) <- paste0(stub,"_",colnames(m))
  tag <- sub("^GSM[0-9]+_","",stub)   # MSC1 / MSC4-Dura / MSC5-BTI
  o <- CreateSeuratObject(m, project=tag, min.cells=3, min.features=200)
  o$sample <- tag; o$gsm <- sub("_.*","",stub); o
}
objs <- lapply(stubs, function(s) tryCatch(read_one(s), error=function(e){cat("skip",s,":",conditionMessage(e),"\n");NULL}))
objs <- objs[!sapply(objs,is.null)]
cat("loaded samples:", length(objs), "\n")
obj <- if(length(objs)>1) merge(objs[[1]], objs[-1]) else objs[[1]]
# 区室/患者元数据：MSC{n} / -Dura / -BTI
obj$compartment <- ifelse(grepl("Dura",obj$sample),"Dura", ifelse(grepl("BTI",obj$sample),"BTI","Tumor"))
obj$patient <- sub("-(Dura|BTI)$","", obj$sample)
cat("cells before QC:", ncol(obj), "\n")

## 2. QC（脑膜瘤肿瘤细胞容忍较高 mt）
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^MT-")
obj <- subset(obj, subset = nFeature_RNA>200 & nFeature_RNA<7500 & percent.mt<20)
cat("cells after QC:", ncol(obj), "\n")

## 3. 标准化 + HVG + PCA + harmony 整合（按 sample）
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) obj <- JoinLayers(obj)  # Seurat v5 合并 counts 层
obj <- NormalizeData(obj) |> FindVariableFeatures(nfeatures=2000) |> ScaleData() |> RunPCA(npcs=30, verbose=FALSE)
obj <- RunHarmony(obj, group.by.vars="sample", verbose=FALSE)
obj <- RunUMAP(obj, reduction="harmony", dims=1:30) |>
       FindNeighbors(reduction="harmony", dims=1:30) |> FindClusters(resolution=0.5)
cat("clusters:", nlevels(Idents(obj)), "\n")

## 4. 广谱注释（canonical markers）
mk <- list(Meningioma=c("MN1","NF2","PTGDS","CLDN1","EPCAM"), Endothelial=c("PECAM1","VWF","CLDN5"),
           Mural=c("RGS5","ACTA2","PDGFRB"), Macrophage=c("CD68","CD163","C1QB","LYZ"),
           Tcell=c("CD3D","CD3E","CD8A"), Bcell=c("CD79A","MS4A1"), Fibroblast=c("DCN","LUM","COL1A1"))
saveRDS(obj, "results/scrna/GSE183655_processed.rds")

## 5. 侵袭性程序投影（核心：哪些细胞态高表达侵袭程序）
P <- readRDS("results/deg/GSE16581_program_validation.rds")  # up_sym/dn_sym（symbol）
up <- intersect(P$up_sym, rownames(obj)); dn <- intersect(P$dn_sym, rownames(obj))
obj <- AddModuleScore(obj, features=list(up), name="AggrUp")
obj <- AddModuleScore(obj, features=list(dn), name="AggrDn")
obj$AggrScore <- obj$AggrUp1 - obj$AggrDn1
cat("\n[PROGRAM PROJECTION] mean AggrScore by cluster:\n"); print(round(tapply(obj$AggrScore, Idents(obj), mean),3))
write.csv(data.frame(cluster=names(tapply(obj$AggrScore,Idents(obj),mean)),
                     AggrScore=tapply(obj$AggrScore,Idents(obj),mean)),
          "results/scrna/aggr_score_by_cluster.csv")
saveRDS(obj, "results/scrna/GSE183655_processed.rds")
cat("\nsaved results/scrna/GSE183655_processed.rds\n")
