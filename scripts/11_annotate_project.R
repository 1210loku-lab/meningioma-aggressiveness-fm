# 11_annotate_project.R — 细胞类型注释 + 侵袭性程序细胞分辨投影（reviewer-grade）
suppressMessages({library(Seurat); library(dplyr); library(ggplot2); library(patchwork)})
set.seed(42)
o <- readRDS("results/scrna/GSE183655_processed.rds")
DefaultAssay(o) <- "RNA"

## 1. canonical markers per lineage（脑膜瘤 TME）
mk <- list(
  Immune=c("PTPRC"),
  Myeloid=c("CD68","C1QB","C1QA","LYZ","AIF1","CD163","ITGAM"),
  Tcell=c("CD3D","CD3E","CD2","IL7R","CD8A"),
  Bplasma=c("CD79A","MS4A1","IGHG1","MZB1"),
  Endothelial=c("PECAM1","VWF","CLDN5","FLT1"),
  Mural=c("RGS5","ACTA2","PDGFRB","NOTCH3"),
  Fibroblast=c("DCN","LUM","COL1A1","COL1A2"),
  Meningioma=c("PTGDS","CLDN1","NNAT","CRABP2","AGT","S100B"))
for(nm in names(mk)){ g<-intersect(mk[[nm]],rownames(o)); if(length(g)) o<-AddModuleScore(o,list(g),name=paste0("sc_",nm)) }
sccols <- grep("^sc_",colnames(o@meta.data),value=TRUE)

## 2. cluster-level: marker scores + patient mixing (肿瘤细胞常按患者聚, 免疫跨患者混)
cl <- levels(Idents(o))
clmean <- sapply(sccols, function(c) tapply(o@meta.data[[c]], Idents(o), mean))
pmix <- sapply(cl, function(k){ p<-o$patient[Idents(o)==k]; t<-sort(table(p),decreasing=TRUE); as.numeric(t[1])/sum(t) }) # 最大患者占比
# 指派：每群取最高 lineage 评分
lin <- sub("^sc_","",sub("[0-9]+$","",sccols))
assign_cl <- sapply(cl, function(k){
  v <- clmean[k,]; names(v)<-lin
  # 免疫细分：若 Immune 高，再在 Myeloid/T/B 里取最高
  top <- names(which.max(v))
  top
})
# 细化免疫
for(k in cl){ if(assign_cl[k]=="Immune"){ sub<-clmean[k,grep("Myeloid|Tcell|Bplasma",lin)]; assign_cl[k]<-sub("[0-9]+$","",sub("^sc_","",names(which.max(sub)))) } }
ct <- assign_cl[as.character(Idents(o))]; names(ct) <- colnames(o)
o$celltype <- factor(ct)
cat("=== cell type counts ===\n"); print(table(o$celltype))
cat("\n=== cluster patient-purity (高=肿瘤样) ===\n"); print(round(sort(pmix,decreasing=TRUE),2))

## 3. 侵袭性程序 by 细胞类型
cat("\n=== AggrScore by cell type ===\n"); print(round(sort(tapply(o$AggrScore,o$celltype,mean),decreasing=TRUE),3))
write.csv(data.frame(celltype=names(tapply(o$AggrScore,o$celltype,mean)),
                     AggrScore=as.numeric(tapply(o$AggrScore,o$celltype,mean)),
                     n=as.numeric(table(o$celltype)[names(tapply(o$AggrScore,o$celltype,mean))])),
          "results/scrna/aggr_score_by_celltype.csv", row.names=FALSE)

## 4. 肿瘤细胞内：侵袭性程序梯度（核心 — 是否存在高侵袭肿瘤态）
tum <- subset(o, celltype=="Meningioma")
cat("\n=== tumor cells:",ncol(tum),"; AggrScore quantiles ===\n"); print(round(quantile(tum$AggrScore,c(.1,.5,.9)),3))
# 肿瘤细胞按 AggrScore 三分，看高 vs 低 DEG（驱动基因）
tum$aggr_grp <- cut(tum$AggrScore, quantile(tum$AggrScore,c(0,.33,.67,1)), labels=c("Low","Mid","High"), include.lowest=TRUE)
Idents(tum)<-"aggr_grp"
deg <- tryCatch(FindMarkers(tum,ident.1="High",ident.2="Low",only.pos=FALSE,logfc.threshold=0.25), error=function(e)NULL)
if(!is.null(deg)){ write.csv(head(deg[order(deg$p_val_adj),],100),"results/scrna/tumor_highAggr_vs_low_DEG.csv"); cat("top tumor High-vs-Low genes:\n"); print(head(rownames(deg[order(deg$p_val_adj),]),20)) }

## 5. 关键图：UMAP cell type + AggrScore + by compartment
p1<-DimPlot(o,group.by="celltype",label=TRUE,repel=TRUE)+ggtitle("Cell types")
p2<-FeaturePlot(o,"AggrScore")+ggtitle("Aggressiveness program score")
p3<-DimPlot(o,group.by="compartment")+ggtitle("Compartment")
ggsave("results/scrna/fig_umap_celltype_aggr.png", p1|p2|p3, width=18,height=5.5,dpi=150)
ggsave("results/scrna/fig_umap_celltype_aggr.pdf", p1|p2|p3, width=18,height=5.5)
saveRDS(o, "results/scrna/GSE183655_annotated.rds")
cat("\nsaved GSE183655_annotated.rds + aggr_score_by_celltype.csv + figs\n")
