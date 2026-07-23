# 11b_reannotate.R — 修正注释：脑膜瘤肿瘤细胞=非免疫/内皮/周细胞的患者特异间充质细胞
suppressMessages({library(Seurat); library(dplyr); library(ggplot2); library(patchwork)}); set.seed(42)
o <- readRDS("results/scrna/GSE183655_processed.rds"); DefaultAssay(o)<-"RNA"
mk <- list(Immune="PTPRC", Myeloid=c("CD68","C1QB","C1QA","LYZ","AIF1","CD163"),
           Tcell=c("CD3D","CD3E","CD2","IL7R"), Bplasma=c("CD79A","MS4A1","IGHG1","MZB1"),
           Endothelial=c("PECAM1","VWF","CLDN5","FLT1"), Mural=c("RGS5","ACTA2","PDGFRB","NOTCH3"))
for(nm in names(mk)){ g<-intersect(mk[[nm]],rownames(o)); if(length(g)) o<-AddModuleScore(o,list(g),name=paste0("sc_",nm)) }
cl<-levels(Idents(o)); sccols<-paste0("sc_",names(mk),"1")
cm<-sapply(sccols,function(c) tapply(o@meta.data[[c]],Idents(o),mean)); colnames(cm)<-names(mk)
pmix<-sapply(cl,function(k){p<-o$patient[Idents(o)==k];max(table(p))/sum(table(p))})
# 规则：免疫(Immune>0.1 且高于内皮/周细胞)→细分；否则内皮/周细胞；否则=Meningioma 肿瘤
assign_cl<-sapply(cl,function(k){
  v<-cm[k,]
  if(v["Immune"]>0.1 & v["Immune"]>v["Endothelial"] & v["Immune"]>v["Mural"]){
    names(which.max(v[c("Myeloid","Tcell","Bplasma")]))
  } else if(v["Endothelial"]>0.3 & v["Endothelial"]>v["Mural"]) "Endothelial"
  else if(v["Mural"]>0.3) "Mural"
  else "Meningioma"
})
ct<-assign_cl[as.character(Idents(o))]; names(ct)<-colnames(o); o$celltype<-factor(ct)
cat("=== corrected cell types ===\n"); print(table(o$celltype))
cat("=== tumor purity check (Meningioma clusters should be patient-specific) ===\n")
cat("Meningioma clusters mean patient-purity:", round(mean(pmix[assign_cl=="Meningioma"]),2),
    "| Immune clusters:", round(mean(pmix[assign_cl %in% c("Myeloid","Tcell","Bplasma")]),2),"\n")
cat("\n=== AggrScore by corrected cell type ===\n"); print(round(sort(tapply(o$AggrScore,o$celltype,mean),decreasing=TRUE),3))
write.csv(data.frame(celltype=names(table(o$celltype)),n=as.numeric(table(o$celltype)),
          AggrScore=round(as.numeric(tapply(o$AggrScore,o$celltype,mean)[names(table(o$celltype))]),3)),
          "results/scrna/aggr_score_by_celltype_v2.csv",row.names=FALSE)
p1<-DimPlot(o,group.by="celltype",label=TRUE,repel=TRUE)+ggtitle("Meningioma scRNA cell types (GSE183655, n=63,628)")
p2<-FeaturePlot(o,"AggrScore")+ggtitle("Bulk aggressiveness program score")
p3<-DimPlot(o,group.by="patient")+ggtitle("Patient")
fig <- (p1|p2|p3)+plot_annotation(tag_levels="a") &
  theme(text=element_text(family="Arial"), plot.tag=element_text(family="Arial",face="bold",size=14))
ggsave("results/scrna/fig_umap_celltype_v2.png",fig,width=19,height=5.5,dpi=300)
ggsave("results/scrna/fig_umap_celltype_v2.pdf",fig,width=19,height=5.5)
saveRDS(o,"results/scrna/GSE183655_annotated.rds"); cat("\nsaved GSE183655_annotated.rds + v2 outputs\n")
