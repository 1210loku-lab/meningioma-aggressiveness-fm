# Rebuild Figure 2 from raw 10x matrices without the offloaded 1-GB Seurat RDS.
# Writes only compact plotting metadata and publication figures.
suppressMessages({
  library(Seurat); library(harmony); library(Matrix); library(ggplot2); library(patchwork)
})
set.seed(42)
options(future.globals.maxSize=8e9)
root <- normalizePath(".")
raw <- file.path(root,"data","raw","10x")
gsm <- read.csv(file.path(root,"data","raw","GSE183655_gsm_list.csv"),stringsAsFactors=FALSE)
tag <- sub("_barcodes.tsv.gz$","",basename(gsm$supp))
tag <- sub("^GSM[0-9]+_","",tag)
names(tag) <- gsm$gsm

read_one <- function(g){
  d <- file.path(raw,g)
  m <- ReadMtx(file.path(d,"matrix.mtx.gz"),file.path(d,"barcodes.tsv.gz"),
               file.path(d,"features.tsv.gz"),feature.column=2)
  colnames(m) <- paste0(g,"_",colnames(m))
  o <- CreateSeuratObject(m,project=tag[g],min.cells=3,min.features=200)
  o$sample <- unname(tag[g]); o$gsm <- g; o
}
objs <- lapply(gsm$gsm,read_one)
obj <- merge(objs[[1]],objs[-1])
obj$compartment <- ifelse(grepl("Dura",obj$sample),"Dura",
                          ifelse(grepl("BTI",obj$sample),"BTI","Tumour"))
obj$patient <- sub("-(Dura|BTI)$","",obj$sample)
obj[["percent.mt"]] <- PercentageFeatureSet(obj,pattern="^MT-")
obj <- subset(obj,subset=nFeature_RNA>200 & nFeature_RNA<7500 & percent.mt<20)
if("JoinLayers" %in% getNamespaceExports("SeuratObject")) obj <- JoinLayers(obj)
obj <- NormalizeData(obj,verbose=FALSE) |>
  FindVariableFeatures(nfeatures=2000,verbose=FALSE) |>
  ScaleData(verbose=FALSE) |>
  RunPCA(npcs=30,verbose=FALSE)
obj <- RunHarmony(obj,group.by.vars="sample",verbose=FALSE)
obj <- RunUMAP(obj,reduction="harmony",dims=1:30,seed.use=42,verbose=FALSE) |>
  FindNeighbors(reduction="harmony",dims=1:30,verbose=FALSE) |>
  FindClusters(resolution=0.5,verbose=FALSE)

program <- read.csv(file.path(root,"docs","Table_S2_aggressiveness_program_genes.csv"),stringsAsFactors=FALSE)
up <- intersect(program$symbol[grepl("^up",program$direction)],rownames(obj))
dn <- intersect(program$symbol[grepl("^down",program$direction)],rownames(obj))
obj <- AddModuleScore(obj,list(up),name="AggrUp",seed=42)
obj <- AddModuleScore(obj,list(dn),name="AggrDn",seed=42)
obj$AggrScore <- obj$AggrUp1-obj$AggrDn1

mk <- list(Immune="PTPRC",Myeloid=c("CD68","C1QB","C1QA","LYZ","AIF1","CD163"),
           Tcell=c("CD3D","CD3E","CD2","IL7R"),Bplasma=c("CD79A","MS4A1","IGHG1","MZB1"),
           Endothelial=c("PECAM1","VWF","CLDN5","FLT1"),Mural=c("RGS5","ACTA2","PDGFRB","NOTCH3"))
for(nm in names(mk)){
  g <- intersect(mk[[nm]],rownames(obj)); if(length(g)) obj <- AddModuleScore(obj,list(g),name=paste0("sc_",nm),seed=42)
}
cl <- levels(Idents(obj)); sccols <- paste0("sc_",names(mk),"1")
cm <- sapply(sccols,function(c) tapply(obj@meta.data[[c]],Idents(obj),mean)); colnames(cm)<-names(mk)
assign_cl <- sapply(cl,function(k){
  v<-cm[k,]
  if(v["Immune"]>0.1 && v["Immune"]>v["Endothelial"] && v["Immune"]>v["Mural"])
    names(which.max(v[c("Myeloid","Tcell","Bplasma")]))
  else if(v["Endothelial"]>0.3 && v["Endothelial"]>v["Mural"]) "Endothelial"
  else if(v["Mural"]>0.3) "Mural" else "Meningioma"
})
obj$celltype <- factor(assign_cl[as.character(Idents(obj))])

um <- Embeddings(obj,"umap")
plotdata <- data.frame(cell=colnames(obj),UMAP1=um[,1],UMAP2=um[,2],
                       celltype=obj$celltype,patient=obj$patient,compartment=obj$compartment,
                       AggrScore=obj$AggrScore,stringsAsFactors=FALSE)
outdir <- file.path(root,"results","audit_submission"); dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
write.csv(plotdata,file.path(outdir,"GSE183655_Figure2_source_data.csv"),row.names=FALSE)
write.csv(data.frame(celltype=names(table(obj$celltype)),n=as.numeric(table(obj$celltype)),
                     AggrScore=as.numeric(tapply(obj$AggrScore,obj$celltype,mean)[names(table(obj$celltype))])),
          file.path(root,"results","scrna","aggr_score_by_celltype_v2.csv"),row.names=FALSE)

theme_set(theme_classic(base_size=11,base_family="Arial"))
p1 <- DimPlot(obj,group.by="celltype",label=TRUE,repel=TRUE,raster=TRUE)+
  labs(title=sprintf("GSE183655 cell types (%s QC cells)",format(ncol(obj),big.mark=",")),x="UMAP 1",y="UMAP 2")+
  theme(legend.position="none")
p2 <- FeaturePlot(obj,"AggrScore",raster=TRUE)+
  labs(title="Bulk-derived program score (descriptive)",x="UMAP 1",y="UMAP 2")
p3 <- DimPlot(obj,group.by="patient",raster=TRUE)+
  labs(title="Patient of origin",x="UMAP 1",y="UMAP 2")
fig <- (p1|p2|p3)+plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(family="Arial",face="bold",size=14))
ggsave(file.path(root,"results","scrna","fig_umap_celltype_v2.png"),fig,width=19,height=5.5,dpi=300,bg="white")
ggsave(file.path(root,"results","scrna","fig_umap_celltype_v2.pdf"),fig,width=19,height=5.5,bg="white")
cat("Figure 2 rebuilt from raw matrices; QC cells=",ncol(obj),"\n",sep="")
print(table(obj$celltype))
