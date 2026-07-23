# 15_GSE206647_grade.R — grade 分层 scRNA：侵袭性程序在肿瘤细胞内是否随 grade 升高（核心闭环）
suppressMessages({library(Seurat); library(Matrix); library(harmony); library(dplyr); library(ggplot2); library(patchwork)})
set.seed(42); options(future.globals.maxSize=10e9)
raw<-"data/raw"; ex<-file.path(raw,"GSE206647_ex"); dir.create(ex,showWarnings=FALSE)

## 患者肿瘤 + 脑膜样本（排除 IOMM-LEE 细胞系 GSM62599*）
smp <- read.csv(file.path(raw,"GSE206647_samples.csv"),stringsAsFactors=FALSE)
smp <- smp[grepl("Mening",smp$src),]   # Meningoma + Meninges, 去细胞系
smp$grade2 <- ifelse(smp$src=="Meninges","Normal", gsub("grade ","",smp$grade))
cat("patient samples:",nrow(smp),"\n"); print(table(smp$grade2))

## 仅解压所需样本文件
patt <- paste(smp$gsm, collapse="|")
system(sprintf("tar -xf %s/GSE206647/GSE206647_RAW.tar -C %s", raw, ex))
mtxs <- list.files(ex, pattern="_matrix.mtx.gz$", full.names=TRUE)
mtxs <- mtxs[grepl(patt, basename(mtxs))]
cat("matrices to load:",length(mtxs),"\n")
read_one<-function(mx){
  stub<-sub("_matrix.mtx.gz$","",mx); gsm<-sub("_.*","",basename(stub))
  m<-ReadMtx(mtx=mx, cells=paste0(stub,"_barcodes.tsv.gz"), features=paste0(stub,"_features.tsv.gz"), feature.column=2)
  colnames(m)<-paste0(gsm,"_",colnames(m))
  o<-CreateSeuratObject(m,project=gsm,min.cells=3,min.features=200)
  gr<-smp$grade2[match(gsm,smp$gsm)]; o$grade<-gr; o$gsm<-gsm; o
}
objs<-lapply(mtxs,function(x)tryCatch(read_one(x),error=function(e){cat("skip",basename(x),conditionMessage(e),"\n");NULL}))
objs<-objs[!sapply(objs,is.null)]
o<-merge(objs[[1]],objs[-1]); rm(objs); gc()
o[["percent.mt"]]<-PercentageFeatureSet(o,pattern="^MT-")
o<-subset(o,subset=nFeature_RNA>200 & nFeature_RNA<7500 & percent.mt<20)
cat("cells after QC:",ncol(o)," | by grade:\n"); print(table(o$grade))
if("JoinLayers"%in%getNamespaceExports("SeuratObject")) o<-JoinLayers(o)
o<-NormalizeData(o)|>FindVariableFeatures(nfeatures=2000)|>ScaleData()|>RunPCA(npcs=30,verbose=FALSE)
o<-RunHarmony(o,group.by.vars="gsm",verbose=FALSE)
o<-RunUMAP(o,reduction="harmony",dims=1:30)|>FindNeighbors(reduction="harmony",dims=1:30)|>FindClusters(resolution=0.5)

## 注释（同 11b 逻辑：非免疫/内皮/周细胞 = 肿瘤）
mk<-list(Immune="PTPRC",Myeloid=c("CD68","C1QB","LYZ","AIF1","CD163"),Tcell=c("CD3D","CD3E","IL7R"),
         Bplasma=c("CD79A","MS4A1","MZB1"),Endothelial=c("PECAM1","VWF","CLDN5"),Mural=c("RGS5","ACTA2","PDGFRB"))
for(nm in names(mk)){g<-intersect(mk[[nm]],rownames(o));if(length(g))o<-AddModuleScore(o,list(g),name=paste0("sc_",nm))}
cl<-levels(Idents(o)); cm<-sapply(paste0("sc_",names(mk),"1"),function(c)tapply(o@meta.data[[c]],Idents(o),mean)); colnames(cm)<-names(mk)
ass<-sapply(cl,function(k){v<-cm[k,]
  if(v["Immune"]>0.1 & v["Immune"]>v["Endothelial"] & v["Immune"]>v["Mural"]) names(which.max(v[c("Myeloid","Tcell","Bplasma")]))
  else if(v["Endothelial"]>0.3 & v["Endothelial"]>v["Mural"]) "Endothelial"
  else if(v["Mural"]>0.3) "Mural" else "Meningioma"})
ct<-ass[as.character(Idents(o))];names(ct)<-colnames(o); o$celltype<-factor(ct)
cat("=== cell types ===\n"); print(table(o$celltype))

## 侵袭性程序投影
P<-readRDS("results/deg/GSE16581_program_validation.rds")
up<-intersect(P$up_sym,rownames(o));dn<-intersect(P$dn_sym,rownames(o))
o<-AddModuleScore(o,list(up),name="AggrUp"); o<-AddModuleScore(o,list(dn),name="AggrDn"); o$AggrScore<-o$AggrUp1-o$AggrDn1

## ⭐ 核心：肿瘤细胞内 AggrScore 随 grade
tum<-subset(o,celltype=="Meningioma")
tum$grade<-factor(tum$grade,levels=c("Normal","I","II","III"))
cat("\n=== ⭐ AggrScore in TUMOR cells by grade ===\n"); print(round(tapply(tum$AggrScore,tum$grade,mean),3))
# 统计：肿瘤细胞 AggrScore ~ grade (有序，排除 Normal)
td<-data.frame(score=tum$AggrScore,grade=tum$grade,gsm=tum$gsm)
tg<-td[td$grade%in%c("I","II","III"),]; tg$grade<-droplevels(tg$grade)
# sample-level mean（避免伪重复）
sm<-aggregate(score~gsm+grade,tg,mean)
cat("sample-level tumor AggrScore by grade:\n"); print(sm[order(sm$grade),])
ct_test<-cor.test(as.numeric(tg$grade),tg$score,method="spearman")
cat(sprintf("cell-level Spearman(grade,score) rho=%.3f p=%.2e\n",ct_test$estimate,ct_test$p.value))
sm_aov<-summary(aov(score~grade,sm)); cat("sample-level ANOVA p=",signif(sm_aov[[1]]$`Pr(>F)`[1],3),"\n")

write.csv(sm,"results/scrna/GSE206647_tumor_aggr_by_grade_sample.csv",row.names=FALSE)
## 轻量 per-cell 元数据 + UMAP（图可独立重生，截断不影响）
um<-Embeddings(o,"umap"); cmd<-data.frame(cell=colnames(o),celltype=o$celltype,grade=o$grade,gsm=o$gsm,
  AggrScore=o$AggrScore,UMAP1=um[,1],UMAP2=um[,2])
write.csv(cmd,"results/scrna/GSE206647_cellmeta.csv",row.names=FALSE)
## 图（在重 rds 之前先存，确保不被截断丢失）
p1<-DimPlot(o,group.by="celltype",label=TRUE,repel=TRUE)+ggtitle("GSE206647 cell types (163,897 cells)")+NoLegend()
p2<-VlnPlot(tum,"AggrScore",group.by="grade",pt.size=0)+geom_boxplot(width=0.12,outlier.shape=NA,fill="white",alpha=0.6)+
    ggtitle("Tumor-cell aggressiveness program\nby WHO grade")+NoLegend()+labs(x="WHO grade",y="Aggressiveness score")
smf<-sm; smf$grade<-factor(smf$grade,levels=c("I","II","III"))
p.sm<-sm_aov[[1]]$`Pr(>F)`[1]
p3<-ggplot(smf,aes(grade,score))+geom_boxplot(width=0.5,outlier.shape=NA,fill="grey90")+
    geom_jitter(width=0.12,size=2.2,color="#a14b3d")+theme_classic(base_size=13)+
    labs(x="WHO grade",y="Sample mean tumor AggrScore",
         title=sprintf("Sample-level (n=16)\nANOVA p=%.2g",p.sm))
ggsave("results/scrna/fig_GSE206647_grade_program.png",(p1|p2|p3)+plot_layout(widths=c(1.3,1,1)),width=17,height=5.2,dpi=150)
ggsave("results/scrna/fig_GSE206647_grade_program.pdf",(p1|p2|p3)+plot_layout(widths=c(1.3,1,1)),width=17,height=5.2)
cat("saved cellmeta + grade program fig + csv\n")
saveRDS(o,"results/scrna/GSE206647_processed.rds")
cat("\nALL-DONE saved GSE206647_processed.rds\n")
