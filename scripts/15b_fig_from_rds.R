# 15b_fig_from_rds.R — 从已存 GSE206647_processed.rds 补 grade 程序图（脚本15漏 patchwork）
suppressMessages({library(Seurat); library(ggplot2); library(patchwork); library(dplyr)})
o<-readRDS("results/scrna/GSE206647_processed.rds")
tum<-subset(o,celltype=="Meningioma"); tum$grade<-factor(tum$grade,levels=c("Normal","I","II","III"))
# sample-level 均值（无伪重复）用于叠点
sm<-read.csv("results/scrna/GSE206647_tumor_aggr_by_grade_sample.csv")
sm$grade<-factor(sm$grade,levels=c("I","II","III"))

p1<-DimPlot(o,group.by="celltype",label=TRUE,repel=TRUE)+ggtitle("GSE206647 cell types (163,897 cells)")+NoLegend()
p2<-VlnPlot(tum,"AggrScore",group.by="grade",pt.size=0)+
    geom_boxplot(width=0.12,outlier.shape=NA,fill="white",alpha=0.6)+
    ggtitle("Tumor-cell aggressiveness program\nby WHO grade")+NoLegend()+
    labs(x="WHO grade",y="Aggressiveness score")
# sample-level dotplot (primary, non-pseudoreplicated)
p3<-ggplot(sm,aes(grade,score))+
    geom_boxplot(width=0.5,outlier.shape=NA,fill="grey90")+
    geom_jitter(width=0.12,size=2.2,color="#a14b3d")+
    theme_classic(base_size=13)+
    labs(x="WHO grade",y="Sample mean tumor AggrScore",
         title="Sample-level (n=16)\nANOVA p=2.3e-4")
ggsave("results/scrna/fig_GSE206647_grade_program.png",(p1|p2|p3)+plot_layout(widths=c(1.3,1,1)),width=17,height=5.2,dpi=150)
ggsave("results/scrna/fig_GSE206647_grade_program.pdf",(p1|p2|p3)+plot_layout(widths=c(1.3,1,1)),width=17,height=5.2)
cat("saved fig_GSE206647_grade_program.{png,pdf}\n")
