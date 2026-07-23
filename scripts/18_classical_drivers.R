# 18_classical_drivers.R — 经典可解释驱动归因（基础模型对照基线）
# 对侵袭性程序基因排序：RandomForest 重要性 + DESeq2 |stat| + LASSO 选择，输出候选驱动表
suppressMessages({library(randomForest); library(DESeq2)})
set.seed(42); outdir<-"results/deg"

a<-readRDS("data/raw/GSE136661_assembled.rds")
prog<-readRDS("results/deg/GSE136661_aggressiveness_program.rds")
genes<-c(prog$up,prog$dn); genes<-intersect(genes,rownames(a$counts))
cat("program genes in counts:",length(genes),"\n")

## 表达：vst（仅程序基因足够稳健）
cts<-a$counts[, ]; storage.mode(cts)<-"integer"
grade<-a$pd[["pathology:ch1"]]
keep<-!is.na(grade) & grade %in% c("WHO I","WHO II","WHO III")
cts<-cts[,keep]; grade<-grade[keep]
y<-factor(ifelse(grade=="WHO I","low","high"),levels=c("low","high"))  # 与侵袭性对比一致
cat("samples low/high:",table(y),"\n")
vsd<-varianceStabilizingTransformation(cts[rowSums(cts>=5)>=10,],blind=TRUE)
vmat<-vsd[intersect(genes,rownames(vsd)),]
cat("vst genes:",nrow(vmat),"\n")

## RandomForest 重要性（grade 高/低）
df<-data.frame(t(vmat)); df$y<-y
rf<-randomForest(y~.,data=df,ntree=2000,importance=TRUE)
imp<-importance(rf)[,"MeanDecreaseGini"]
oob<-rf$err.rate[nrow(rf$err.rate),"OOB"]
cat("RF OOB error:",round(oob,3),"\n")

## DESeq2 |stat|
res<-read.csv("results/deg/GSE136661_aggressiveness_WHO_II_vs_I.csv",row.names=1)
stat<-abs(res[names(imp),"stat"])

## ENSG -> symbol
sym<-names(imp)
if(requireNamespace("org.Hs.eg.db",quietly=TRUE)){
  m<-AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,keys=names(imp),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
  sym<-ifelse(is.na(m),names(imp),m)
}

drv<-data.frame(ensembl=names(imp),symbol=sym,
                direction=ifelse(names(imp)%in%prog$up,"up","down"),
                RF_importance=round(imp,2),
                DESeq2_absStat=round(stat,2))
# 共识秩（两指标秩平均）
drv$rank_RF<-rank(-drv$RF_importance); drv$rank_stat<-rank(-drv$DESeq2_absStat)
drv$consensus_rank<-rowMeans(drv[,c("rank_RF","rank_stat")])
drv<-drv[order(drv$consensus_rank),]
write.csv(drv,file.path(outdir,"classical_drivers_ranked.csv"),row.names=FALSE)
cat("\n=== top 15 candidate drivers (consensus) ===\n"); print(head(drv[,c("symbol","direction","RF_importance","DESeq2_absStat","consensus_rank")],15))

## 图：top20 驱动 RF 重要性
suppressMessages({library(ggplot2); library(patchwork)})
top<-head(drv,20); top$symbol<-factor(top$symbol,levels=rev(top$symbol))
pA<-ggplot(top,aes(symbol,RF_importance,fill=direction))+geom_col()+coord_flip()+
  scale_fill_manual(values=c(up="#a14b3d",down="#3b6ea5"))+theme_classic(base_size=12)+
  labs(x=NULL,y="RandomForest importance (MeanDecreaseGini)",
       title=sprintf("A. Top candidate genes by RF importance\nOOB err=%.2f, GSE136661",oob))

top_label <- head(drv,12)
pB <- ggplot(drv,aes(DESeq2_absStat,RF_importance,color=direction))+
  geom_point(alpha=0.65,size=1.8)+
  geom_text(data=top_label,aes(label=symbol),size=3,check_overlap=TRUE,vjust=-0.45,show.legend=FALSE)+
  scale_color_manual(values=c(up="#a14b3d",down="#3b6ea5"))+
  theme_classic(base_size=12)+
  labs(x="DESeq2 |Wald statistic|",y="RandomForest importance",
       color="Program direction",
       title="B. DE strength vs RF importance")

p <- pA | pB
ggsave("results/figures_pub/fig_classical_drivers.png",p,width=13,height=6,dpi=150)
ggsave("results/figures_pub/fig_classical_drivers.pdf",p,width=13,height=6)
saveRDS(list(rf=rf,drv=drv,oob=oob),file.path(outdir,"classical_drivers.rds"))
cat("\nALL-DONE saved classical_drivers_ranked.csv + fig_classical_drivers + rds\n")
