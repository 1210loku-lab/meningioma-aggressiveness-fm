# 16b_fig_GSE74385.R — GSE74385 复发外验图（程序评分 × outcome 全级别 + grade-I 特异）
suppressMessages({library(ggplot2); library(patchwork)})
source("scripts/_fig_style.R")
m<-read.csv("results/deg/GSE74385_program_recurrence.csv")
m$outcome<-factor(m$outcome,levels=c("NR","M","R"))
mk<-m[!is.na(m$outcome),]
p_known<-wilcox.test(score~recur,m)$p.value
# panel A: all grades by outcome
pA<-ggplot(mk,aes(outcome,score,fill=outcome))+geom_boxplot(width=0.6,outlier.shape=NA)+
  geom_jitter(width=0.12,size=1.6,alpha=0.7)+
  scale_fill_manual(values=c(NR="#3b6ea5",M="#caa54b",R="#a14b3d"))+
  theme_classic(base_size=12)+theme(legend.position="none")+
  labs(x="Outcome (NR / Malignant-prog / Recurrent)",y="Aggressiveness score",
       title=sprintf("GSE74385 known outcomes (n=%d)\nR/M vs NR p=%.1e",nrow(mk),p_known))
# panel B: grade-I — 诚实揭示 batch 混杂（1R 全 batch2，1NR/1M 全 batch1）
suppressMessages(library(GEOquery))
gg<-getGEO(filename="data/raw/GSE74385/GSE74385_series_matrix.txt.gz",getGPL=FALSE); ph<-pData(gg)
bt<-gsub("batch: ","",ph[,grep("batch",colnames(ph))[1]]); names(bt)<-ph$title
g1<-m[m$grade==1,]; g1$recur<-factor(g1$recur,levels=c("NonRecur","Recur/Malig"))
g1$batch<-factor(bt[g1$title])
pB<-ggplot(g1,aes(recur,score))+geom_boxplot(width=0.55,outlier.shape=NA,fill="grey92")+
  geom_jitter(aes(color=batch),width=0.12,size=2.4,alpha=0.9)+
  scale_color_manual(values=c("1"="#3b6ea5","2"="#d98a00"))+
  theme_classic(base_size=12)+
  labs(x="WHO grade I only",y="Aggressiveness score",color="Batch",
       title="Grade-I subset: CONFOUNDED by batch\n(all recurrent=batch2) — hypothesis-generating only")
figS1 <- (pA|pB)+plot_annotation(tag_levels="A") & fig_style
save_fig(figS1, "results/figures_pub/fig_GSE74385_recurrence", width=11, height=4.6)
cat("saved fig_GSE74385_recurrence.{png,pdf}\n")
