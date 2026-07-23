# 26_purity_control.R — P1 反驳"单细胞 grade 关联只是肿瘤纯度/测序深度":
# 检验肿瘤细胞内 AggrScore 是否被 nCount 驱动；残差化 nCount 后 grade 关联是否存活
suppressMessages({library(Seurat)})
sink("results/deg/P1_purity_control.txt", split=TRUE)
o <- readRDS("results/scrna/GSE206647_processed.rds")
o <- subset(o, celltype=="Meningioma" & grade %in% c("I","II","III"))
df <- data.frame(score=o$AggrScore, nCount=o$nCount_RNA, nFeat=o$nFeature_RNA,
                 grade=factor(o$grade,levels=c("I","II","III")), gsm=o$gsm)
cat("tumor cells:",nrow(df),"\n")

cat("\n--- 每细胞 AggrScore vs 测序深度 ---\n")
cat("Spearman(score, nCount) =", round(cor(df$score, df$nCount, method="spearman"),3),"\n")
cat("Spearman(score, nFeature)=", round(cor(df$score, df$nFeat, method="spearman"),3),"\n")

cat("\n--- 残差化：score ~ log10(nCount) 取残差，再看 grade 关联 ---\n")
df$resid <- residuals(lm(score ~ log10(nCount) + log10(nFeat), data=df))
# sample-level（无伪重复）
sm0 <- aggregate(score~gsm+grade, df, mean)
smr <- aggregate(resid~gsm+grade, df, mean)
cat("原始 sample-level ANOVA p =", signif(summary(aov(score~grade, sm0))[[1]]$`Pr(>F)`[1],3),"\n")
cat("纯度校正后 sample-level ANOVA p =", signif(summary(aov(resid~grade, smr))[[1]]$`Pr(>F)`[1],3),"\n")
cat("\n原始 sample 均值:\n"); print(round(tapply(sm0$score,sm0$grade,mean),3))
cat("校正后 sample 均值:\n"); print(round(tapply(smr$resid,smr$grade,mean),3))
cat("\n解读：若校正深度后 grade 的 ANOVA 仍显著，则 grade 关联非纯度/深度假象。\n")
sink(); cat("\nALL-DONE -> results/deg/P1_purity_control.txt\n")
