# 31_pseudobulk_cellcycle.R — R3-B
# 回应独立审稿 R2 #13：单细胞 "tumour-intrinsic" 证据不足 —— 深度(nCount/nFeature)≠肿瘤纯度，
# 且 cell-level score 存在同患者 pseudoreplication。
# 修正：① per-patient tumour pseudobulk 程序评分 → sample-level grade 检验(n=患者)
#       ② cell-cycle(S/G2M)/核糖体/线粒体 作协变量的 cell-level mixed model(患者随机效应)
# 独立重算，不复用旧 per-cell ANOVA。
suppressMessages({library(Seurat); library(Matrix)})
HAS_LME <- requireNamespace("lme4", quietly=TRUE)
sink("results/deg/R3B_pseudobulk_cellcycle.txt", split=TRUE)
cat("=== R3-B per-patient pseudobulk + cell-cycle control ===\n")
cat("lme4 available:", HAS_LME, "\n")

obj <- readRDS("results/scrna/GSE206647_processed.rds")
cat("loaded RDS; cells=", ncol(obj), " genes=", nrow(obj), "\n")
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data
ctcol <- if("celltype" %in% colnames(md)) "celltype" else grep("celltype|cell_type|annotation", colnames(md), value=TRUE)[1]
gcol  <- if("grade" %in% colnames(md)) "grade" else grep("grade", colnames(md), value=TRUE)[1]
scol  <- if("gsm" %in% colnames(md)) "gsm" else grep("gsm|sample|orig", colnames(md), value=TRUE)[1]
cat("using cols: celltype=",ctcol," grade=",gcol," sample=",scol,"\n")

## 程序基因（symbol）
P <- readRDS("results/deg/GSE16581_program_validation.rds")
up <- intersect(P$up_sym, rownames(obj)); dn <- intersect(P$dn_sym, rownames(obj))
cat("program genes in scRNA: up",length(up)," dn",length(dn),"\n")

## cell-cycle 评分（Seurat 内置 S/G2M 基因）
obj <- CellCycleScoring(obj, s.features=cc.genes.updated.2019$s.genes,
                        g2m.features=cc.genes.updated.2019$g2m.genes, set.ident=FALSE)
## 核糖体/线粒体比例
obj[["pct_ribo"]] <- PercentageFeatureSet(obj, pattern="^RP[SL]")
if(!"percent.mt" %in% colnames(obj@meta.data)) obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^MT-")

## ---------- 肿瘤细胞子集（Meningioma） ----------
is_tumor <- md[[ctcol]] %in% c("Meningioma","Tumor","Tumour")
cat("tumour cells:", sum(is_tumor), "/", nrow(md), "\n")
tum <- subset(obj, cells = rownames(md)[is_tumor])

## ---------- ① per-patient pseudobulk（仅肿瘤细胞，原始 counts 求和→CPM→log）----------
cnt <- GetAssayData(tum, layer="counts")
samp <- tum@meta.data[[scol]]; grd <- tum@meta.data[[gcol]]
samp_grade <- tapply(as.character(grd), samp, function(x) names(sort(table(x),decreasing=TRUE))[1])
sl <- split(colnames(cnt), samp)
pb <- sapply(sl, function(cells) Matrix::rowSums(cnt[,cells,drop=FALSE]))
cpm <- t(t(pb)/colSums(pb))*1e6; logcpm <- log1p(cpm)
z <- t(scale(t(logcpm)))
score_pb <- colMeans(z[up,,drop=FALSE],na.rm=TRUE) - colMeans(z[dn,,drop=FALSE],na.rm=TRUE)
pbdf <- data.frame(gsm=names(score_pb), grade=samp_grade[names(score_pb)], score=score_pb,
                   ncells=sapply(sl,length)[names(score_pb)])
pbdf <- pbdf[!is.na(pbdf$grade),]
cat("\n[pseudobulk] per-patient tumour program score (n samples=",nrow(pbdf),"):\n")
print(pbdf[order(pbdf$grade),], row.names=FALSE)
write.csv(pbdf, "results/deg/R3B_pseudobulk_score_per_patient.csv", row.names=FALSE)
## 是否仍随 grade 递增（normal 排除后做 ordered test）
pt <- pbdf[pbdf$grade %in% c("I","II","III"),]
if(nrow(pt)>=4 && length(unique(pt$grade))>=2){
  pt$go <- as.integer(factor(pt$grade, levels=c("I","II","III")))
  ct <- suppressWarnings(cor.test(pt$go, pt$score, method="spearman"))
  cat(sprintf("\n[pseudobulk sample-level] Spearman score~grade-ordinal rho=%.3f p=%.4g (n=%d patients)\n",
      ct$estimate, ct$p.value, nrow(pt)))
  cat(sprintf("ANOVA(tumour pseudobulk score ~ grade): p=%.4g\n",
      summary(aov(score~factor(grade), data=pt))[[1]][["Pr(>F)"]][1]))
}

## ---------- ② cell-level mixed model：grade 信号在控制 cell-cycle/depth/ribo/mito 后是否仍在 ----------
cm <- tum@meta.data
cm$AggrScore <- tum@meta.data[["AggrScore"]]
if(is.null(cm$AggrScore)){
  tum <- AddModuleScore(tum, features=list(up), name="UP"); tum<-AddModuleScore(tum,features=list(dn),name="DN")
  cm$AggrScore <- tum$UP1 - tum$DN1
}
cm$grade_ord <- as.integer(factor(as.character(cm[[gcol]]), levels=c("Normal","I","II","III")))
cm$samp <- factor(cm[[scol]])
keep <- !is.na(cm$grade_ord)
cm <- cm[keep,]
cat("\n[cell-level] n cells=",nrow(cm)," patients=",nlevels(droplevels(cm$samp)),"\n")
if(HAS_LME){
  ## 朴素：AggrScore ~ grade_ord（cell-level，pseudoreplication）
  m0 <- lme4::lmer(AggrScore ~ grade_ord + (1|samp), data=cm)
  ## 控制 cell-cycle/depth/ribo/mito
  m1 <- lme4::lmer(AggrScore ~ grade_ord + S.Score + G2M.Score + nFeature_RNA + pct_ribo + percent.mt + (1|samp), data=cm)
  cf0 <- summary(m0)$coefficients["grade_ord",]
  cf1 <- summary(m1)$coefficients["grade_ord",]
  cat(sprintf("mixed model grade_ord beta (random intercept=patient):\n  unadjusted   = %.4f (t=%.2f)\n  cc/depth-adj = %.4f (t=%.2f)\n",
      cf0[1], cf0[3], cf1[1], cf1[3]))
  cat("→ 若调整后 beta 同号且 |t|>2，则 grade 关联非 cell-cycle/深度假象（但显著性须以 sample-level pseudobulk 为准，避免 pseudoreplication）。\n")
} else cat("lme4 缺失，跳过 mixed model；以 sample-level pseudobulk 为主结论。\n")

sink(); cat("\nALL-DONE -> results/deg/R3B_pseudobulk_cellcycle.txt\n")
