# Submission-stage independent headline audit.
#
# This script writes only to results/audit_submission/. It deliberately does not
# source analysis scripts or read their headline result tables/RDS objects. Bulk
# statistics are recomputed from assembled/raw GEO inputs. The GSE206647 audit
# re-aggregates raw 10x counts using the frozen cell identity labels in the saved
# cell metadata; it does not reuse the old pseudobulk scores or target matrix.
#
# ============================================================================
# PARTIALLY SUPERSEDED (2026-07-17 review-response round). DO NOT backfill the
# two fields below into the manuscript — they reflect pre-revision methodology:
#
#   * gse74385_firth_or_per_sd (== 3.36 here) uses the OLD recurrence-OR-malignant-
#     -progression COMPOSITE endpoint, grade-adjusted only (no batch). The current
#     PRIMARY endpoint is recurrence-only, grade+batch-adjusted -> OR 5.18
#     (95% CI 1.40-25.83, p=0.0118). Canonical source: scripts/53_GSE74385_endpoint_sensitivity.R
#     -> results/audit_submission/GSE74385_endpoint_sensitivity.csv
#
#   * geneformer_best_pc_abs_rho_grade / _selection_corrected_p (== 0.539 / 5e-4 here)
#     are the OLD CELL-LEVEL PCA with cell-label permutation (n=2000), which is
#     pseudoreplicated (grade is a 16-patient label). The current analysis is
#     patient-level: primary p=0.0516, sensitivity p=0.00070, 100,000 patient-label
#     permutations. Canonical source: scripts/51_repair_geneformer_patient_mapping.R
#     + scripts/52_geneformer_patient_level_audit.R
#     -> results/audit_submission/geneformer_patient_level_axis_summary.csv
#
# STILL CURRENT and manuscript-consistent: gse16581_firth_* (OR 14.74), grade/
# recurrence p-values, discovery DEG/program membership, gse206647_pseudobulk_*,
# measurable/aligned gene counts, PI3/PITX1/NF2/LTK rho+FDR.
# See results/audit_submission/independent_headline_audit_SUPERSEDED_NOTE.md.
# ============================================================================

suppressMessages({
  library(DESeq2)
  library(GEOquery)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(illuminaHumanv4.db)
  library(Matrix)
  library(logistf)
  library(glmnet)
  library(pROC)
})

set.seed(20260715)
root <- normalizePath(".")
outdir <- file.path(root, "results", "audit_submission")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
headline <- list()

## 1. Recompute the fixed programme from GSE136661 counts.
a <- readRDS(file.path(root, "data", "raw", "GSE136661_assembled.rds"))
counts <- a$counts; key <- a$key
rownames(key) <- key$gsm; key <- key[colnames(counts),,drop=FALSE]
grade <- factor(gsub("WHO ", "", key[["pathology:ch1"]]), levels=c("I","II","III"))
cohort <- key[["cohort:ch1"]]
primary <- cohort %in% c("Discovery","Validation") & grade %in% c("I","II")
mp <- counts[,primary,drop=FALSE]
kp <- data.frame(grade=droplevels(grade[primary]), row.names=colnames(mp))
storage.mode(mp) <- "integer"
mp <- mp[rowSums(mp >= 5) >= 10,,drop=FALSE]
dds <- DESeqDataSetFromMatrix(mp, kp, ~grade)
dds <- DESeq(dds, quiet=TRUE)
de <- as.data.frame(results(dds, contrast=c("grade","II","I")))
de <- de[order(de$padj),]
up <- rownames(head(de[!is.na(de$padj) & de$padj<0.05 & de$log2FoldChange>1,],100))
dn <- rownames(head(de[!is.na(de$padj) & de$padj<0.05 & de$log2FoldChange< -1,],100))
headline$discovery_degs_fdr05 <- sum(de$padj<0.05,na.rm=TRUE)
headline$program_up <- length(up); headline$program_down <- length(dn)
program_saved <- read.csv(file.path(root,"docs","Table_S2_aggressiveness_program_genes.csv"),stringsAsFactors=FALSE)
headline$program_membership_match <- setequal(c(up,dn), program_saved$ensembl)

up.sym <- program_saved$symbol[program_saved$direction=="up_in_WHO_II_vs_I"]
dn.sym <- program_saved$symbol[program_saved$direction=="down_in_WHO_II_vs_I"]

## 2. GSE16581 from local series matrix + local GPL annotation.
g165 <- getGEO(filename=file.path(root,"data","raw","GSE16581_series_matrix.txt.gz"),getGPL=FALSE)
x165 <- exprs(g165); p165 <- pData(g165)
if(max(x165,na.rm=TRUE)>100) x165 <- log2(x165+1)
gpl <- getGEO(filename=file.path(root,"data","raw","GPL570.soft.gz"))
tab <- Table(gpl)
symcol <- grep("gene symbol|symbol",colnames(tab),ignore.case=TRUE,value=TRUE)[1]
probe.map <- setNames(as.character(tab[[symcol]]),as.character(tab$ID))
sym <- probe.map[rownames(x165)]
valid165 <- !is.na(sym) & sym!=""
x165 <- x165[valid165,,drop=FALSE]; sym <- sym[valid165]
ord165 <- order(rowMeans(x165),decreasing=TRUE)
x165 <- x165[ord165,,drop=FALSE]; sym <- sym[ord165]
dedup165 <- !duplicated(sym); x165 <- x165[dedup165,,drop=FALSE]; rownames(x165) <- sym[dedup165]
z165 <- t(scale(t(x165)))
score165 <- colMeans(z165[intersect(up.sym,rownames(z165)),,drop=FALSE],na.rm=TRUE)-
            colMeans(z165[intersect(dn.sym,rownames(z165)),,drop=FALSE],na.rm=TRUE)
g165.grade <- factor(p165[["who grade:ch1"]])
rf <- suppressWarnings(as.integer(p165[["recurrence_frequency:ch1"]]))
rec165 <- as.integer(rf>0)
zscore165 <- as.numeric(scale(score165))
fit165 <- logistf(rec165 ~ zscore165 + g165.grade)
headline$gse16581_n <- length(rec165); headline$gse16581_events <- sum(rec165)
headline$gse16581_firth_or_per_sd <- exp(coef(fit165)["zscore165"])
headline$gse16581_firth_ci_low <- exp(fit165$ci.lower["zscore165"])
headline$gse16581_firth_ci_high <- exp(fit165$ci.upper["zscore165"])
headline$gse16581_firth_p <- fit165$prob["zscore165"]
headline$gse16581_grade_anova_p <- anova(lm(score165~g165.grade))$`Pr(>F)`[1]
headline$gse16581_recurrence_t_p <- t.test(score165~factor(rec165))$p.value

## 3. GSE74385 from raw normalized matrix and phenotype.
g743 <- getGEO(filename=file.path(root,"data","raw","GSE74385","GSE74385_series_matrix.txt.gz"),getGPL=FALSE)
p743 <- pData(g743)
subtype <- gsub("subtype: ","",p743[,grep("subtype",colnames(p743),ignore.case=TRUE)[1]])
meta743 <- data.frame(title=as.character(p743$title),subtype=subtype,
                      grade=substr(subtype,1,1),outcome=gsub("[0-9]","",subtype),
                      stringsAsFactors=FALSE)
meta743$batch <- factor(gsub("batch: ","",p743[,grep("batch",colnames(p743),ignore.case=TRUE)[1]]))
raw743 <- read.delim(gzfile(file.path(root,"data","raw","GSE74385","GSE74385_normalized.txt.gz")),
                     check.names=FALSE,stringsAsFactors=FALSE)
rownames(raw743) <- raw743[,1]; raw743 <- raw743[,-1]
raw743 <- raw743[,!grepl("Detection",colnames(raw743)),drop=FALSE]
raw743 <- raw743[,meta743$title,drop=FALSE]
map743 <- mapIds(illuminaHumanv4.db,rownames(raw743),"SYMBOL","PROBEID",multiVals="first")
valid <- !is.na(map743) & map743!=""
means <- rowMeans(raw743[valid,,drop=FALSE])
ord <- order(means,decreasing=TRUE)
idx <- which(valid)[ord]
idx <- idx[!duplicated(as.character(map743[idx]))]
x743 <- as.matrix(raw743[idx,,drop=FALSE]); rownames(x743) <- as.character(map743[idx])
z743 <- t(scale(t(x743)))
up743 <- intersect(up.sym,rownames(z743)); dn743 <- intersect(dn.sym,rownames(z743))
score743 <- colMeans(z743[up743,,drop=FALSE],na.rm=TRUE)-colMeans(z743[dn743,,drop=FALSE],na.rm=TRUE)
meta743$score <- score743[meta743$title]
# SUPERSEDED endpoint: R and M pooled = composite (grade-adjusted only) -> OR 3.36.
# Current primary is recurrence-only + grade+batch (OR 5.18); see scripts/53. Header note above.
meta743$recbin <- ifelse(meta743$outcome=="NR",0,ifelse(meta743$outcome %in% c("R","M"),1,NA))
known <- meta743[!is.na(meta743$recbin),]
headline$gse74385_total_n <- nrow(meta743); headline$gse74385_known_n <- nrow(known)
headline$gse74385_unknown_n <- sum(is.na(meta743$recbin))
headline$gse74385_wilcoxon_p <- wilcox.test(score~factor(recbin),known)$p.value
fit.batch <- glm(recbin~score+batch,known,family=binomial)
headline$gse74385_batch_adjusted_score_p <- summary(fit.batch)$coef["score","Pr(>|z|)"]
fit743 <- logistf(recbin~scale(score)+factor(grade),data=known)
coef.name <- grep("scale",names(coef(fit743)),value=TRUE)[1]
headline$gse74385_firth_or_per_sd <- exp(coef(fit743)[coef.name])
headline$gse74385_firth_ci_low <- exp(fit743$ci.lower[coef.name])
headline$gse74385_firth_ci_high <- exp(fit743$ci.upper[coef.name])
headline$gse74385_firth_p <- fit743$prob[coef.name]
write.csv(meta743,file.path(outdir,"independent_GSE74385_scores.csv"),row.names=FALSE)

## 4. Re-aggregate raw GSE206647 10x counts using frozen cell identities.
cellmeta <- read.csv(file.path(root,"results","scrna","GSE206647_cellmeta.csv"),stringsAsFactors=FALSE)
program.sym <- unique(c(up.sym,dn.sym))
samples <- sort(unique(cellmeta$gsm))
pb.tum <- matrix(0,nrow=length(program.sym),ncol=length(samples),dimnames=list(program.sym,samples))
lib.tum <- setNames(rep(0,length(samples)),samples)
total.tum <- setNames(rep(0,length(program.sym)),program.sym)
total.other <- total.tum; n.tum <- 0; n.other <- 0
for(gsm in samples){
  prefix <- list.files(file.path(root,"data","raw","GSE206647_ex"),pattern=paste0("^",gsm,"_.*_matrix\\.mtx\\.gz$"),full.names=TRUE)
  stopifnot(length(prefix)==1)
  base <- sub("_matrix\\.mtx\\.gz$","",prefix)
  feat <- read.delim(gzfile(paste0(base,"_features.tsv.gz")),header=FALSE,stringsAsFactors=FALSE)
  bar <- readLines(gzfile(paste0(base,"_barcodes.tsv.gz")))
  mat.all <- readMM(gzfile(prefix))
  fullcell <- paste0(gsm,"_",bar)
  md <- cellmeta[cellmeta$gsm==gsm,]
  cols <- match(md$cell,fullcell); stopifnot(!anyNA(cols))
  sym.raw <- as.character(feat[[2]])
  rows <- which(sym.raw %in% program.sym)
  is.tum <- md$celltype=="Meningioma"
  lib.tum[gsm] <- sum(mat.all[,cols[is.tum],drop=FALSE])
  mat <- mat.all[rows,cols,drop=FALSE]; sym.raw <- sym.raw[rows]
  collapse <- function(cols.use){
    if(!length(cols.use)) return(setNames(rep(0,length(program.sym)),program.sym))
    v <- Matrix::rowSums(mat[,cols.use,drop=FALSE])
    vv <- rowsum(v,sym.raw,reorder=FALSE)[,1]
    out <- setNames(rep(0,length(program.sym)),program.sym); out[names(vv)] <- vv; out
  }
  vt <- collapse(which(is.tum)); vo <- collapse(which(!is.tum))
  pb.tum[,gsm] <- vt; total.tum <- total.tum+vt; total.other <- total.other+vo
  n.tum <- n.tum+sum(is.tum); n.other <- n.other+sum(!is.tum)
  rm(mat,mat.all); gc(verbose=FALSE)
}
logcpm <- log1p(t(t(pb.tum)/lib.tum)*1e6)
zpb <- t(scale(t(logcpm)))
pb.score <- colMeans(zpb[intersect(up.sym,rownames(zpb)),,drop=FALSE],na.rm=TRUE)-
            colMeans(zpb[intersect(dn.sym,rownames(zpb)),,drop=FALSE],na.rm=TRUE)
sample.grade <- sapply(samples,function(g) names(sort(table(cellmeta$grade[cellmeta$gsm==g]),decreasing=TRUE))[1])
tumour <- sample.grade %in% c("I","II","III")
grade.ord <- setNames(c(1,2,3),c("I","II","III"))[sample.grade[tumour]]
ct <- suppressWarnings(cor.test(pb.score[tumour],grade.ord,method="spearman",exact=FALSE))
headline$gse206647_tumour_patients <- sum(tumour)
headline$gse206647_pseudobulk_rho <- unname(ct$estimate)
headline$gse206647_pseudobulk_p <- ct$p.value
headline$gse206647_pseudobulk_anova_p <- summary(aov(pb.score[tumour]~factor(sample.grade[tumour])))[[1]][["Pr(>F)"]][1]

rho <- pval <- setNames(rep(NA_real_,length(program.sym)),program.sym)
for(g in program.sym){
  x <- logcpm[g,tumour]
  if(length(unique(x))>1){
    tt <- suppressWarnings(cor.test(x,grade.ord,method="spearman",exact=FALSE))
    rho[g] <- unname(tt$estimate); pval[g] <- tt$p.value
  }
}
fdr <- p.adjust(pval,"BH")
direction <- setNames(c(rep("up",length(up.sym)),rep("down",length(dn.sym))),c(up.sym,dn.sym))
aligned <- ifelse(direction[names(rho)]=="up",rho>0,rho<0)
headline$gse206647_measurable_program_genes <- sum(!is.na(rho))
headline$gse206647_aligned_fdr10 <- sum(aligned & fdr<0.10,na.rm=TRUE)
specificity <- log2((total.tum/n.tum+0.01)/(total.other/n.other+0.01))
target <- data.frame(symbol=names(rho),direction=direction[names(rho)],rho=rho,p=pval,fdr=fdr,
                     aligned=aligned,tumour_specificity_log2=specificity[names(rho)])
write.csv(target,file.path(outdir,"independent_GSE206647_target_evidence.csv"),row.names=FALSE)
write.csv(data.frame(gsm=samples,grade=sample.grade,score=pb.score),
          file.path(outdir,"independent_GSE206647_pseudobulk_scores.csv"),row.names=FALSE)
for(g in c("PI3","PITX1","NF2","LTK")){
  headline[[paste0(tolower(g),"_rho")]] <- rho[g]
  headline[[paste0(tolower(g),"_fdr")]] <- fdr[g]
}

## 5. Independent Geneformer PC-selection permutation from saved embeddings.
## SUPERSEDED: cell-level PCA + cell-label permutation (pseudoreplicated).
## Current analysis is patient-level (scripts/51 + 52); see header note above.
emb <- read.csv(file.path(root,"results","scrna","gf_out","tumor_emb.csv"),check.names=FALSE)
ecols <- grep("^[0-9]+$",colnames(emb),value=TRUE)
xx <- as.matrix(emb[,ecols])
pc <- prcomp(xx,center=TRUE,scale.=FALSE,rank.=10)$x
go <- setNames(c(1,2,3),c("I","II","III"))[as.character(emb$grade)]
rhos <- sapply(1:10,function(i) abs(cor(pc[,i],go,method="spearman")))
best <- which.max(rhos); obs <- rhos[best]
same.axis <- cor(pc[,best],emb$AggrScore,method="spearman")
set.seed(20260715)
null <- replicate(2000,max(sapply(1:10,function(i) abs(cor(pc[,i],sample(go),method="spearman")))))
headline$geneformer_best_pc_abs_rho_grade <- obs
headline$geneformer_same_axis_abs_rho_program <- abs(same.axis)
headline$geneformer_selection_corrected_p <- (sum(null>=obs)+1)/(length(null)+1)

## Write machine-readable and human-readable summaries.
values <- data.frame(metric=names(headline),value=vapply(headline,function(x) paste(x,collapse=";"),character(1)))
write.csv(values,file.path(outdir,"independent_headline_audit.csv"),row.names=FALSE)
report <- c(
  "Independent headline audit (raw/assembled inputs; old headline outputs not read)",
  paste(values$metric,values$value,sep=": ")
)
writeLines(report,file.path(outdir,"independent_headline_audit.txt"))
cat(paste(report,collapse="\n"),"\n")
