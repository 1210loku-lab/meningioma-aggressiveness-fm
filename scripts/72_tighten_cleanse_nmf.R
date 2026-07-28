# 72_tighten_cleanse_nmf.R
# One principled tightening pass: the first cleanse (script 71) left residual
# myeloid signal (a CD74/CCL3/CCL4/TYROBP metaprogram survived). Tighten the
# PREDEFINED, non-circular lineage filter using MYELOID/IMMUNE-SPECIFIC markers
# at >=1 (genes a bona fide meningioma tumour cell should not strongly express),
# while avoiding ACTA2/TAGLN (which meningioma mesenchymal cells do express).
# Re-run per-patient NMF metaprogram discovery ONLY (core results already shown
# robust in script 71). Decide: is a proliferation MP the clean grade-rising top?
# Run: Rscript scripts/72_tighten_cleanse_nmf.R

suppressMessages({library(Seurat); library(RcppML); library(Matrix); library(dplyr)})
set.seed(1)
mpdir <- "results/scrna/metaprograms_cleansed2"; dir.create(mpdir, showWarnings=FALSE, recursive=TRUE)
obj <- readRDS("results/scrna/GSE206647_processed.rds")
dat <- GetAssayData(obj, layer="data"); genes <- rownames(dat); md <- obj@meta.data
tum_idx <- which(md$celltype == "Meningioma")

# lineage-specific markers; threshold = >=N genes at lognorm>1
specs <- list(
  Immune_spec = list(g=c("PTPRC","TYROBP","FCER1G","C1QA","C1QB","C1QC","LYZ","AIF1",
                          "CD68","CSF1R","ITGAM","MS4A6A","LST1","CCL3","CCL4","SRGN","CD14"), n=1),
  Tcell = list(g=c("CD3D","CD3E","CD2","CD8A","NKG7","GZMB","GZMK","IL7R","CCL5"), n=1),
  Mast  = list(g=c("TPSAB1","TPSB2","CPA3","MS4A2"), n=1),
  Bplasma = list(g=c("CD79A","CD79B","MS4A1","MZB1","IGHG1"), n=1),
  Endoth = list(g=c("PECAM1","VWF","CLDN5","FLT1","EGFL7"), n=1),
  Mural  = list(g=c("MYH11","RGS5","PDGFRB","NOTCH3","KCNJ8"), n=1),   # NOT ACTA2/TAGLN
  Glial  = list(g=c("PLP1","MBP","MAG","MOG","OLIG1","OLIG2","PDGFRA","GFAP"), n=1)
)
contam <- rep(FALSE, length(tum_idx))
for (nm in names(specs)) {
  g <- intersect(specs[[nm]]$g, genes); if (!length(g)) next
  npos <- Matrix::colSums(dat[g, tum_idx, drop=FALSE] > 1)
  contam <- contam | (npos >= specs[[nm]]$n)
}
keep <- !contam
clean_cells <- colnames(obj)[tum_idx[keep]]
gsm_all <- as.character(md$gsm[tum_idx]); grade_all <- as.character(md$grade[tum_idx])
cat(sprintf("Tight cleanse: flagged %d (%.1f%%), retained %d clean tumour cells\n",
            sum(contam), 100*mean(contam), sum(keep)))
pre <- table(gsm_all); post <- table(gsm_all[keep])
patc <- data.frame(gsm=names(pre), grade=grade_all[match(names(pre),gsm_all)],
                   n_pre=as.integer(pre), n_post=as.integer(post[names(pre)]))
patc$n_post[is.na(patc$n_post)] <- 0L
patc$removed_pct <- round(100*(1-patc$n_post/patc$n_pre),1)
print(patc[order(patc$grade,-patc$removed_pct),], row.names=FALSE)

# quick core re-check on tighter set: AggrScore ~ grade
aggr <- md$AggrScore[tum_idx]
newp <- tapply(aggr[keep], gsm_all[keep], mean)
gm <- setNames(grade_all[match(names(newp),gsm_all)], names(newp))
ok <- gm %in% c("I","II","III")
ord <- setNames(c(1,2,3),c("I","II","III"))[gm[ok]]
ct <- suppressWarnings(cor.test(as.numeric(newp[ok]), ord, method="spearman", exact=FALSE))
cat(sprintf("AggrScore~grade (tight cleanse): rho=%.3f p=%.3g\n", ct$estimate, ct$p.value))

MIN_CELLS <- 100L
patients <- patc$gsm[patc$n_post >= MIN_CELLS & patc$grade %in% c("I","II","III")]
cat("patients kept for NMF:", length(patients), "\n")

## NMF metaprograms (script-70 logic)
TOPG <- 50L; KS <- 4:9; prog_list <- list(); prog_meta <- list()
for (p in patients) {
  cells <- clean_cells[as.character(md[clean_cells,"gsm"])==p]
  if (length(cells) < 200) next
  sub <- obj[, cells]
  sub <- FindVariableFeatures(sub, selection.method="vst", nfeatures=2000, verbose=FALSE)
  X <- as.matrix(GetAssayData(sub, layer="data")[VariableFeatures(sub), ]); X <- X-rowMeans(X); X[X<0]<-0
  X <- X[rowSums(X)>0,,drop=FALSE]
  for (k in KS) {
    fit <- tryCatch(RcppML::nmf(X,k=k,verbose=FALSE,seed=1), error=function(e) NULL); if (is.null(fit)) next
    W <- fit$w; rownames(W) <- rownames(X)
    for (j in seq_len(ncol(W))) {
      id <- sprintf("%s_k%d_f%d",p,k,j); prog_list[[id]] <- names(sort(W[,j],decreasing=TRUE))[1:TOPG]
      prog_meta[[id]] <- data.frame(id=id, patient=p, grade=patc$grade[patc$gsm==p], k=k, factor=j)
    }
  }
}
prog_meta <- bind_rows(prog_meta)
ids <- names(prog_list); n <- length(ids); J <- matrix(0,n,n,dimnames=list(ids,ids))
for (i in 1:(n-1)) for (j in (i+1):n){jac<-length(intersect(prog_list[[i]],prog_list[[j]]))/length(union(prog_list[[i]],prog_list[[j]]));J[i,j]<-jac;J[j,i]<-jac}
cl <- cutree(hclust(as.dist(1-J),method="average"), h=0.9); prog_meta$cluster <- cl[prog_meta$id]
cl_tab <- prog_meta %>% group_by(cluster) %>% summarise(n_prog=n(), n_pat=n_distinct(patient), .groups="drop") %>%
  filter(n_pat>=3) %>% arrange(desc(n_pat))
mp_sig <- list()
for (cc in cl_tab$cluster){members<-prog_meta$id[prog_meta$cluster==cc];tg<-table(unlist(prog_list[members]));thr<-max(2,ceiling(0.25*length(members)));s<-names(sort(tg[tg>=thr],decreasing=TRUE));if(length(s)<10)s<-names(sort(tg,decreasing=TRUE))[1:min(30,length(tg))];mp_sig[[as.character(cc)]]<-head(s,50)}
names(mp_sig) <- paste0("MP", seq_along(mp_sig))
tsc <- AddModuleScore(obj[,clean_cells], features=mp_sig, name="MP", seed=1)
sc_cols <- paste0("MP", seq_along(mp_sig)); mdc <- tsc@meta.data
aggr_cor <- sapply(sc_cols, function(s) cor(mdc[[s]], mdc$AggrScore, method="spearman"))
pa <- mdc %>% group_by(gsm) %>% summarise(across(all_of(sc_cols),mean), grade=grade[1], .groups="drop")
pa <- pa[pa$grade %in% c("I","II","III"),]; pa$go <- setNames(c(1,2,3),c("I","II","III"))[pa$grade]
grho <- sapply(sc_cols, function(s) suppressWarnings(cor(pa[[s]],pa$go,method="spearman")))
gp <- sapply(sc_cols, function(s) suppressWarnings(cor.test(pa[[s]],pa$go,method="spearman",exact=FALSE)$p.value))
summ <- data.frame(metaprogram=sc_cols, n_patients=cl_tab$n_pat, aggr_spearman=round(aggr_cor,3),
                   grade_rho=round(grho,3), grade_p=signif(gp,3)) %>% arrange(desc(aggr_spearman))
sig_df <- do.call(rbind, lapply(names(mp_sig), function(m) data.frame(metaprogram=m, rank=seq_along(mp_sig[[m]]), gene=mp_sig[[m]])))
write.csv(summ, file.path(mpdir,"metaprogram_summary.csv"), row.names=FALSE)
write.csv(sig_df, file.path(mpdir,"metaprogram_signatures.csv"), row.names=FALSE)
saveRDS(list(clean_cells=clean_cells, patc=patc), file.path(mpdir,"clean_cells2.rds"))
cc_genes <- c("MKI67","TOP2A","CDK1","TPX2","BIRC5","CENPF","UBE2C","NUSAP1","PTTG1","ASPM","STMN1","TYMS")
cat("\n=== TIGHT-CLEANSE metaprogram summary ===\n")
for (m in summ$metaprogram) {
  g <- head(sig_df$gene[sig_df$metaprogram==m], 20); r <- summ[summ$metaprogram==m,]
  cat(sprintf("%-5s npat=%2d aggr=%+.2f gradeRho=%+.2f p=%.3g cc=%d\n  %s\n", m, r$n_patients, r$aggr_spearman, r$grade_rho, r$grade_p, sum(g %in% cc_genes), paste(g,collapse=", ")))
}
cat("\nDONE (metaprograms_cleansed2)\n")
