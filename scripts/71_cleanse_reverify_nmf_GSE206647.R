# 71_cleanse_reverify_nmf_GSE206647.R
# Plan C (post SR desk-reject): stringently cleanse the GSE206647 "Meningioma"
# annotation using PREDEFINED lineage markers (independent of AggrScore/NMF, i.e.
# non-circular), then INDEPENDENTLY re-verify every single-cell main result on the
# cleansed tumour set and re-run the per-patient NMF metaprogram discovery.
# Reports OLD vs NEW transparently; the new numbers replace the old ones.
#
# Inputs : results/scrna/GSE206647_processed.rds (celltype, grade, gsm, AggrScore)
#          docs/Table_S2_aggressiveness_program_genes.csv (200-gene program)
# Outputs: results/scrna/cleansed/
#   cleanse_per_patient.csv        pre/post counts, removed fraction, dropped flag
#   reverify_aggrscore_grade.csv   OLD vs NEW per-patient AggrScore + Spearman/ANOVA
#   reverify_gene_evidence.csv     script-36 gene-evidence recomputed on cleansed
#   mixedmodel_slopes.txt          cell-level lme4 slopes (unadj/adj)
#   metaprograms_cleansed/*.csv    re-run NMF metaprograms on cleansed cells
# Run: Rscript scripts/71_cleanse_reverify_nmf_GSE206647.R

suppressMessages({library(Seurat); library(RcppML); library(Matrix); library(dplyr); library(lme4)})
set.seed(1)
outdir <- "results/scrna/cleansed"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
mpdir  <- "results/scrna/metaprograms_cleansed"; dir.create(mpdir, showWarnings = FALSE, recursive = TRUE)

message("Loading GSE206647 ...")
obj <- readRDS("results/scrna/GSE206647_processed.rds")
dat <- GetAssayData(obj, layer = "data")
genes <- rownames(dat)
md <- obj@meta.data
is_tum0 <- md$celltype == "Meningioma"
message(sprintf("Meningioma-annotated cells: %d", sum(is_tum0)))

## ---- 1. PREDEFINED, non-circular contamination flag (same markers as recon) ----
markers <- list(
  Tcell   = c("CD3D","CD3E","CD2","NKG7","GZMK"),
  Myeloid = c("LYZ","C1QA","C1QB","AIF1","TYROBP","LST1","FCER1G"),
  Mast    = c("TPSAB1","TPSB2","CPA3","MS4A2"),
  Bplasma = c("CD79A","MS4A1","MZB1"),
  Endoth  = c("PECAM1","VWF","CLDN5","FLT1"),
  Mural   = c("ACTA2","MYH11","RGS5","PDGFRB"),
  Glial   = c("PLP1","MBP","MAG","OLIG1","OLIG2","PDGFRA")
)
tum_idx <- which(is_tum0)
contam <- rep(FALSE, length(tum_idx))
lineage_hit <- rep("", length(tum_idx))
for (nm in names(markers)) {
  g <- intersect(markers[[nm]], genes)
  if (!length(g)) next
  npos <- Matrix::colSums(dat[g, tum_idx, drop = FALSE] > 1)
  hit <- npos >= 2
  contam <- contam | hit
  lineage_hit[hit & lineage_hit == ""] <- nm
}
# pan-immune backstop: PTPRC(CD45) strongly positive
ptprc <- if ("PTPRC" %in% genes) dat["PTPRC", tum_idx] > 1 else rep(FALSE, length(tum_idx))
contam <- contam | ptprc
keep_cell <- !contam
message(sprintf("Flagged contaminant: %d (%.1f%%); retained clean tumour: %d",
                sum(contam), 100*mean(contam), sum(keep_cell)))

clean_cells <- colnames(obj)[tum_idx[keep_cell]]
gsm_all <- as.character(md$gsm[tum_idx]); grade_all <- as.character(md$grade[tum_idx])

## ---- 2. per-patient cleanse log + dropped patients ----------------------------
pre  <- table(gsm_all)
post <- table(gsm_all[keep_cell])
pat_tab <- data.frame(gsm = names(pre),
                      grade = grade_all[match(names(pre), gsm_all)],
                      n_pre = as.integer(pre),
                      n_post = as.integer(post[names(pre)]))
pat_tab$n_post[is.na(pat_tab$n_post)] <- 0L
pat_tab$removed_pct <- round(100*(1 - pat_tab$n_post/pat_tab$n_pre), 1)
MIN_CELLS <- 100L                      # patients below this drop out of grade tests
pat_tab$dropped <- pat_tab$n_post < MIN_CELLS
write.csv(pat_tab, file.path(outdir, "cleanse_per_patient.csv"), row.names = FALSE)
message("\n=== per-patient cleanse ==="); print(pat_tab, row.names = FALSE)
dropped <- pat_tab$gsm[pat_tab$dropped]
message(sprintf("Patients dropped (< %d clean cells): %s", MIN_CELLS,
                if (length(dropped)) paste(dropped, collapse=", ") else "none"))

## ---- 3. RE-VERIFY per-patient AggrScore vs grade (OLD vs NEW) ------------------
# per-cell AggrScore is fixed; OLD averages over all Meningioma cells, NEW over
# clean cells only -> isolates the contamination effect on the patient-level test.
aggr <- md$AggrScore[tum_idx]
old_pat <- tapply(aggr, gsm_all, mean)
new_pat <- tapply(aggr[keep_cell], gsm_all[keep_cell], mean)
gmap <- setNames(grade_all[match(names(old_pat), gsm_all)], names(old_pat))
rv <- data.frame(gsm = names(old_pat), grade = gmap[names(old_pat)],
                 aggr_old = round(as.numeric(old_pat),3),
                 aggr_new = round(as.numeric(new_pat[names(old_pat)]),3),
                 n_post = pat_tab$n_post[match(names(old_pat), pat_tab$gsm)])
rv$grade_ord <- setNames(c(1,2,3), c("I","II","III"))[rv$grade]
write.csv(rv, file.path(outdir, "reverify_aggrscore_grade.csv"), row.names = FALSE)

tst <- function(score, ord, keeprows) {
  k <- keeprows & !is.na(ord) & is.finite(score)
  ct <- suppressWarnings(cor.test(score[k], ord[k], method="spearman", exact=FALSE))
  av <- summary(aov(score[k] ~ factor(ord[k])))[[1]]$`Pr(>F)`[1]
  c(rho = unname(ct$estimate), p = ct$p.value, anova_p = av, n = sum(k))
}
tumr <- !is.na(rv$grade_ord)                      # tumour grades I/II/III only
old_s <- tst(rv$aggr_old, rv$grade_ord, tumr)
new_all <- tst(rv$aggr_new, rv$grade_ord, tumr)
new_drop <- tst(rv$aggr_new, rv$grade_ord, tumr & !(rv$gsm %in% dropped))
message("\n=== AggrScore ~ grade (patient-level Spearman) ===")
message(sprintf("OLD (all Meningioma, n=%d): rho=%.3f p=%.3g ANOVA=%.3g", old_s["n"], old_s["rho"], old_s["p"], old_s["anova_p"]))
message(sprintf("NEW (cleansed, n=%d):        rho=%.3f p=%.3g ANOVA=%.3g", new_all["n"], new_all["rho"], new_all["p"], new_all["anova_p"]))
message(sprintf("NEW (cleansed, drop<%d, n=%d): rho=%.3f p=%.3g ANOVA=%.3g", MIN_CELLS, new_drop["n"], new_drop["rho"], new_drop["p"], new_drop["anova_p"]))

## ---- 4. RE-VERIFY gene-evidence matrix (script-36 logic on cleansed) -----------
program <- read.csv("docs/Table_S2_aggressiveness_program_genes.csv", check.names=FALSE)
counts_all <- GetAssayData(obj, assay="RNA", layer="counts")
pgenes <- intersect(program$symbol, rownames(counts_all))
scells <- split(clean_cells, as.character(md[clean_cells, "gsm"]))
scells <- scells[sapply(scells, length) >= MIN_CELLS]      # drop tiny patients
lib <- vapply(scells, function(cc) sum(counts_all[, cc, drop=FALSE]), numeric(1))
pb  <- sapply(scells, function(cc) Matrix::rowSums(counts_all[pgenes, cc, drop=FALSE]))
logcpm <- log1p(t(t(pb)/lib) * 1e6)
sgrade <- vapply(scells, function(cc) names(sort(table(as.character(md[cc,"grade"])),decreasing=TRUE))[1], character(1))
gord <- setNames(c(1,2,3), c("I","II","III"))[sgrade]
rho <- pval <- setNames(rep(NA_real_, length(pgenes)), pgenes)
for (g in pgenes) {
  k <- !is.na(gord) & is.finite(logcpm[g,])
  if (sum(k) >= 6 && length(unique(logcpm[g,k])) > 1) {
    tt <- suppressWarnings(cor.test(logcpm[g,k], gord[k], method="spearman", exact=FALSE))
    rho[g] <- unname(tt$estimate); pval[g] <- tt$p.value
  }
}
fdr <- p.adjust(pval, "BH")
ge <- data.frame(symbol=pgenes, rho=round(rho,3), pval=signif(pval,3), fdr=signif(fdr,3))
ge <- ge[order(ge$fdr),]
write.csv(ge, file.path(outdir, "reverify_gene_evidence.csv"), row.names=FALSE)
n_meas <- sum(!is.na(rho)); n_sig <- sum(fdr < 0.10, na.rm=TRUE)
message(sprintf("\n=== gene-evidence (cleansed, %d patients) ===", length(scells)))
message(sprintf("measurable=%d ; FDR<0.10=%d (OLD: 159 measurable, 56 sig)", n_meas, n_sig))
for (g in c("PI3","PITX1","NF2","LTK")) if (g %in% pgenes)
  message(sprintf("  %-6s rho=%.3f FDR=%.3g", g, rho[g], fdr[g]))

## ---- 5. cell-level mixed model (script-31 style, cleansed) ---------------------
tsub <- obj[, clean_cells]
tsub <- CellCycleScoring(tsub, s.features = cc.genes.updated.2019$s.genes,
                         g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
tsub[["pct_ribo"]] <- PercentageFeatureSet(tsub, pattern = "^RP[SL]")
dd <- tsub@meta.data
dd <- dd[dd$grade %in% c("I","II","III") & !(dd$gsm %in% dropped), ]
dd$grade_ord <- setNames(c(1,2,3), c("I","II","III"))[as.character(dd$grade)]
m_un <- lmer(AggrScore ~ grade_ord + (1|gsm), data=dd, REML=FALSE)
m_aj <- lmer(AggrScore ~ grade_ord + S.Score + G2M.Score + nFeature_RNA + pct_ribo + percent.mt + (1|gsm), data=dd, REML=FALSE)
sl_un <- summary(m_un)$coefficients["grade_ord",]
sl_aj <- summary(m_aj)$coefficients["grade_ord",]
writeLines(c(
  sprintf("cleansed cell-level mixed model (n=%d cells, %d patients)", nrow(dd), length(unique(dd$gsm))),
  sprintf("unadjusted grade_ord slope=%.3f  t=%.2f", sl_un["Estimate"], sl_un["t value"]),
  sprintf("adjusted   grade_ord slope=%.3f  t=%.2f", sl_aj["Estimate"], sl_aj["t value"]),
  "(OLD: unadjusted 0.051 t=5.9 ; adjusted 0.040 t=4.7)"
), file.path(outdir, "mixedmodel_slopes.txt"))
message(sprintf("\n=== mixed model (cleansed) unadj slope=%.3f t=%.2f | adj slope=%.3f t=%.2f",
                sl_un["Estimate"], sl_un["t value"], sl_aj["Estimate"], sl_aj["t value"]))

## ---- 6. RE-RUN NMF metaprograms on cleansed cells (script-70 logic) ------------
message("\n=== re-running per-patient NMF on cleansed cells ===")
TOPG <- 50L; KS <- 4:9
prog_list <- list(); prog_meta <- list()
patients <- pat_tab$gsm[!pat_tab$dropped & pat_tab$grade %in% c("I","II","III")]
for (p in patients) {
  cells <- clean_cells[as.character(md[clean_cells,"gsm"]) == p]
  if (length(cells) < 200) next
  sub <- obj[, cells]
  sub <- FindVariableFeatures(sub, selection.method="vst", nfeatures=2000, verbose=FALSE)
  hvg <- VariableFeatures(sub)
  X <- as.matrix(GetAssayData(sub, layer="data")[hvg, ]); X <- X - rowMeans(X); X[X<0] <- 0
  X <- X[rowSums(X) > 0, , drop=FALSE]
  for (k in KS) {
    fit <- tryCatch(RcppML::nmf(X, k=k, verbose=FALSE, seed=1), error=function(e) NULL)
    if (is.null(fit)) next
    W <- fit$w; rownames(W) <- rownames(X)
    for (j in seq_len(ncol(W))) {
      id <- sprintf("%s_k%d_f%d", p, k, j)
      prog_list[[id]] <- names(sort(W[,j], decreasing=TRUE))[1:TOPG]
      prog_meta[[id]] <- data.frame(id=id, patient=p, grade=pat_tab$grade[pat_tab$gsm==p], k=k, factor=j)
    }
  }
}
prog_meta <- bind_rows(prog_meta)
ids <- names(prog_list); n <- length(ids)
J <- matrix(0,n,n,dimnames=list(ids,ids))
for (i in 1:(n-1)) for (j in (i+1):n) {
  jac <- length(intersect(prog_list[[i]],prog_list[[j]]))/length(union(prog_list[[i]],prog_list[[j]]))
  J[i,j] <- jac; J[j,i] <- jac
}
cl <- cutree(hclust(as.dist(1-J), method="average"), h=0.9)
prog_meta$cluster <- cl[prog_meta$id]
cl_tab <- prog_meta %>% group_by(cluster) %>% summarise(n_prog=n(), n_pat=n_distinct(patient), .groups="drop") %>%
  filter(n_pat >= 3) %>% arrange(desc(n_pat))
mp_sig <- list()
for (cc in cl_tab$cluster) {
  members <- prog_meta$id[prog_meta$cluster==cc]
  tabg <- table(unlist(prog_list[members])); thr <- max(2, ceiling(0.25*length(members)))
  sig <- names(sort(tabg[tabg>=thr], decreasing=TRUE)); if (length(sig)<10) sig <- names(sort(tabg,decreasing=TRUE))[1:min(30,length(tabg))]
  mp_sig[[as.character(cc)]] <- head(sig, 50)
}
names(mp_sig) <- paste0("MP", seq_along(mp_sig))
tsc <- AddModuleScore(obj[, clean_cells], features=mp_sig, name="MP", seed=1)
sc_cols <- paste0("MP", seq_along(mp_sig))
mdc <- tsc@meta.data
aggr_cor <- sapply(sc_cols, function(s) cor(mdc[[s]], mdc$AggrScore, method="spearman"))
pat_act <- mdc %>% group_by(gsm) %>% summarise(across(all_of(sc_cols), mean), grade=grade[1], .groups="drop")
pat_act <- pat_act[pat_act$grade %in% c("I","II","III") & !(pat_act$gsm %in% dropped), ]
pat_act$grade_ord <- setNames(c(1,2,3),c("I","II","III"))[pat_act$grade]
grho <- sapply(sc_cols, function(s) suppressWarnings(cor(pat_act[[s]], pat_act$grade_ord, method="spearman")))
gp   <- sapply(sc_cols, function(s) suppressWarnings(cor.test(pat_act[[s]], pat_act$grade_ord, method="spearman", exact=FALSE)$p.value))
summ <- data.frame(metaprogram=sc_cols,
                   n_patients=cl_tab$n_pat[match(setNames(cl_tab$cluster,paste0("MP",seq_len(nrow(cl_tab))))[sc_cols], cl_tab$cluster)],
                   aggr_spearman=round(aggr_cor,3), grade_rho=round(grho,3), grade_p=signif(gp,3)) %>%
  arrange(desc(aggr_spearman))
sig_df <- do.call(rbind, lapply(names(mp_sig), function(m) data.frame(metaprogram=m, rank=seq_along(mp_sig[[m]]), gene=mp_sig[[m]])))
write.csv(summ, file.path(mpdir,"metaprogram_summary.csv"), row.names=FALSE)
write.csv(sig_df, file.path(mpdir,"metaprogram_signatures.csv"), row.names=FALSE)
message("\n=== CLEANSED metaprogram summary (top by AggrScore corr) ===")
print(head(summ, 8), row.names=FALSE)
top_mp <- summ$metaprogram[1]
message(sprintf("Top MP (%s) genes: %s", top_mp, paste(head(mp_sig[[top_mp]],25), collapse=", ")))
saveRDS(list(clean_cells=clean_cells, dropped=dropped), file.path(outdir,"clean_cells.rds"))
message("\nDONE. Outputs in ", outdir, " and ", mpdir)
