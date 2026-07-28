# 43b_fig2_from_annotated_rds.R — regenerate Figure 2 (GSE183655) by plotting the
# already-annotated Seurat object, bypassing the raw-matrix rebuild in scripts/43
# (whose 10x .mtx inputs are offloaded/unavailable in this environment).
# Reads the annotated RDS via the results/scrna symlink (external volume must be
# mounted). Also caches a compact source-data CSV so future replots need no RDS.
# Typography identical to the other submission figures (scripts/_fig_style.R).
suppressMessages({library(Seurat); library(ggplot2); library(patchwork)})
source("scripts/_fig_style.R")
set.seed(42)
root <- normalizePath(".")

obj <- readRDS("results/scrna/GSE183655_annotated.rds")

# Cache compact source data (UMAP coords + labels) for reproducible replotting.
um <- Embeddings(obj, "umap")
plotdata <- data.frame(cell = colnames(obj), UMAP1 = um[, 1], UMAP2 = um[, 2],
                       celltype = obj$celltype, patient = obj$patient,
                       AggrScore = obj$AggrScore, stringsAsFactors = FALSE)
dir.create(file.path(root, "results", "audit_submission"), recursive = TRUE, showWarnings = FALSE)
write.csv(plotdata, file.path(root, "results", "audit_submission", "GSE183655_Figure2_source_data.csv"),
          row.names = FALSE)

theme_set(theme_classic(base_size = 11, base_family = "Arial"))
p1 <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE, raster = TRUE) +
  labs(title = sprintf("GSE183655 cell types (%s QC cells)", format(ncol(obj), big.mark = ",")),
       x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "none")
p2 <- FeaturePlot(obj, "AggrScore", raster = TRUE) +
  labs(title = "Bulk-derived program score (descriptive)", x = "UMAP 1", y = "UMAP 2")
p3 <- DimPlot(obj, group.by = "patient", raster = TRUE) +
  labs(title = "Patient of origin", x = "UMAP 1", y = "UMAP 2")
fig <- (p1 | p2 | p3) + plot_annotation(tag_levels = "A") & fig_style
save_fig(fig, file.path(root, "results", "scrna", "fig_umap_celltype_v2"), width = 19, height = 5.5)
cat("Figure 2 regenerated from annotated RDS; cells=", ncol(obj), "\n", sep = "")
