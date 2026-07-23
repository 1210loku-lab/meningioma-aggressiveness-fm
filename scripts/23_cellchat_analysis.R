suppressMessages({library(Seurat); library(CellChat); library(patchwork)})
outdir <- "results/descriptive"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
cat("--- Descriptive CellChat analysis (GSE183655) ---\n")

o <- readRDS("results/scrna/GSE183655_annotated.rds")
Idents(o) <- "celltype"
o <- subset(o, downsample=500)
data.input <- GetAssayData(o, assay="RNA", layer="data")
labels <- Idents(o)
meta <- data.frame(labels=labels, row.names=names(labels))
cellchat <- createCellChat(object=data.input, meta=meta, group.by="labels")
cellchat@DB <- subsetDB(CellChatDB.human, search="Secreted Signaling")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells=10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file.path(outdir, "23_cellchat_results.rds"))

pdf(file.path(outdir, "23_fig_cellchat_network.pdf"), width=8, height=8, family="Arial")
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow=c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight=groupSize, weight.scale=TRUE,
                 label.edge=FALSE, title.name="Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight=groupSize, weight.scale=TRUE,
                 label.edge=FALSE, title.name="Interaction weights/strength")
dev.off()

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name="netP")
pdf(file.path(outdir, "23_fig_cellchat_pathways.pdf"), width=12, height=8, family="Arial")
for (pathway in head(cellchat@netP$pathways, 3)) {
  netVisual_aggregate(cellchat, signaling=pathway, layout="circle")
}
dev.off()
