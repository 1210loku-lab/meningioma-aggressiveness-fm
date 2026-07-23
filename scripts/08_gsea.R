suppressMessages({library(clusterProfiler);library(org.Hs.eg.db);library(AnnotationDbi)}); set.seed(42)
P<-readRDS("results/deg/GSE136661_aggressiveness_program.rds"); res<-P$res
res<-res[!is.na(res$stat),]
ent<-mapIds(org.Hs.eg.db,rownames(res),"ENTREZID","ENSEMBL")
rk<-res$stat; names(rk)<-ent; rk<-rk[!is.na(names(rk))]; rk<-sort(rk[!duplicated(names(rk))],decreasing=TRUE)
cat("ranked genes:",length(rk),"\n")
gse<-gseGO(rk,OrgDb=org.Hs.eg.db,ont="BP",pvalueCutoff=0.1,eps=0,verbose=FALSE)
df<-as.data.frame(gse)
cat("\n=== GSEA GO:BP — top ACTIVATED in high-grade (NES>0) ===\n")
print(head(df[df$NES>0,c("Description","NES","p.adjust")],12))
cat("\n=== top SUPPRESSED in high-grade (NES<0) ===\n")
print(head(df[df$NES<0,c("Description","NES","p.adjust")],12))
write.csv(df,"results/deg/GSEA_GOBP_aggressiveness.csv")
cat("\nsaved results/deg/GSEA_GOBP_aggressiveness.csv (",nrow(df),"terms )\n")
