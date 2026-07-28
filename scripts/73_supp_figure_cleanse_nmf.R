# 73_supp_figure_cleanse_nmf.R
# Minimal supplementary figure for the IJMS resubmission (post SR desk-reject).
# Uses the MODERATE cleanse (script 71; all 16 patients retained, 21% removed).
# Three panels:
#  A) core robustness: per-patient tumour AggrScore vs WHO grade survives stringent
#     lineage-marker purification (OLD all-Meningioma vs cleansed).
#  B) reproducible NMF metaprograms: grade association vs aggressiveness-score
#     correlation; the proliferation program is the stable positive-positive one.
#  C) proliferation metaprogram per-patient activity rises with grade.
# All Arial per project style. Outputs PDF+PNG to results/figures_pub/.
suppressMessages({library(ggplot2); library(dplyr); library(patchwork); library(tidyr)})
th <- theme_classic(base_size=11) +
  theme(text=element_text(family="Arial"), plot.title=element_text(size=10, hjust=0.5),
        axis.title=element_text(size=10))
outdir <- "results/figures_pub"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
cl <- "results/scrna/cleansed"; mp <- "results/scrna/metaprograms_cleansed"

## Panel A: robustness of AggrScore-grade to cleansing
rv <- read.csv(file.path(cl,"reverify_aggrscore_grade.csv"))
rv <- rv[rv$grade %in% c("I","II","III"),]
rv$grade <- factor(rv$grade, levels=c("I","II","III"))
rho_old <- cor(rv$aggr_old, as.integer(rv$grade), method="spearman")
rho_new <- cor(rv$aggr_new, as.integer(rv$grade), method="spearman")
la <- rv %>% select(gsm,grade,aggr_old,aggr_new) %>%
  pivot_longer(c(aggr_old,aggr_new), names_to="set", values_to="score")
la$set <- factor(ifelse(la$set=="aggr_old","All annotated","Lineage-filtered"), levels=c("All annotated","Lineage-filtered"))
pA <- ggplot(la, aes(grade, score, fill=set)) +
  geom_boxplot(width=0.6, outlier.shape=NA, alpha=0.5, position=position_dodge(0.7)) +
  geom_point(aes(color=set), position=position_dodge(0.7), size=1.6) +
  scale_fill_manual(values=c("grey65","#a14b3d")) + scale_color_manual(values=c("grey40","#7d2f22")) +
  labs(x="WHO grade", y="Patient-level tumour\naggressiveness score", fill=NULL, color=NULL,
       title=sprintf("Grade association after lineage-marker exclusion\nPatient-level AddModuleScore: all annotated rho=%.2f; lineage-filtered rho=%.2f", rho_old, rho_new)) +
  th + theme(legend.position="top", legend.text=element_text(size=8))

## Panel B: reproducible metaprograms — grade rho vs PATIENT-LEVEL AggrScore correlation (n=15)
summ <- read.csv(file.path(mp,"metaprogram_summary_patientlevel.csv"))
sig  <- read.csv(file.path(mp,"metaprogram_signatures.csv"))
ccg <- c("MKI67","TOP2A","CDK1","TPX2","BIRC5","CENPF","UBE2C","NUSAP1","PTTG1","TYMS","STMN1")
imm <- c("PTPRC","CD74","HLA-DRA","CCL3","CCL4","TYROBP","IL1B","SRGN","CD3D","NKG7","TPSAB1")
lab <- sapply(summ$metaprogram, function(m){
  g <- sig$gene[sig$metaprogram==m]
  if (sum(g %in% ccg) >= 4) "Proliferation"
  else if (sum(g %in% imm) >= 3) "Immune-like (filter-sensitive)"
  else "Other"
})
summ$class <- lab
ann <- summ[summ$class %in% c("Proliferation","Immune-like (filter-sensitive)"),]
ann$labtxt <- ifelse(ann$class=="Proliferation", sprintf("BH-FDR=%.3f", ann$grade_fdr),
                     sprintf("p=%.3f; BH-FDR=%.3f", ann$grade_p, ann$grade_fdr))
pB <- ggplot(summ, aes(grade_rho, aggr_patient_rho)) +
  geom_hline(yintercept=0, linetype=3, color="grey70") + geom_vline(xintercept=0, linetype=3, color="grey70") +
  geom_point(aes(color=class, size=n_patients)) +
  ggrepel::geom_text_repel(data=ann, aes(label=labtxt), size=2.6, family="Arial", min.segment.length=0, seed=1) +
  scale_color_manual(values=c("Proliferation"="#c0392b","Immune-like (filter-sensitive)"="#2c7fb8","Other"="grey60")) +
  scale_size_continuous(range=c(2,5), name="Patients with\nprogram recovery") +
  labs(x="Grade association (Spearman rho)", y="Patient-level correlation with\naggressiveness score (Spearman rho, n=15)", color=NULL,
       title="Metaprograms in lineage-filtered\ntumour/tumour-like cells") +
  th + theme(legend.position="right", legend.text=element_text(size=8), legend.title=element_text(size=8))

prolif_mp <- summ$metaprogram[sapply(summ$metaprogram, function(m) sum(sig$gene[sig$metaprogram==m] %in% ccg) >= 4)][1]

fig <- (pA | pB) + plot_layout(widths=c(1,1.15)) + plot_annotation(tag_levels="A")
ggsave(file.path(outdir,"FigS5_cleanse_metaprograms.pdf"), fig, width=11, height=4.6, device=cairo_pdf)
ggsave(file.path(outdir,"FigS5_cleanse_metaprograms.png"), fig, width=11, height=4.6, dpi=600)
cat("Proliferation MP:", prolif_mp, " grade_rho=", summ$grade_rho[summ$metaprogram==prolif_mp], "\n")
cat("saved results/figures_pub/FigS5_cleanse_metaprograms.{pdf,png}\n")
