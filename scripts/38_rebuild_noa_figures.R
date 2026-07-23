# Rebuild the two figures affected by the patient-level evidence audit.
# Uses only saved CSV outputs; does not modify raw data.
suppressMessages({
  library(ggplot2)
  library(patchwork)
})

dir.create("results/figures_pub", showWarnings = FALSE, recursive = TRUE)
theme_set(theme_classic(base_size = 12, base_family = "Arial"))

## Figure 3: make the inferential panel use the primary summed-count pseudobulk score.
cm <- read.csv("results/scrna/GSE206647_cellmeta.csv", check.names = FALSE)
pb <- read.csv("results/deg/R3B_pseudobulk_score_per_patient.csv", check.names = FALSE)

cell_cols <- c(
  Bplasma = "#F8766D", Endothelial = "#C49A00", Meningioma = "#00BA38",
  Mural = "#00BFC4", Myeloid = "#619CFF", Tcell = "#F564E3"
)
cent <- aggregate(cbind(UMAP1, UMAP2) ~ celltype, cm, median)
p3a <- ggplot(cm, aes(UMAP1, UMAP2, colour = celltype)) +
  geom_point(size = 0.08, alpha = 0.55) +
  geom_text(data = cent, aes(label = celltype), colour = "black", size = 3.5,
            fontface = "bold", check_overlap = TRUE) +
  scale_colour_manual(values = cell_cols) +
  coord_equal() +
  labs(title = "GSE206647 cell types (163,897 cells)", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "none")

tum <- cm[cm$celltype == "Meningioma", ]
tum$grade <- factor(tum$grade, levels = c("Normal", "I", "II", "III"))
p3b <- ggplot(tum, aes(grade, AggrScore, fill = grade)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.75) +
  scale_fill_manual(values = c(Normal = "#F8766D", I = "#7CAE00", II = "#00BFC4", III = "#C77CFF")) +
  labs(title = "Cell-level score (descriptive)", x = "WHO grade", y = "Aggressiveness score") +
  theme(legend.position = "none")

pb_t <- pb[pb$grade %in% c("I", "II", "III"), ]
pb_t$grade <- factor(pb_t$grade, levels = c("I", "II", "III"))
ct <- suppressWarnings(cor.test(as.numeric(pb_t$grade), pb_t$score,
                                method = "spearman", exact = FALSE))
ap <- summary(aov(score ~ grade, pb_t))[[1]]["grade", "Pr(>F)"]
p3c <- ggplot(pb_t, aes(grade, score)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey92") +
  geom_jitter(width = 0.10, height = 0, size = 2.5, colour = "#A14B3D") +
  labs(
    title = sprintf("Patient pseudobulk (n=16)\nSpearman rho=%.2f, p=%.1e; ANOVA p=%.1e",
                    unname(ct$estimate), ct$p.value, ap),
    x = "WHO grade", y = "Summed-count pseudobulk program score"
  )

fig3 <- (p3a | p3b | p3c) + plot_layout(widths = c(1.25, 1, 1.05)) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(family="Arial", face="bold", size=14))
ggsave("results/scrna/fig_GSE206647_grade_program.png", fig3,
       width = 17, height = 5.4, dpi = 300)
ggsave("results/scrna/fig_GSE206647_grade_program.pdf", fig3,
       width = 17, height = 5.4)

## Figure 4: replace the discovery-only RF plot with patient-level cross-modal evidence.
ev <- read.csv("results/drug_repurposing/target_evidence_matrix.csv", check.names = FALSE)
ev$direction_short <- ifelse(grepl("^up", ev$direction), "up", "down")
ev$aligned_fdr10 <- !is.na(ev$scrna_grade_fdr) & ev$aligned_scrna_grade & ev$scrna_grade_fdr < 0.10

top <- head(ev[order(-ev$internal_evidence_score), ], 15)
top$symbol <- factor(top$symbol, levels = rev(top$symbol))
p4a <- ggplot(top, aes(internal_evidence_score, symbol, colour = direction_short)) +
  geom_segment(aes(x = 0, xend = internal_evidence_score, yend = symbol),
               colour = "grey75", linewidth = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(up = "#A14B3D", down = "#3B6EA5")) +
  labs(title = "Top cross-modal internal-evidence scores",
       x = "Internal evidence score (ranking only)", y = NULL, colour = "Program direction")

meas <- ev[!is.na(ev$scrna_grade_rho) & !is.na(ev$scrna_grade_fdr), ]
meas$minuslog10_fdr <- -log10(pmax(meas$scrna_grade_fdr, 1e-300))
sel_names <- c("PI3", "PITX1", "NF2", "LTK")
sel <- meas[meas$symbol %in% sel_names, ]
p4b <- ggplot(meas, aes(scrna_grade_rho, minuslog10_fdr)) +
  geom_hline(yintercept = -log10(0.10), linetype = 2, colour = "grey55") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70") +
  geom_point(aes(colour = aligned_fdr10, size = pmax(scrna_tumor_specificity_log2, 0)), alpha = 0.72) +
  geom_text(data = sel, aes(label = symbol), colour = "black", fontface = "bold",
            size = 3.4, vjust = -0.8, check_overlap = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#A14B3D", `FALSE` = "grey70"),
                      labels = c(`TRUE` = "Aligned, FDR<0.10", `FALSE` = "Other")) +
  scale_size_continuous(range = c(1.4, 5), name = "Tumour specificity\n(log2)") +
  labs(title = "Patient-level grade association (159 measurable genes)",
       subtitle = sprintf("%d directionally aligned genes at BH-FDR<0.10", sum(meas$aligned_fdr10)),
       x = "GSE206647 pseudobulk Spearman rho", y = "-log10(BH-FDR)", colour = NULL)

selected <- ev[match(sel_names, ev$symbol), ]
fmt_fdr <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "g", digits = 2))
tiles <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) {
  z <- selected[i, ]
  direction_ok <- isTRUE(z$aligned_scrna_grade)
  data.frame(
    symbol = z$symbol,
    metric = factor(c("Bulk log2FC", "scRNA rho", "scRNA FDR", "Tumour specificity", "External recovery"),
                    levels = c("Bulk log2FC", "scRNA rho", "scRNA FDR", "Tumour specificity", "External recovery")),
    label = c(sprintf("%.2f", z$log2FoldChange), sprintf("%.2f", z$scrna_grade_rho),
              fmt_fdr(z$scrna_grade_fdr), sprintf("%.2f", z$scrna_tumor_specificity_log2),
              sprintf("%d/2", sum(z$recovered_GSE16581, z$recovered_GSE74385))),
    support = c("support", if (direction_ok) "support" else "discordant",
                if (!is.na(z$scrna_grade_fdr) && z$scrna_grade_fdr < 0.10) "support" else "weak",
                if (!is.na(z$scrna_tumor_specificity_log2) && z$scrna_tumor_specificity_log2 > 1) "support" else "weak",
                "support")
  )
}))
tiles$symbol <- factor(tiles$symbol, levels = rev(sel_names))
p4c <- ggplot(tiles, aes(metric, symbol, fill = support)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(support = "#C9E2D0", weak = "#E6E6E6", discordant = "#F3C4C0")) +
  labs(title = "Exact evidence fields for selected candidates", x = NULL, y = NULL, fill = NULL) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")

fig4 <- (p4a | p4b) / p4c + plot_layout(heights = c(1.35, 0.65)) +
  plot_annotation(title = "Cross-modal prioritisation of candidate meningioma aggressiveness markers", tag_levels="a") &
  theme(plot.tag=element_text(family="Arial", face="bold", size=14))
ggsave("results/figures_pub/fig_classical_drivers.png", fig4,
       width = 14, height = 9.5, dpi = 300)
ggsave("results/figures_pub/fig_classical_drivers.pdf", fig4,
       width = 14, height = 9.5)

cat("Rebuilt Figure 3 and Figure 4 from saved result tables.\n")
