#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

root <- normalizePath(getwd(), mustWork = TRUE)
setwd(root)
out <- file.path(root, "results", "figures_pub")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.tag = element_text(family = "Arial", face = "bold", size = 14),
    legend.title = element_text(face = "bold")
  )
grade_cols <- c(I = "#3B6EA5", II = "#D19A2A", III = "#A14B3D")

axis <- read.csv("results/audit_submission/geneformer_patient_level_axis_scores.csv")
axis$grade <- factor(axis$grade, levels = c("I", "II", "III"))
summ <- read.csv("results/audit_submission/geneformer_patient_level_axis_summary.csv")
lopo <- read.csv("results/audit_submission/geneformer_patient_level_lopo_per_repeat.csv")

p5a <- ggplot(axis, aes(grade, primary_axis, colour = grade)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "grey35") +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.2) +
  scale_colour_manual(values = grade_cols) + theme_pub + theme(legend.position = "none") +
  labs(title = "Patient-mean embedding before PCA", x = "WHO grade", y = "Primary PC2 score")

p5b <- ggplot(axis, aes(grade, sensitivity_axis, colour = grade)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "grey35") +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.2) +
  scale_colour_manual(values = grade_cols) + theme_pub + theme(legend.position = "none") +
  labs(title = "Cell PCA before patient averaging", x = "WHO grade", y = "Sensitivity PC4 score")

perm <- rbind(
  data.frame(analysis = "Patient mean\nthen PCA", statistic = c("Observed", "Null 95th"),
             value = c(abs(summ$spearman_grade[1]), summ$null_q95_abs_rho_max[1]), p = summ$empirical_p[1]),
  data.frame(analysis = "Cell PCA\nthen patient mean", statistic = c("Observed", "Null 95th"),
             value = c(abs(summ$spearman_grade[2]), summ$null_q95_abs_rho_max[2]), p = summ$empirical_p[2])
)
perm$analysis <- factor(perm$analysis, levels = unique(perm$analysis))
perm$statistic <- factor(perm$statistic, levels = c("Null 95th", "Observed"))
p5c <- ggplot(perm, aes(analysis, value, fill = statistic)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(data = perm[perm$statistic == "Observed", ],
            aes(label = sprintf("patient-permutation\np = %.4f", p), y = pmin(value + 0.12, 0.98)),
            position = position_dodge(width = 0.72), size = 3.2, family = "Arial") +
  scale_fill_manual(values = c("Null 95th" = "#A7B0B8", "Observed" = "#A14B3D")) +
  coord_cartesian(ylim = c(0, 1.02), clip = "off") + theme_pub +
  labs(title = "Best-of-top-10 patient-level permutation", x = NULL, y = "Maximum |Spearman rho|", fill = NULL)

lopo$model <- factor(lopo$feature_set,
                     levels = c("Geneformer_patient_mean_embedding", "Program_gene_patient_mean_expression"),
                     labels = c("Geneformer", "Program genes"))
p5d <- ggplot(lopo, aes(model, lopo_auc, group = repeat_id)) +
  geom_line(colour = "grey75", linewidth = 0.5) +
  geom_point(aes(colour = model), size = 2.5) +
  stat_summary(aes(group = model), fun = mean, geom = "crossbar", width = 0.48,
               colour = "black", linewidth = 0.55) +
  scale_colour_manual(values = c(Geneformer = "#6F8FB8", `Program genes` = "#B8574D")) +
  coord_cartesian(ylim = c(0.90, 1.005)) + theme_pub + theme(legend.position = "none") +
  labs(title = "Repeated nested leave-one-patient-out", x = NULL, y = "Patient-level AUC")

fig5 <- (p5a | p5b) / (p5c | p5d) +
  plot_annotation(title = "Patient-level Geneformer representation benchmark", tag_levels = "a") &
  theme(text = element_text(family = "Arial"),
        plot.title = element_text(family = "Arial", face = "bold"),
        plot.tag = element_text(family = "Arial", face = "bold", size = 14))
ggsave(file.path(out, "fig_geneformer_patient_level.pdf"), fig5, width = 13, height = 9, device = cairo_pdf)
ggsave(file.path(out, "fig_geneformer_patient_level.png"), fig5, width = 13, height = 9, dpi = 400, bg = "white")

score <- read.csv("results/audit_submission/signature_scoring_sensitivity.csv")
endp <- read.csv("results/audit_submission/GSE74385_endpoint_sensitivity.csv")

fixed <- score[score$discovery_gene_set == "fixed_200_original" &
                 score$outcome %in% c("recurrence", "recurrence_only") &
                 grepl("__(z|rank)__", score$specification), ]
fixed$cohort <- factor(fixed$cohort, levels = c("GSE16581", "GSE74385"))
fixed$score_label <- paste(ifelse(fixed$scoring_method == "rank", "Within-sample rank", "Within-cohort z"),
                           ifelse(fixed$common_platform_genes_only, "common 139", "available genes"), sep = "\n")
pSa <- ggplot(fixed, aes(score_label, recurrence_auc, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey55") +
  scale_fill_manual(values = c(GSE16581 = "#3B6EA5", GSE74385 = "#A14B3D")) +
  coord_cartesian(ylim = c(0.5, 1.0)) + theme_pub +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "Fixed-program scoring sensitivity", x = NULL, y = "Recurrence AUC", fill = NULL)

sizes <- score[score$scoring_method == "rank" & !score$common_platform_genes_only &
                 score$outcome %in% c("recurrence", "recurrence_only") &
                 grepl("^balanced_", score$discovery_gene_set), ]
sizes$n_genes <- as.numeric(sub("balanced_", "", sizes$discovery_gene_set))
pSb <- ggplot(sizes, aes(n_genes, recurrence_auc, colour = cohort)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
  scale_colour_manual(values = c(GSE16581 = "#3B6EA5", GSE74385 = "#A14B3D")) +
  coord_cartesian(ylim = c(0.65, 1.0)) + theme_pub +
  labs(title = "Program-size sensitivity (rank score)", x = "Discovery-program size", y = "Recurrence AUC", colour = NULL)

endp$label <- factor(endp$scenario,
                     levels = rev(endp$scenario),
                     labels = rev(c("Recurrence only", "Malignant progression only", "Composite", "Composite, no grade III")))
pSc <- ggplot(endp, aes(firth_score_OR, label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = firth_score_CI_low, xmax = firth_score_CI_high),
                orientation = "y", width = 0.18, colour = "#3B6EA5") +
  geom_point(size = 2.6, colour = "#A14B3D") +
  scale_x_log10() + theme_pub +
  labs(title = "GSE74385 endpoint sensitivity", x = "Firth OR per 1-SD score (95% CI)", y = NULL)

map <- read.csv("results/audit_submission/geneformer_patient_mapping_audit.csv")
map_text <- readLines("results/audit_submission/geneformer_patient_mapping_audit.txt", warn = FALSE)
direct_line <- grep("Direct row-wise grade agreement", map_text, value = TRUE)
direct_agreement <- as.numeric(sub(".*: ", "", direct_line))
map_stat <- data.frame(
  check = factor(c("Direct row-order\ngrade agreement", "Mapped grade\nagreement", "Unique one-to-one\nmapping"),
                 levels = c("Direct row-order\ngrade agreement", "Mapped grade\nagreement", "Unique one-to-one\nmapping")),
  value = c(direct_agreement,
            as.numeric(length(unique(map$grade)) > 0 && all(map$score_delta <= 1e-12)),
            nrow(unique(map[c("embedding_row", "cells_csv_row")])) / nrow(map))
)
pSd <- ggplot(map_stat, aes(check, value, fill = check)) +
  geom_col(width = 0.62) + geom_text(aes(label = sprintf("%.1f%%", 100 * value)), vjust = -0.4, family = "Arial") +
  coord_cartesian(ylim = c(0, 1.08)) + scale_fill_manual(values = c("#B8574D", "#4E9A67", "#4E9A67")) +
  theme_pub + theme(legend.position = "none") + labs(title = "Historical embedding mapping audit", x = NULL, y = "Proportion")

figS2 <- (pSa | pSb) / (pSc | pSd) +
  plot_annotation(title = "Endpoint, scoring and Geneformer-mapping sensitivity analyses", tag_levels = "a") &
  theme(text = element_text(family = "Arial"),
        plot.title = element_text(family = "Arial", face = "bold"),
        plot.tag = element_text(family = "Arial", face = "bold", size = 14))
ggsave(file.path(out, "fig_review_sensitivity.pdf"), figS2, width = 13, height = 9, device = cairo_pdf)
ggsave(file.path(out, "fig_review_sensitivity.png"), figS2, width = 13, height = 9, dpi = 400, bg = "white")

cat("Wrote revised Fig5 and FigS2 to", out, "\n")
