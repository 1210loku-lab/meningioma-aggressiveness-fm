# _fig_style.R — shared submission-figure typography + export.
# Sourced by the figure scripts so every panel uses one spec:
#   - panel letters (patchwork tags): UPPERCASE, Arial Bold, 12 pt   (set tag_levels = "A" at call site)
#   - per-panel titles: Arial, 10 pt, horizontally centred over the panel
#   - clean Arial rendering: PDF via cairo_pdf (editable text, real kerning),
#     PNG via ragg::agg_png (systemfonts Arial, high DPI) -- fixes the loose
#     character spacing that the default grDevices pdf()/png() devices produced.
# Letter-spacing (tracking) is intentionally NOT applied here: ggplot2/grid has
# no letter-spacing control. The 1 pt tracking is applied during manual layout
# in the vector editor, where it is a real, uniform control.
suppressMessages({library(ggplot2); library(systemfonts); library(ragg)})

FIG_TAG_PT   <- 12   # panel letter size (pt)
FIG_TITLE_PT <- 10   # per-panel title size (pt)

# Applied to a patchwork with `& fig_style`: centres and resizes every panel
# title and sets the uppercase tag style. Combine with plot_annotation(tag_levels = "A").
fig_style <- theme(
  text          = element_text(family = "Arial"),
  plot.title    = element_text(family = "Arial", size = FIG_TITLE_PT, hjust = 0.5),
  plot.subtitle = element_text(family = "Arial", size = FIG_TITLE_PT - 1, hjust = 0.5),
  plot.tag      = element_text(family = "Arial", face = "bold", size = FIG_TAG_PT)
)

# Export one plot as an editable cairo_pdf plus a crisp ragg PNG preview.
# `stem` is a path without extension. width/height in inches.
save_fig <- function(plot, stem, width, height, dpi = 400) {
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width, height = height, family = "Arial")
  print(plot); grDevices::dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width, height = height, units = "in", res = dpi)
  print(plot); grDevices::dev.off()
  invisible(stem)
}
