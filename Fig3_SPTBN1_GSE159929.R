# ============================================================
# Figure 3: SPTBN1 Expression Across Organs and Bone-Related Cell Types
# Dataset: Single-cell RNA-seq (GSE159929)
# Author: Jongyun Jung
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
})

set.seed(42)
outdir <- "../04_Figure"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ============================================================
# 1. Load pre-processed organ-level Seurat objects
# ============================================================
organ_list <- list(
  Bladder = subset_Bladder, Blood = subset_Blood, `Common bile duct` = subset_Common,
  Esophagus = subset_Esophagus, Heart = subset_Heart, `Small intestine` = subset_Intestine,
  Liver = subset_Liver, `Lymph node` = subset_Lymph, Marrow = subset_Marrow,
  Muscle = subset_muscle, Rectum = subset_Rectum, Skin = subset_Skin,
  Spleen = subset_Spleen, Stomach = subset_Stomach, Trachea = subset_Trachea
)

# ============================================================
# 2. Global theme
# ============================================================
common_theme <- theme(
  plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
  axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
  axis.text.y = element_text(size = 12),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  legend.position = "none",
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
  panel.grid = element_blank(),
  plot.margin = margin(t = 5, r = 5, b = 5, l = 35, unit = "pt")
)

# ============================================================
# 3. Identify bone-related cell types
# ============================================================
bone_keywords <- c("osteoblast", "osteoclast", "chondro", "stromal", "macrophage", "MSC", "BMSC")

# ============================================================
# 4. Generate violin plots and expression summaries
# ============================================================
plot_list <- list()
expr_summary <- data.frame()

for (org in names(organ_list)) {
  obj <- organ_list[[org]]
  df <- FetchData(obj, vars = c("SPTBN1", "HCA"))
  colnames(df) <- c("Expression", "CellType")
  df$Organ <- org
  
  top_ct <- df %>%
    group_by(CellType) %>%
    summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
    filter(mean_expr > 0.5 | grepl(paste(bone_keywords, collapse = "|"), CellType, ignore.case = TRUE)) %>%
    slice_max(mean_expr, n = 10) %>%
    pull(CellType)
  
  # Always retain osteoblast/osteoclast if present
  top_ct <- unique(c(top_ct, intersect(c("osteoblast", "osteoclast"), unique(df$CellType))))
  df_filt <- df %>% filter(CellType %in% top_ct)
  
  p <- ggplot(df_filt, aes(x = CellType, y = Expression, fill = CellType)) +
    geom_violin(trim = TRUE) +
    geom_jitter(size = 0.4, alpha = 0.5, width = 0.1) +
    ggtitle(org) +
    ylab("Expression Level") +
    common_theme +
    theme(axis.title.y = element_blank())
  plot_list[[org]] <- p
  
  # Expression summary
  expr_summary <- bind_rows(
    expr_summary,
    df_filt %>%
      group_by(Organ, CellType) %>%
      summarise(
        MeanExpr = mean(Expression, na.rm = TRUE),
        MedianExpr = median(Expression, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
  )
}

# ============================================================
# 5. Combine representative organs (top 6)
# ============================================================
selected_organs <- c("Marrow", "Muscle", "Liver", "Heart", "Skin", "Trachea")
combined_plot <- plot_grid(plotlist = plot_list[selected_organs], ncol = 3, align = "hv")
combined_plot <- ggdraw(combined_plot) +
  draw_label("Expression Level", x = 0.015, y = 0.5, vjust = 0.5, angle = 90,
             size = 18, fontface = "bold")

# Save violin panel
ggsave(file.path(outdir, "Fig3_SPTBN1_scRNA_violin.png"),
       combined_plot, width = 16, height = 10, dpi = 600)
ggsave(file.path(outdir, "Fig3_SPTBN1_scRNA_violin.pdf"),
       combined_plot, width = 16, height = 10, dpi = 600, device = cairo_pdf)
cat("✅ Violin plots saved.\n")

# ============================================================
# 6. Save summary table
# ============================================================
write.csv(expr_summary, file.path(outdir, "Fig3_SPTBN1_expression_summary.csv"), row.names = FALSE)
