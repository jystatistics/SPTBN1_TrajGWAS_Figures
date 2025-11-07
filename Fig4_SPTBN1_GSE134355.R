# ============================================================
# Figure 4: SPTBN1 Expression and Functional Enrichment
# Dataset: GSE134355 (single-cell RNA-seq)
# Author: Jongyun Jung
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggplotify)
  # library(clusterProfiler)
  # library(Seurat)
})

set.seed(42)

# ----------------------------
# 1. Global Parameters & Paths
# ----------------------------
top_n_cat <- 10
outdir <- "../04_Figure"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ----------------------------
# 2. Custom Transformations & Themes
# ----------------------------
neg_log10_trans <- scales::trans_new(
  "neg_log10",
  transform = function(x) -log10(x),
  inverse   = function(x) 10^(-x),
  breaks    = scales::log_breaks(base = 10)
)

color_scale_custom <- scale_color_continuous(
  trans = neg_log10_trans,
  name = expression(-log[10](p.adj)),
  breaks = c(2, 4, 6, 8, 10),
  labels = c(2, 4, 6, 8, 10),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = 0.5,
    barheight = 6,
    label.theme = element_text(size = 9)
  )
)

plot_margin_fix <- theme(plot.margin = unit(c(5.5, 5.5, 5.5, 30), "pt"))

# ============================================================
# 3. GO Enrichment Plots
# ============================================================

p1 <- dotplot(ego_bp, showCategory = top_n_cat, color = "p.adjust") +
  ggtitle("GO BP Enrichment – Correlated Genes") +
  color_scale_custom +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  plot_margin_fix

p2 <- dotplot(ego_de, showCategory = top_n_cat, color = "p.adjust") +
  ggtitle("GO BP Enrichment – High vs. Low SPTBN1") +
  color_scale_custom +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  plot_margin_fix

# ============================================================
# 4. Seurat Visualizations (p3–p6)
# ============================================================

p3 <- VlnPlot(merged_obj, features = "SPTBN1", group.by = "Tissue", pt.size = 0) +
  ggtitle("SPTBN1 Expression Across Tissues") +
  xlab(NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )

p4 <- VlnPlot(merged_obj, features = "SPTBN1", group.by = "Stage", pt.size = 0) +
  ggtitle("SPTBN1 Expression by Developmental Stage") +
  xlab(NULL) +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )

p5 <- DimPlot(merged_obj, group.by = "Stage", reduction = "umap", raster = FALSE) +
  ggtitle("UMAP by Stage") +
  theme(plot.title = element_text(size = 16, face = "bold"))

p6 <- FeaturePlot(merged_obj, features = "SPTBN1", reduction = "umap", raster = FALSE) +
  ggtitle("SPTBN1 Expression on UMAP") +
  theme(plot.title = element_text(size = 16, face = "bold"))

# ============================================================
# 5. Combine Panels (3 × 2 Layout)
# ============================================================

p_combined <- (
  (as.ggplot(p1) + labs(tag = "A")) +
    (as.ggplot(p2) + labs(tag = "B")) +
    (as.ggplot(p3) + labs(tag = "C")) +
    (as.ggplot(p4) + labs(tag = "D")) +
    (as.ggplot(p5) + labs(tag = "E")) +
    (as.ggplot(p6) + labs(tag = "F"))
) +
  plot_layout(ncol = 3, widths = c(0.9, 1, 1)) &
  theme(
    plot.tag = element_text(size = 20, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

# ============================================================
# 6. Save Outputs (PNG + PDF)
# ============================================================
ggsave(file.path(outdir, "Fig4_SPTBN1_GSE134355.png"),
       p_combined, width = 24, height = 12, dpi = 600)
ggsave(file.path(outdir, "Fig4_SPTBN1_GSE134355.pdf"),
       p_combined, width = 24, height = 12, dpi = 600, device = cairo_pdf)
