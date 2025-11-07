# ============================================================
# Figure 1: Highlighted Manhattan Plot
# Trait: Hip BMD (Trajectory GWAS)
# Population: Black (WHI)
# Author: Jongyun Jung
# ============================================================

# -----------------------------
# 1. Load Required Packages
# -----------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
})

set.seed(42)

# -----------------------------
# 2. User-Defined Parameters
# -----------------------------
pop <- "Black"
trait <- "Hip_BMD"

# Define I/O (generic, portable)
input_file  <- paste0("input_data/", pop, "_TrajGWAS_", trait, "_pval.txt")
output_file <- paste0("Fig1_ManhattanPlot_", pop, "_", trait, ".png")

if (!dir.exists("input_data")) dir.create("input_data")
if (!dir.exists("output_figures")) dir.create("output_figures")

# -----------------------------
# 3. Read and Prepare GWAS Data
# -----------------------------
gwasResults <- read_delim(input_file, delim = "\t", show_col_types = FALSE)

# Ensure numeric
gwasResults <- gwasResults %>%
  mutate(across(c(betapval, taupval, jointpval), as.numeric))

# Compute cumulative base-pair positions for chromosomes
don <- gwasResults %>%
  group_by(chr) %>%
  summarise(chr_len = max(pos, na.rm = TRUE)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  left_join(gwasResults, ., by = "chr") %>%
  arrange(chr, pos) %>%
  mutate(BPcum = pos + tot)

axisdf <- don %>%
  group_by(chr) %>%
  summarize(center = (max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE)) / 2)

# -----------------------------
# 4. Define Significance Line
# -----------------------------
gwas_threshold <- 5e-8
y_threshold <- -log10(gwas_threshold)
max_y_value <- max(-log10(don$betapval[don$betapval > 0]), na.rm = TRUE) + 0.5

# -----------------------------
# 5. Annotate Top SNPs (Example)
# -----------------------------
snp_map <- data.frame(
  snpid = c(
    "19:57601562:T:C", "7:90007258:T:A", "19:53333615:C:T",
    "21:46006705:C:T", "4:167305701:A:C", "2:10017268:G:T",
    "9:78263860:C:T", "13:95396127:T:G", "1:237271291:T:C", "10:76803567:T:G"
  ),
  rsID = c(
    "rs7251124", "rs7779628", "rs1974831", "rs113083240",
    "rs1565649", "rs4233949", "rs744339", "rs7334743",
    "rs74701136", "rs16931938"
  ),
  stringsAsFactors = FALSE
)

don <- don %>% left_join(snp_map, by = "snpid")
highlight_data <- don %>% filter(rsID %in% snp_map$rsID)

# -----------------------------
# 6. Manhattan Plot
# -----------------------------
manhattan_plot <- ggplot(don, aes(x = BPcum, y = -log10(betapval))) +
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 1.2) +
  geom_point(data = highlight_data %>% filter(betadir == 1),
             aes(x = BPcum, y = -log10(betapval)),
             shape = 24, size = 4, color = "black", fill = "red") +
  geom_point(data = highlight_data %>% filter(betadir == -1),
             aes(x = BPcum, y = -log10(betapval)),
             shape = 25, size = 4, color = "black", fill = "blue") +
  geom_text_repel(data = highlight_data, aes(label = rsID),
                  size = 5, color = "black", max.overlaps = Inf, nudge_y = 0.5) +
  geom_hline(yintercept = y_threshold, linetype = "dotted", color = "red", linewidth = 0.8) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max_y_value)) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P))),
       title = paste("Manhattan Plot –", trait, "(", pop, ")")) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

# -----------------------------
# 7. Save Output
# -----------------------------
ggsave(filename = file.path("output_figures", paste0("Fig1_Manhattan_", pop, "_", trait, ".png")),
       plot = manhattan_plot, width = 14, height = 8, dpi = 600)
ggsave(filename = file.path("output_figures", paste0("Fig1_Manhattan_", pop, "_", trait, ".pdf")),
       plot = manhattan_plot, width = 14, height = 8, dpi = 600, device = cairo_pdf)

message("✅ Figure 1 Manhattan plot successfully saved in 'output_figures/' directory.")
