# ============================================================
# Figure 2: Highlighted Manhattan Plot
# Trait: Spine BMD (Trajectory GWAS)
# Population: Black (WHI)
# Author: Jongyun Jung
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
})

set.seed(42)

# ------------------------------------------------------------
# 1. Parameters and I/O Paths (generic)
# ------------------------------------------------------------
pop    <- "Black"
trait  <- "Spine_BMD"

input_file  <- paste0("input_data/", pop, "_TrajGWAS_", trait, "_pval.txt")
output_dir  <- "output_figures"
if (!dir.exists("input_data"))  dir.create("input_data",  recursive = TRUE)
if (!dir.exists(output_dir))    dir.create(output_dir,    recursive = TRUE)

# ------------------------------------------------------------
# 2. Load GWAS Results
# ------------------------------------------------------------
gwasResults <- read_delim(input_file, delim = "\t", show_col_types = FALSE) %>%
  mutate(across(c(betapval, taupval, jointpval), as.numeric))

# Compute cumulative chromosome positions
don <- gwasResults %>%
  group_by(chr) %>%
  summarise(chr_len = max(pos, na.rm = TRUE)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  left_join(gwasResults, by = "chr") %>%
  arrange(chr, pos) %>%
  mutate(BPcum = pos + tot)

axisdf <- don %>%
  group_by(chr) %>%
  summarise(center = (max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE)) / 2)

# ------------------------------------------------------------
# 3. Significance Thresholds
# ------------------------------------------------------------
gwas_threshold <- 5e-8
y_threshold    <- -log10(gwas_threshold)
max_y_value    <- max(-log10(don$betapval[don$betapval > 0]), na.rm = TRUE) + 0.5

# ------------------------------------------------------------
# 4. Highlighted SNPs (SPTBN1 & top loci)
# ------------------------------------------------------------
snp_map <- data.frame(
  snpid = c(
    "3:55631684:T:C", "1:222653820:C:T", "4:85530951:C:A", "14:37349831:C:T",
    "3:55222393:G:A", "4:116802297:T:A", "3:55605068:A:G", "2:141798761:C:T",
    "5:10671703:G:T", "3:55585153:G:T"
  ),
  rsID = c(
    "rs4955814", "rs112455822", "rs11939816", "rs1367028",
    "rs112331282", "rs12512743", "rs7621448", "rs4233949", # SPTBN1
    "rs4702719", "rs62251475"
  ),
  stringsAsFactors = FALSE
)

don <- don %>% left_join(snp_map, by = "snpid")
highlight_data <- don %>% filter(rsID %in% snp_map$rsID)

# ------------------------------------------------------------
# 5. Manhattan Plot
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 6. Save Outputs (PNG + PDF)
# ------------------------------------------------------------
ggsave(file.path(output_dir, paste0("Fig2_Manhattan_", pop, "_", trait, ".png")),
       plot = manhattan_plot, width = 14, height = 8, dpi = 600)
ggsave(file.path(output_dir, paste0("Fig2_Manhattan_", pop, "_", trait, ".pdf")),
       plot = manhattan_plot, width = 14, height = 8, dpi = 600, device = cairo_pdf)

message("✅ Figure 2 Manhattan plot successfully saved in 'output_figures/' directory.")
