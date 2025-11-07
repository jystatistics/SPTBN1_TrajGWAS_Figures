# ============================================================
# Figure 5: Mendelian Randomization Forest Plots (MOF & Hip Fx)
# Author: Jongyun Jung
# ============================================================

# --- 1. Packages ---
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# --- 1A. Compute Pooled (IVW) Causal Effects -------------------------------
# Function to compute pooled HR using inverse-variance weighted method
ivw_pool <- function(hr, lower, upper) {
  # convert HR to log scale
  beta <- log(hr)
  se <- (log(upper) - log(lower)) / (2 * 1.96)
  w <- 1 / se^2
  
  beta_ivw <- sum(beta * w) / sum(w)
  se_ivw <- sqrt(1 / sum(w))
  
  ci_low <- exp(beta_ivw - 1.96 * se_ivw)
  ci_high <- exp(beta_ivw + 1.96 * se_ivw)
  hr_ivw <- exp(beta_ivw)
  
  return(data.frame(HR = hr_ivw, Lower = ci_low, Upper = ci_high))
}

# --- Compute pooled HRs for each sex and overall ---------------------------
female_ivw_mof   <- ivw_pool(
  hr    = c(1.065, 1.045, 1.038),
  lower = c(1.030, 1.010, 1.005),
  upper = c(1.100, 1.080, 1.070)
)
male_ivw_mof     <- ivw_pool(
  hr    = c(1.015, 1.012, 1.010),
  lower = c(0.985, 0.980, 0.970),
  upper = c(1.045, 1.045, 1.050)
)
overall_ivw_mof  <- ivw_pool(
  hr    = c(female_ivw_mof$HR, male_ivw_mof$HR),
  lower = c(female_ivw_mof$Lower, male_ivw_mof$Lower),
  upper = c(female_ivw_mof$Upper, male_ivw_mof$Upper)
)

female_ivw_hipfx <- ivw_pool(
  hr    = c(1.055, 1.038, 1.032),
  lower = c(1.020, 1.000, 0.998),
  upper = c(1.090, 1.075, 1.065)
)
male_ivw_hipfx   <- ivw_pool(
  hr    = c(1.010, 1.008, 1.006),
  lower = c(0.980, 0.950, 0.960),
  upper = c(1.040, 1.050, 1.050)
)
overall_ivw_hipfx <- ivw_pool(
  hr    = c(female_ivw_hipfx$HR, male_ivw_hipfx$HR),
  lower = c(female_ivw_hipfx$Lower, male_ivw_hipfx$Lower),
  upper = c(female_ivw_hipfx$Upper, male_ivw_hipfx$Upper)
)

# --- 1B. Replace pooled_data table automatically with computed IVW results -----
pooled_data <- data.frame(
  Subgroup = c("Overall Pooled Effect", "Pooled Female Effect", "Pooled Male Effect"),
  Sex = c("Overall", "Female", "Male"),
  
  HR_MOF = c(overall_ivw_mof$HR, female_ivw_mof$HR, male_ivw_mof$HR),
  Lower_MOF = c(overall_ivw_mof$Lower, female_ivw_mof$Lower, male_ivw_mof$Lower),
  Upper_MOF = c(overall_ivw_mof$Upper, female_ivw_mof$Upper, male_ivw_mof$Upper),
  
  HR_HipFx = c(overall_ivw_hipfx$HR, female_ivw_hipfx$HR, male_ivw_hipfx$HR),
  Lower_HipFx = c(overall_ivw_hipfx$Lower, female_ivw_hipfx$Lower, male_ivw_hipfx$Lower),
  Upper_HipFx = c(overall_ivw_hipfx$Upper, female_ivw_hipfx$Upper, male_ivw_hipfx$Upper),
  
  pval = c(0.045, 0.018, 0.27)
)

# --- 2. Helper Functions ------------------------------------------------------

prepare_plot_data <- function(outcome) {
  hr_cols <- paste0(c("HR_", "Lower_", "Upper_"), outcome)
  base <- bind_rows(validation_data, pooled_data)
  base %>%
    mutate(
      HR_Causal = !!sym(hr_cols[1]),
      HR_Lower  = !!sym(hr_cols[2]),
      HR_Upper  = !!sym(hr_cols[3]),
      Sig = ifelse(pval < 0.05, "Significant", "Non-Significant"),
      Label = sprintf("%.3f (%.3f, %.3f)", HR_Causal, HR_Lower, HR_Upper),
      Subgroup = factor(
        Subgroup,
        levels = c(
          "Overall Pooled Effect", "Pooled Female Effect",
          rev(filter(validation_data, Sex == "Female")$Subgroup),
          "Pooled Male Effect", rev(filter(validation_data, Sex == "Male")$Subgroup)
        )
      )
    )
}

add_plot_shapes <- function(df) {
  df %>%
    mutate(
      Plot_Shape = case_when(
        Subgroup == "Overall Pooled Effect" ~ 23,
        Subgroup == "Pooled Female Effect" ~ 21,
        Subgroup == "Pooled Male Effect"   ~ 22,
        TRUE ~ 15
      ),
      Plot_Size = ifelse(Subgroup %in%
                           c("Overall Pooled Effect", "Pooled Female Effect", "Pooled Male Effect"), 4, 2)
    )
}

generate_forest_plot <- function(data, outcome_label, x_limit_max = 1.15) {
  label_x <- x_limit_max * 1.03
  bold_labels <- c("Overall Pooled Effect", "Pooled Female Effect", "Pooled Male Effect")
  
  ggplot(data, aes(y = Subgroup, x = HR_Causal, xmin = HR_Lower, xmax = HR_Upper, color = Sig)) +
    geom_pointrange(aes(shape = factor(Plot_Shape), size = Plot_Size), fatten = 2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_text(aes(label = Label, x = label_x), hjust = 0, size = 5.5, color = "black") +
    scale_color_manual(values = c("Significant" = "#D62728", "Non-Significant" = "#1F77B4")) +
    scale_shape_manual(values = c("23" = 23, "21" = 21, "22" = 22, "15" = 15)) +
    scale_size_identity() +
    scale_x_continuous(
      trans = 'log10',
      limits = c(0.95, x_limit_max),
      breaks = seq(1.0, x_limit_max, by = 0.05),
      minor_breaks = NULL
    ) +
    labs(x = paste("Hazard Ratio for", outcome_label), y = NULL) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(
        face = ifelse(levels(data$Subgroup) %in% bold_labels, "bold", "plain"),
        color = "black", size = 20
      ),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linetype = "dotted", color = "gray80"),
      plot.margin = unit(c(0.5, 5, 0.5, 0.5), "cm")
    )
}

# --- 3. Generate Clean Data for Plots ----------------------------------------
plot_data_mof   <- add_plot_shapes(prepare_plot_data("MOF"))
plot_data_hipfx <- add_plot_shapes(prepare_plot_data("HipFx"))

# --- 4. Build and Combine Figures --------------------------------------------
plot_mof   <- generate_forest_plot(plot_data_mof, "Major Osteoporotic Fracture")
plot_hipfx <- generate_forest_plot(plot_data_hipfx, "Hip Fracture")

combined_plot <- (plot_mof | plot_hipfx) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 22, colour = "black"),
      plot.tag.position = c(0.02, 1.05)
    )
  )

# --- 5. Save Outputs (PNG + PDF) ---------------------------------------------
output_dir <- "../04_Figure"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "Fig5_MR_ForestPlot.png"),
       plot = combined_plot, width = 22, height = 8, dpi = 600, units = "in")
ggsave(file.path(output_dir, "Fig5_MR_ForestPlot.pdf"),
       plot = combined_plot, width = 22, height = 8, dpi = 600, units = "in", device = cairo_pdf)

cat("✅ Figure 5 saved successfully to:", output_dir, "\n")
