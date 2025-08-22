################################################################################
# SIMPLE MR VISUALIZATION SCRIPT (No External Dependencies)
# EndoMR-PREG: Association of Genetic Liability to Endometriosis with Pregnancy Outcomes
# 
# Purpose: Generate publication-ready figures with basic R + ggplot2
# Author: Jonas Vibert, MD
# Date: August 2025
# 
# Requirements: Only TwoSampleMR, dplyr, ggplot2, stringr (standard packages)
################################################################################

# Load only essential libraries
library(TwoSampleMR)
library(dplyr)
library(ggplot2)
library(stringr)

# Simple helper functions
clean_name <- function(x) {
  str_replace_all(x, "[^A-Za-z0-9]", "_")
}

################################################################################
# CONFIGURATION
################################################################################

message("📊 Starting simple MR visualization...")

# Check if harmonized data exists
if (!exists("dat")) {
  if (file.exists("results/harmonised_rahmioglu_bpo.csv")) {
    dat <- read.csv("results/harmonised_rahmioglu_bpo.csv")
    message("✅ Loaded harmonized data from file")
  } else {
    stop("❌ Harmonized data 'dat' not found. Please run harmonization script first.")
  }
}

# Create organized results directory structure
results_base_dir <- "visual_results"
plots_dir <- file.path(results_base_dir, "plots")
tables_dir <- file.path(results_base_dir, "tables")

# Create directories
if (!dir.exists(results_base_dir)) dir.create(results_base_dir, recursive = TRUE)
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)

message("📁 Created organized directory structure:")
message("   📊 ", results_base_dir, "/")
message("     📈 plots/")
message("     📋 tables/")

# Bonferroni-corrected significant outcomes (P < 0.001)
bonferroni_significant <- c(
  "finngen_R12_N14_FEMALEINFERT_filtered",     # Female infertility
  "finngen_R12_O15_PLAC_PRAEVIA_filtered"      # Placenta praevia
)

# Nominally significant outcomes (P < 0.05 but > 0.001)
nominally_significant <- c(
  "5sBPKR",                                    # Premature rupture of membranes
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered", # Placental abruption
  "mxO2D8",                                    # Elective caesarean delivery
  "z2noFK",                                    # Low Apgar score at 1 minute
  "O2sj1s"                                     # Preterm birth (all)
)

# All outcomes of interest
all_outcomes <- c(bonferroni_significant, nominally_significant)

# Clean outcome labels
outcome_labels <- list(
  "finngen_R12_N14_FEMALEINFERT_filtered" = "Female infertility",
  "finngen_R12_O15_PLAC_PRAEVIA_filtered" = "Placenta praevia",
  "5sBPKR" = "Premature rupture of membranes",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered" = "Placental abruption",
  "mxO2D8" = "Elective caesarean delivery",
  "z2noFK" = "Low Apgar score at 1 minute", 
  "O2sj1s" = "Preterm birth (all)"
)

################################################################################
# RUN MR ANALYSES
################################################################################

message("🔬 Running MR analyses...")

# Main MR analysis
mr_results <- mr(dat)

# Single SNP analysis
single_snp_results <- mr_singlesnp(dat)

# Leave-one-out analysis  
loo_results <- mr_leaveoneout(dat)

################################################################################
# PUBLICATION THEME
################################################################################

theme_clean <- theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey95", size = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 15)),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

# Color scheme
colors_significance <- c(
  "Bonferroni significant" = "#D73027",
  "Nominally significant" = "#FDB863"
)

################################################################################
# SCATTER PLOTS
################################################################################

message("📈 Creating scatter plots...")

for (outcome_id in all_outcomes) {
  
  # Check if outcome exists in data
  outcome_data <- dat[dat$id.outcome == outcome_id, ]
  outcome_mr <- mr_results[mr_results$id.outcome == outcome_id, ]
  
  if (nrow(outcome_data) == 0 || nrow(outcome_mr) == 0) {
    message("⚠️ Skipping ", outcome_id, " - no data")
    next
  }
  
  # Clean method names
  outcome_mr$method_clean <- case_when(
    outcome_mr$method == "Inverse variance weighted" ~ "IVW",
    outcome_mr$method == "MR Egger" ~ "MR Egger",
    outcome_mr$method == "Weighted median" ~ "Weighted median",
    TRUE ~ outcome_mr$method
  )
  
  # Get significance level
  sig_level <- ifelse(outcome_id %in% bonferroni_significant, 
                     "Bonferroni significant", "Nominally significant")
  
  # Get clean outcome name
  outcome_name <- outcome_labels[[outcome_id]]
  if (is.null(outcome_name)) outcome_name <- outcome_id
  
  # Create scatter plot
  p <- ggplot(outcome_data, aes(x = beta.exposure, y = beta.outcome)) +
    
    # Add SNP points
    geom_point(aes(size = 1/se.outcome^2), 
               color = "grey50", alpha = 0.7, shape = 21, fill = "white") +
    
    # Add regression lines
    geom_abline(data = outcome_mr, 
                aes(intercept = 0, slope = b, color = method_clean, linetype = method_clean),
                size = 1, alpha = 0.9) +
    
    # Scales
    scale_color_manual(values = c("IVW" = "#2166AC", "MR Egger" = "#D73027", 
                                 "Weighted median" = "#762A83"),
                      name = "MR Method") +
    scale_linetype_manual(values = c("IVW" = "solid", "MR Egger" = "dashed", 
                                   "Weighted median" = "dotted"),
                         name = "MR Method") +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    
    # Labels
    labs(
      title = paste("MR Analysis:", outcome_name),
      subtitle = paste("SNP effects on endometriosis vs. outcome (", sig_level, ")", sep = ""),
      x = "SNP effect on endometriosis (β)",
      y = "SNP effect on outcome (β)",
      caption = "Point size ∝ inverse variance"
    ) +
    
    theme_clean
  
  # Save plot
  filename_base <- clean_name(outcome_name)
  
  ggsave(
    filename = file.path(plots_dir, paste0("scatter_", filename_base, ".png")),
    plot = p, width = 10, height = 8, dpi = 300, bg = "white"
  )
  
  ggsave(
    filename = file.path(plots_dir, paste0("scatter_", filename_base, ".pdf")),
    plot = p, width = 10, height = 8, bg = "white"
  )
  
  message("💾 Saved scatter plot: ", outcome_name)
}

################################################################################
# SUMMARY FOREST PLOT
################################################################################

message("🌳 Creating summary forest plot...")

# Prepare data for forest plot
forest_data <- mr_results %>%
  filter(id.outcome %in% all_outcomes, method == "Inverse variance weighted") %>%
  mutate(
    outcome_clean = sapply(id.outcome, function(x) {
      label <- outcome_labels[[x]]
      if (is.null(label)) x else label
    }),
    OR = exp(b),
    OR_lower = exp(b - 1.96 * se),
    OR_upper = exp(b + 1.96 * se),
    significance_level = ifelse(id.outcome %in% bonferroni_significant,
                               "Bonferroni significant", "Nominally significant"),
    CI_text = sprintf("%.2f (%.2f-%.2f)", OR, OR_lower, OR_upper),
    p_text = case_when(
      pval < 0.001 ~ sprintf("P = %.1e", pval),
      pval < 0.01 ~ sprintf("P = %.3f", pval),
      TRUE ~ sprintf("P = %.2f", pval)
    )
  ) %>%
  arrange(desc(significance_level == "Bonferroni significant"), pval)

# Reorder factor levels
forest_data$outcome_clean <- factor(forest_data$outcome_clean, 
                                   levels = forest_data$outcome_clean)

################################################################################
# SAVE SUMMARY TABLES
################################################################################

message("📋 Creating and saving summary tables...")

# 1. Main results summary table
main_results_table <- forest_data %>%
  select(
    Outcome = outcome_clean,
    "Significance Level" = significance_level,
    "Number of SNPs" = nsnp,
    "Effect Size (β)" = b,
    "Standard Error" = se,
    "Odds Ratio" = OR,
    "95% CI Lower" = OR_lower,
    "95% CI Upper" = OR_upper,
    "P-value" = pval,
    "OR (95% CI)" = CI_text
  ) %>%
  mutate(
    "Effect Size (β)" = round(`Effect Size (β)`, 4),
    "Standard Error" = round(`Standard Error`, 4),
    "Odds Ratio" = round(`Odds Ratio`, 3),
    "95% CI Lower" = round(`95% CI Lower`, 3),
    "95% CI Upper" = round(`95% CI Upper`, 3),
    "P-value" = case_when(
      `P-value` < 0.001 ~ "< 0.001",
      `P-value` < 0.01 ~ sprintf("%.3f", `P-value`),
      TRUE ~ sprintf("%.2f", `P-value`)
    )
  )

write.csv(main_results_table, file.path(tables_dir, "main_results_summary.csv"), 
          row.names = FALSE)
message("💾 Saved main results table: main_results_summary.csv")

# 2. Full MR methods comparison table
full_mr_table <- mr_results %>%
  filter(id.outcome %in% all_outcomes) %>%
  mutate(
    outcome_clean = sapply(id.outcome, function(x) {
      label <- outcome_labels[[x]]
      if (is.null(label)) x else label
    }),
    OR = exp(b),
    OR_lower = exp(b - 1.96 * se),
    OR_upper = exp(b + 1.96 * se),
    OR_CI = sprintf("%.2f (%.2f-%.2f)", OR, OR_lower, OR_upper)
  ) %>%
  select(
    Outcome = outcome_clean,
    Method = method,
    "Number of SNPs" = nsnp,
    "Effect Size (β)" = b,
    "Standard Error" = se,
    "OR (95% CI)" = OR_CI,
    "P-value" = pval
  ) %>%
  mutate(
    "Effect Size (β)" = round(`Effect Size (β)`, 4),
    "Standard Error" = round(`Standard Error`, 4),
    "P-value" = case_when(
      `P-value` < 0.001 ~ "< 0.001",
      `P-value` < 0.01 ~ sprintf("%.3f", `P-value`),
      TRUE ~ sprintf("%.2f", `P-value`)
    )
  ) %>%
  arrange(Outcome, Method)

write.csv(full_mr_table, file.path(tables_dir, "all_methods_comparison.csv"), 
          row.names = FALSE)
message("💾 Saved full methods table: all_methods_comparison.csv")

# 3. Significance summary table  
significance_summary <- data.frame(
  "Significance Level" = c("Bonferroni significant (P < 0.001)", 
                          "Nominally significant (P < 0.05)"),
  "Number of Outcomes" = c(length(bonferroni_significant), 
                          length(nominally_significant)),
  "Outcomes" = c(paste(sapply(bonferroni_significant, function(x) outcome_labels[[x]]), collapse = "; "),
                paste(sapply(nominally_significant, function(x) outcome_labels[[x]]), collapse = "; ")),
  stringsAsFactors = FALSE
)

write.csv(significance_summary, file.path(tables_dir, "significance_summary.csv"), 
          row.names = FALSE)
message("💾 Saved significance summary: significance_summary.csv")

# Create forest plot
forest_plot <- ggplot(forest_data, aes(y = outcome_clean, x = OR)) +
  
  # Reference line
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  
  # Confidence intervals
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper, color = significance_level),
                height = 0.2, size = 1, alpha = 0.8) +
  
  # Point estimates
  geom_point(aes(color = significance_level, size = significance_level, 
                shape = significance_level), 
            fill = "white", stroke = 1.5) +
  
  # Effect size text
  geom_text(aes(x = OR_upper + 0.15, label = CI_text),
            hjust = 0, size = 3.5, fontface = "bold") +
  
  # P-value text
  geom_text(aes(x = OR_upper + 0.15, label = p_text, 
               y = as.numeric(outcome_clean) - 0.2),
            hjust = 0, size = 3, color = "grey40") +
  
  # Scales
  scale_color_manual(values = colors_significance, name = "Significance Level") +
  scale_size_manual(values = c("Bonferroni significant" = 4, 
                              "Nominally significant" = 3), guide = "none") +
  scale_shape_manual(values = c("Bonferroni significant" = 22, 
                               "Nominally significant" = 21), guide = "none") +
  scale_x_log10(limits = c(0.7, 4.5),
                breaks = c(0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0),
                labels = c("0.8", "1.0", "1.25", "1.5", "2.0", "2.5", "3.0")) +
  
  # Labels
  labs(
    title = "Mendelian Randomization: Endometriosis and Pregnancy Outcomes",
    subtitle = "Causal effect estimates (Inverse Variance Weighted method)",
    x = "Odds Ratio (95% CI)",
    y = "",
    caption = "Red squares = Bonferroni significant (P < 0.001); Orange circles = Nominally significant (P < 0.05)"
  ) +
  
  theme_clean +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "bottom"
  )

# Save forest plot
ggsave(
  filename = file.path(plots_dir, "summary_forest_plot.png"),
  plot = forest_plot, width = 12, height = 8, dpi = 300, bg = "white"
)

ggsave(
  filename = file.path(plots_dir, "summary_forest_plot.pdf"),
  plot = forest_plot, width = 12, height = 8, bg = "white"
)

message("💾 Saved summary forest plot")

################################################################################
# INDIVIDUAL FOREST PLOTS
################################################################################

message("🌲 Creating individual forest plots...")

for (outcome_id in all_outcomes) {
  
  # Filter single SNP data
  snp_data <- single_snp_results[single_snp_results$id.outcome == outcome_id, ]
  
  if (nrow(snp_data) == 0) {
    message("⚠️ Skipping forest plot for ", outcome_id, " - no data")
    next
  }
  
  # Prepare data
  snp_data$SNP_clean <- ifelse(snp_data$SNP == "All", "Summary (IVW)", snp_data$SNP)
  snp_data$OR <- exp(snp_data$b)
  snp_data$OR_lower <- exp(snp_data$b - 1.96 * snp_data$se)
  snp_data$OR_upper <- exp(snp_data$b + 1.96 * snp_data$se)
  snp_data$is_summary <- snp_data$SNP == "All"
  
  # Order data
  snp_data <- snp_data[order(snp_data$is_summary, snp_data$OR), ]
  snp_data$SNP_clean <- factor(snp_data$SNP_clean, levels = snp_data$SNP_clean)
  
  # Get names and significance
  outcome_name <- outcome_labels[[outcome_id]]
  if (is.null(outcome_name)) outcome_name <- outcome_id
  sig_level <- ifelse(outcome_id %in% bonferroni_significant, 
                     "Bonferroni significant", "Nominally significant")
  
  # Create forest plot
  p <- ggplot(snp_data, aes(y = SNP_clean, x = OR)) +
    
    # Reference line
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    
    # Confidence intervals
    geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper, 
                      color = is_summary, size = is_summary),
                  height = 0.3, alpha = 0.8) +
    
    # Points
    geom_point(aes(color = is_summary, size = is_summary, shape = is_summary)) +
    
    # Scales
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#D73027"), guide = "none") +
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2), guide = "none") +
    scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23), guide = "none") +
    scale_x_log10() +
    
    # Labels
    labs(
      title = paste("Forest Plot:", outcome_name),
      subtitle = paste("Individual SNP effects (", sig_level, ")", sep = ""),
      x = "Odds Ratio (95% CI)",
      y = "",
      caption = "Diamond = Summary estimate"
    ) +
    
    theme_clean +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  # Save
  filename_base <- clean_name(outcome_name)
  
  ggsave(
    filename = file.path(plots_dir, paste0("forest_", filename_base, ".png")),
    plot = p, width = 10, height = max(6, nrow(snp_data) * 0.3), dpi = 300, bg = "white"
  )
  
  ggsave(
    filename = file.path(plots_dir, paste0("forest_", filename_base, ".pdf")),
    plot = p, width = 10, height = max(6, nrow(snp_data) * 0.3), bg = "white"
  )
  
  message("💾 Saved forest plot: ", outcome_name)
  
  # Save individual outcome table with SNP details
  snp_table <- snp_data %>%
    filter(SNP != "All") %>%  # Exclude summary row for individual tables
    select(
      SNP = SNP,
      "Effect Size (β)" = b,
      "Standard Error" = se,
      "Odds Ratio" = OR,
      "OR Lower CI" = OR_lower,
      "OR Upper CI" = OR_upper,
      "P-value" = pval
    ) %>%
    mutate(
      "Effect Size (β)" = round(`Effect Size (β)`, 4),
      "Standard Error" = round(`Standard Error`, 4),
      "Odds Ratio" = round(`Odds Ratio`, 3),
      "OR Lower CI" = round(`OR Lower CI`, 3),
      "OR Upper CI" = round(`OR Upper CI`, 3),
      "P-value" = case_when(
        `P-value` < 0.001 ~ "< 0.001",
        `P-value` < 0.01 ~ sprintf("%.3f", `P-value`),
        TRUE ~ sprintf("%.2f", `P-value`)
      )
    )
  
  # Save individual SNP table
  filename_base <- clean_name(outcome_name)
  write.csv(snp_table, file.path(tables_dir, paste0("snp_details_", filename_base, ".csv")), 
            row.names = FALSE)
  message("💾 Saved SNP table: snp_details_", filename_base, ".csv")
}

################################################################################
# FINAL SUMMARY
################################################################################

message("\n🎉 MR visualization with organized output complete!")
message("📁 Results saved in organized structure:")
message("   📊 ", results_base_dir, "/")
message("     📈 plots/ - All visualization files")
message("     📋 tables/ - All summary tables")

message("\n📈 Generated visualizations:")
message("   • Scatter plots (", length(all_outcomes), " outcomes)")
message("   • Individual forest plots (", length(all_outcomes), " outcomes)")  
message("   • Summary forest plot (all outcomes)")
message("   • Bonferroni significant: ", length(bonferroni_significant))
message("   • Nominally significant: ", length(nominally_significant))

message("\n📋 Generated tables:")
message("   • Main results summary")
message("   • All methods comparison") 
message("   • Significance summary")
message("   • Individual SNP details (", length(all_outcomes), " outcomes)")

# List created files
plot_files <- list.files(plots_dir, pattern = "\\.(png|pdf)$")
table_files <- list.files(tables_dir, pattern = "\\.csv$")

message("\n📊 File summary:")
message("   📈 ", length(plot_files), " plot files (PNG + PDF)")
message("   📋 ", length(table_files), " CSV tables")

message("\n📋 Table files created:")
for(file in table_files) {
  message("   • ", file)
}

message("\n📄 All plots available in PNG (300 DPI) and PDF formats")
message("✨ Publication-ready with clean aesthetics and organized structure")
message("\n✅ Ready for manuscript preparation!")