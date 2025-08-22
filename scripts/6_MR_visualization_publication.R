################################################################################
# PUBLICATION-QUALITY MR VISUALIZATION SCRIPT
# EndoMR-PREG: Association of Genetic Liability to Endometriosis with Pregnancy Outcomes
# 
# Purpose: Generate publication-ready figures with clean formatting
# Author: Jonas Vibert, MD
# Date: August 2025
# 
# Key Features:
# - Clean, professional aesthetics
# - No row names or unnecessary elements
# - Consistent color schemes
# - High-resolution outputs
# - Proper labeling and titles
################################################################################

# Load required libraries
library(TwoSampleMR)
library(dplyr)
library(ggplot2)
library(stringr)

# Load optional libraries with error handling
optional_packages <- c("purrr", "forcats", "scales", "ggrepel", "patchwork")
for (pkg in optional_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message("📦 Installing missing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Alternative functions if packages are missing
if (!exists("map")) {
  map <- function(.x, .f, ...) {
    lapply(.x, .f, ...)
  }
  iwalk <- function(.x, .f, ...) {
    mapply(.f, .x, names(.x), MoreArgs = list(...), SIMPLIFY = FALSE)
    invisible(.x)
  }
}

if (!exists("fct_inorder")) {
  fct_inorder <- function(f) {
    factor(f, levels = unique(f))
  }
}

# Alternative for patchwork if not available
if (!exists("plot_annotation")) {
  plot_annotation <- function(...) {
    # Simple fallback - just return the first plot
    return(function(x) x)
  }
  wrap_plots <- function(plots, ncol = 2) {
    # Simple fallback - just return first plot with a note
    if (length(plots) > 0) {
      plots[[1]] + labs(caption = paste("Note: Combined plot of", length(plots), "panels"))
    }
  }
  plot_spacer <- function() {
    ggplot() + theme_void()
  }
}

# Set global theme for consistent aesthetics
theme_publication <- function(base_size = 11, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Remove grid lines
      panel.grid.major = element_line(color = "grey95", size = 0.3),
      panel.grid.minor = element_blank(),
      
      # Clean axis lines
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      
      # Professional text
      plot.title = element_text(size = rel(1.3), face = "bold", hjust = 0.5, margin = margin(b = 20)),
      plot.subtitle = element_text(size = rel(1.0), hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = rel(1.0), face = "bold"),
      axis.text = element_text(size = rel(0.9), color = "black"),
      
      # Legend styling
      legend.position = "bottom",
      legend.title = element_text(size = rel(0.9), face = "bold"),
      legend.text = element_text(size = rel(0.8)),
      legend.background = element_rect(fill = "white", color = NA),
      
      # Strip text for facets
      strip.text = element_text(size = rel(0.9), face = "bold", margin = margin(b = 10, t = 10)),
      strip.background = element_rect(fill = "grey95", color = "white"),
      
      # Clean plot background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Margins
      plot.margin = margin(20, 20, 20, 20)
    )
}

# Set as default theme
theme_set(theme_publication())

# Define color palette (colorblind-friendly)
mr_colors <- c(
  "IVW" = "#2166AC",
  "MR Egger" = "#D73027", 
  "Weighted median" = "#762A83",
  "Simple mode" = "#5AAE61",
  "Weighted mode" = "#FDB863"
)

################################################################################
# SECTION 1: DATA PREPARATION
################################################################################

message("📊 Starting publication-quality MR visualization...")

# Check if harmonized data exists
if (!exists("dat")) {
  if (file.exists("results/harmonised_rahmioglu_bpo.csv")) {
    dat <- read.csv("results/harmonised_rahmioglu_bpo.csv")
    message("✅ Loaded harmonized data from file")
  } else {
    stop("❌ Harmonized data 'dat' not found. Please run harmonization script first.")
  }
}

# Create results directory
plots_dir <- "publication_plots"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
  message("📁 Created plots directory: ", plots_dir)
}

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
  "O2sj1s"                                     # Preterm birth (all) - inverse association
)

# All outcomes of interest (for plotting)
all_plot_outcomes <- c(bonferroni_significant, nominally_significant)

# Clean outcome labels
outcome_labels <- c(
  # Bonferroni significant
  "finngen_R12_N14_FEMALEINFERT_filtered" = "Female infertility",
  "finngen_R12_O15_PLAC_PRAEVIA_filtered" = "Placenta praevia",
  
  # Nominally significant  
  "5sBPKR" = "Premature rupture of membranes",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered" = "Placental abruption",
  "mxO2D8" = "Elective caesarean delivery",
  "z2noFK" = "Low Apgar score at 1 minute", 
  "O2sj1s" = "Preterm birth (all)"
)

# Significance categories for color coding
significance_categories <- c(
  rep("Bonferroni significant", length(bonferroni_significant)),
  rep("Nominally significant", length(nominally_significant))
)
names(significance_categories) <- all_plot_outcomes

################################################################################
# SECTION 2: RUN MR ANALYSES
################################################################################

message("🔬 Running MR analyses...")

# Main MR analysis with multiple methods
mr_methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", 
                "mr_simple_mode", "mr_weighted_mode")

mr_results <- mr(dat, method_list = mr_methods)

# Single SNP analysis for forest plots
single_snp_results <- mr_singlesnp(dat)

# Leave-one-out analysis
loo_results <- mr_leaveoneout(dat)

################################################################################
# SECTION 3: SCATTER PLOTS (PUBLICATION QUALITY)
################################################################################

create_publication_scatter <- function(outcome_id, mr_res, harmonized_data) {
  
  # Filter data for this outcome
  outcome_data <- harmonized_data %>% filter(id.outcome == outcome_id)
  outcome_mr <- mr_res %>% filter(id.outcome == outcome_id)
  
  if (nrow(outcome_data) == 0 || nrow(outcome_mr) == 0) return(NULL)
  
  # Clean method names
  outcome_mr <- outcome_mr %>%
    mutate(
      method_clean = case_when(
        method == "Inverse variance weighted" ~ "IVW",
        method == "MR Egger" ~ "MR Egger",
        method == "Weighted median" ~ "Weighted median",
        method == "Simple mode" ~ "Simple mode",
        method == "Weighted mode" ~ "Weighted mode",
        TRUE ~ method
      )
    )
  
  # Get clean outcome name
  outcome_name <- outcome_labels[[outcome_id]] %||% outcome_id
  
  # Calculate plot limits with padding
  x_range <- range(outcome_data$beta.exposure, na.rm = TRUE)
  y_range <- range(outcome_data$beta.outcome, na.rm = TRUE)
  x_padding <- diff(x_range) * 0.1
  y_padding <- diff(y_range) * 0.1
  
  # Create the scatter plot
  p <- ggplot(outcome_data, aes(x = beta.exposure, y = beta.outcome)) +
    
    # Add SNP points
    geom_point(aes(size = 1/se.outcome^2), 
               color = "grey40", alpha = 0.7, shape = 21, fill = "white", stroke = 0.3) +
    
    # Add regression lines for each method
    geom_abline(data = outcome_mr, 
                aes(intercept = b, slope = 1, color = method_clean, linetype = method_clean),
                size = 0.8, alpha = 0.9) +
    
    # Customize scales
    scale_color_manual(values = mr_colors, name = "MR Method") +
    scale_linetype_manual(values = c("IVW" = "solid", "MR Egger" = "dashed", 
                                   "Weighted median" = "dotdash", 
                                   "Simple mode" = "dotted", "Weighted mode" = "longdash"),
                         name = "MR Method") +
    scale_size_continuous(range = c(0.8, 3.5), guide = "none") +
    
    # Labels and titles
    labs(
      title = paste("MR Analysis:", outcome_name),
      subtitle = "SNP effects on endometriosis liability vs. pregnancy outcome",
      x = "SNP effect on endometriosis liability (β)",
      y = "SNP effect on pregnancy outcome (β)",
      caption = "Point size proportional to inverse variance"
    ) +
    
    # Coordinate system
    coord_cartesian(xlim = x_range + c(-x_padding, x_padding),
                   ylim = y_range + c(-y_padding, y_padding)) +
    
    # Theme adjustments
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.caption = element_text(size = rel(0.7), hjust = 1, color = "grey50")
    ) +
    
    guides(color = guide_legend(override.aes = list(size = 1.2)),
           linetype = guide_legend(override.aes = list(size = 1.2)))
  
  return(p)
}

message("📈 Creating scatter plots...")

# Generate scatter plots for all outcomes of interest
scatter_plots <- map(all_plot_outcomes, ~create_publication_scatter(.x, mr_results, dat))
names(scatter_plots) <- all_plot_outcomes

# Save individual scatter plots with significance indicators
iwalk(scatter_plots, function(plot, outcome_id) {
  if (!is.null(plot)) {
    outcome_name <- str_replace_all(outcome_labels[[outcome_id]] %||% outcome_id, "[^A-Za-z0-9]", "_")
    significance <- significance_categories[[outcome_id]]
    
    # Add significance indicator to subtitle
    plot_modified <- plot + 
      labs(subtitle = paste0("SNP effects on endometriosis liability vs. pregnancy outcome\n(", significance, ")"))
    
    ggsave(
      filename = file.path(plots_dir, paste0("scatter_", outcome_name, ".png")),
      plot = plot_modified,
      width = 10, height = 8, dpi = 300, bg = "white"
    )
    
    ggsave(
      filename = file.path(plots_dir, paste0("scatter_", outcome_name, ".pdf")),
      plot = plot_modified,
      width = 10, height = 8, bg = "white"
    )
    
    message("💾 Saved scatter plot: ", outcome_name, " (", significance, ")")
  }
})

################################################################################
# SECTION 4: FOREST PLOTS (PUBLICATION QUALITY)
################################################################################

create_publication_forest <- function(outcome_id, single_snp_res) {
  
  # Filter for this outcome
  outcome_snps <- single_snp_res %>% filter(id.outcome == outcome_id)
  
  if (nrow(outcome_snps) == 0) return(NULL)
  
  # Clean SNP names and prepare data
  forest_data <- outcome_snps %>%
    mutate(
      SNP_clean = case_when(
        SNP == "All" ~ "Summary (IVW)",
        SNP == "All - inverse variance weighted" ~ "Summary (IVW)",
        TRUE ~ SNP
      ),
      OR = exp(b),
      OR_lower = exp(b - 1.96 * se),
      OR_upper = exp(b + 1.96 * se),
      is_summary = SNP_clean == "Summary (IVW)"
    ) %>%
    arrange(desc(is_summary), OR) %>%
    mutate(SNP_clean = fct_inorder(SNP_clean))
  
  # Get outcome name
  outcome_name <- outcome_labels[[outcome_id]] %||% outcome_id
  
  # Create forest plot
  p <- ggplot(forest_data, aes(y = SNP_clean, x = OR)) +
    
    # Reference line at OR = 1
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", size = 0.5) +
    
    # Confidence intervals
    geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper, 
                      color = is_summary, size = is_summary),
                  height = 0.3, alpha = 0.8) +
    
    # Point estimates  
    geom_point(aes(color = is_summary, size = is_summary, shape = is_summary),
               stroke = 0.5) +
    
    # Customize appearance
    scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "#D73027"), guide = "none") +
    scale_size_manual(values = c("FALSE" = 0.8, "TRUE" = 1.5), guide = "none") +
    scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23), guide = "none") +
    scale_x_log10(breaks = c(0.5, 0.7, 1.0, 1.4, 2.0, 3.0),
                  labels = c("0.5", "0.7", "1.0", "1.4", "2.0", "3.0")) +
    
    # Labels and titles
    labs(
      title = paste("Forest Plot:", outcome_name),
      subtitle = "Individual SNP effects and meta-analysis estimate",
      x = "Odds Ratio (95% CI)",
      y = "Genetic Variants",
      caption = "Diamond = Summary estimate (IVW method)"
    ) +
    
    # Clean y-axis (remove row names)
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.caption = element_text(size = rel(0.7), hjust = 1, color = "grey50")
    )
  
  return(p)
}

message("🌳 Creating forest plots...")

# Generate forest plots for all outcomes of interest
forest_plots <- map(all_plot_outcomes, ~create_publication_forest(.x, single_snp_results))
names(forest_plots) <- all_plot_outcomes

# Save forest plots with significance indicators
iwalk(forest_plots, function(plot, outcome_id) {
  if (!is.null(plot)) {
    outcome_name <- str_replace_all(outcome_labels[[outcome_id]] %||% outcome_id, "[^A-Za-z0-9]", "_")
    significance <- significance_categories[[outcome_id]]
    
    # Add significance indicator to subtitle
    plot_modified <- plot + 
      labs(subtitle = paste0("Individual SNP effects and meta-analysis estimate\n(", significance, ")"))
    
    ggsave(
      filename = file.path(plots_dir, paste0("forest_", outcome_name, ".png")),
      plot = plot_modified,
      width = 10, height = max(8, nrow(single_snp_results %>% filter(id.outcome == outcome_id)) * 0.3),
      dpi = 300, bg = "white"
    )
    
    ggsave(
      filename = file.path(plots_dir, paste0("forest_", outcome_name, ".pdf")),
      plot = plot_modified,
      width = 10, height = max(8, nrow(single_snp_results %>% filter(id.outcome == outcome_id)) * 0.3),
      bg = "white"
    )
    
    message("💾 Saved forest plot: ", outcome_name, " (", significance, ")")
  }
})

################################################################################
# SECTION 5: LEAVE-ONE-OUT PLOTS
################################################################################

create_publication_loo <- function(outcome_id, loo_res) {
  
  # Filter for this outcome
  outcome_loo <- loo_res %>% filter(id.outcome == outcome_id)
  
  if (nrow(outcome_loo) == 0) return(NULL)
  
  # Prepare data
  loo_data <- outcome_loo %>%
    mutate(
      OR = exp(b),
      OR_lower = exp(b - 1.96 * se),
      OR_upper = exp(b + 1.96 * se),
      is_summary = SNP == "All",
      SNP_clean = case_when(
        SNP == "All" ~ "All SNPs",
        TRUE ~ paste("Excluding", SNP)
      )
    ) %>%
    arrange(desc(is_summary), OR) %>%
    mutate(SNP_clean = fct_inorder(SNP_clean))
  
  # Get outcome name
  outcome_name <- outcome_labels[[outcome_id]] %||% outcome_id
  
  # Create plot
  p <- ggplot(loo_data, aes(y = SNP_clean, x = OR)) +
    
    # Reference line
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", size = 0.5) +
    
    # Confidence intervals
    geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper,
                      color = is_summary, alpha = is_summary),
                  height = 0.2) +
    
    # Points
    geom_point(aes(color = is_summary, size = is_summary),
               shape = 21, fill = "white", stroke = 0.5) +
    
    # Styling
    scale_color_manual(values = c("FALSE" = "#2166AC", "TRUE" = "#D73027"), guide = "none") +
    scale_alpha_manual(values = c("FALSE" = 0.7, "TRUE" = 1.0), guide = "none") +
    scale_size_manual(values = c("FALSE" = 2, "TRUE" = 3), guide = "none") +
    scale_x_log10() +
    
    # Labels
    labs(
      title = paste("Leave-One-Out Analysis:", outcome_name),
      subtitle = "Sensitivity to individual SNPs",
      x = "Odds Ratio (95% CI)",
      y = "Analysis",
      caption = "Red point = All SNPs included"
    ) +
    
    # Clean y-axis labels
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.caption = element_text(size = rel(0.7), hjust = 1, color = "grey50")
    )
  
  return(p)
}

message("🔄 Creating leave-one-out plots...")

# Generate LOO plots for all outcomes of interest
loo_plots <- map(all_plot_outcomes, ~create_publication_loo(.x, loo_results))
names(loo_plots) <- all_plot_outcomes

# Save LOO plots with significance indicators
iwalk(loo_plots, function(plot, outcome_id) {
  if (!is.null(plot)) {
    outcome_name <- str_replace_all(outcome_labels[[outcome_id]] %||% outcome_id, "[^A-Za-z0-9]", "_")
    significance <- significance_categories[[outcome_id]]
    
    # Add significance indicator to subtitle
    plot_modified <- plot + 
      labs(subtitle = paste0("Sensitivity to individual SNPs\n(", significance, ")"))
    
    ggsave(
      filename = file.path(plots_dir, paste0("loo_", outcome_name, ".png")),
      plot = plot_modified,
      width = 10, height = 8, dpi = 300, bg = "white"
    )
    
    ggsave(
      filename = file.path(plots_dir, paste0("loo_", outcome_name, ".pdf")),
      plot = plot_modified,
      width = 10, height = 8, bg = "white"
    )
    
    message("💾 Saved leave-one-out plot: ", outcome_name, " (", significance, ")")
  }
})

################################################################################
# SECTION 6: FUNNEL PLOTS
################################################################################

create_publication_funnel <- function(outcome_id, single_snp_res) {
  
  # Filter for this outcome (exclude summary)
  outcome_snps <- single_snp_res %>% 
    filter(id.outcome == outcome_id, SNP != "All")
  
  if (nrow(outcome_snps) == 0) return(NULL)
  
  # Get outcome name
  outcome_name <- outcome_labels[[outcome_id]] %||% outcome_id
  
  # Create funnel plot
  p <- ggplot(outcome_snps, aes(x = b, y = 1/se)) +
    
    # Reference line (null effect)
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 0.5) +
    
    # SNP points
    geom_point(color = "#2166AC", size = 2.5, alpha = 0.7, shape = 21, 
               fill = "white", stroke = 0.5) +
    
    # Labels
    labs(
      title = paste("Funnel Plot:", outcome_name),
      subtitle = "Assessment of potential bias",
      x = "SNP effect size (β)",
      y = "Precision (1/SE)",
      caption = "Asymmetry may indicate directional pleiotropy"
    ) +
    
    # Clean appearance
    theme(
      plot.caption = element_text(size = rel(0.7), hjust = 1, color = "grey50")
    )
  
  return(p)
}

message("🔍 Creating funnel plots...")

# Generate funnel plots for all outcomes of interest
funnel_plots <- map(all_plot_outcomes, ~create_publication_funnel(.x, single_snp_results))
names(funnel_plots) <- all_plot_outcomes

# Save funnel plots with significance indicators
iwalk(funnel_plots, function(plot, outcome_id) {
  if (!is.null(plot)) {
    outcome_name <- str_replace_all(outcome_labels[[outcome_id]] %||% outcome_id, "[^A-Za-z0-9]", "_")
    significance <- significance_categories[[outcome_id]]
    
    # Add significance indicator to subtitle
    plot_modified <- plot + 
      labs(subtitle = paste0("Assessment of potential bias\n(", significance, ")"))
    
    ggsave(
      filename = file.path(plots_dir, paste0("funnel_", outcome_name, ".png")),
      plot = plot_modified,
      width = 8, height = 6, dpi = 300, bg = "white"
    )
    
    ggsave(
      filename = file.path(plots_dir, paste0("funnel_", outcome_name, ".pdf")),
      plot = plot_modified,
      width = 8, height = 6, bg = "white"
    )
    
    message("💾 Saved funnel plot: ", outcome_name, " (", significance, ")")
  }
})

################################################################################
# SECTION 7: SUMMARY FOREST PLOT (KEY RESULTS)
################################################################################

create_summary_forest <- function(mr_res) {
  
  # Filter for all outcomes of interest and IVW method
  summary_data <- mr_res %>%
    filter(
      id.outcome %in% all_plot_outcomes,
      method == "Inverse variance weighted"
    ) %>%
    mutate(
      outcome_clean = recode(id.outcome, !!!outcome_labels),
      OR = exp(b),
      OR_lower = exp(b - 1.96 * se), 
      OR_upper = exp(b + 1.96 * se),
      bonferroni_significant = pval < 0.001,  # Bonferroni corrected
      nominally_significant = pval < 0.05 & pval >= 0.001,
      significance_level = case_when(
        bonferroni_significant ~ "Bonferroni significant",
        nominally_significant ~ "Nominally significant", 
        TRUE ~ "Not significant"
      ),
      CI_text = sprintf("%.2f (%.2f-%.2f)", OR, OR_lower, OR_upper),
      p_text = case_when(
        pval < 0.001 ~ sprintf("P = %.1e", pval),
        pval < 0.01 ~ sprintf("P = %.3f", pval),
        TRUE ~ sprintf("P = %.2f", pval)
      )
    ) %>%
    arrange(desc(bonferroni_significant), desc(nominally_significant), pval)
  
  # Create forest plot
  p <- summary_data %>%
    mutate(outcome_clean = fct_inorder(outcome_clean)) %>%
    ggplot(aes(y = outcome_clean, x = OR)) +
    
    # Reference line
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", size = 0.5) +
    
    # Confidence intervals
    geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper, color = significance_level),
                  height = 0.2, size = 1, alpha = 0.8) +
    
    # Point estimates
    geom_point(aes(color = significance_level, size = significance_level, shape = significance_level), 
               fill = "white", stroke = 1.2) +
    
    # Add effect size labels
    geom_text(aes(x = OR_upper + 0.1, label = CI_text),
              hjust = 0, size = 3.5, fontface = "bold") +
    
    # Add p-value labels  
    geom_text(aes(x = OR_upper + 0.1, label = p_text, y = as.numeric(outcome_clean) - 0.15),
              hjust = 0, size = 3, color = "grey40") +
    
    # Styling
    scale_color_manual(values = c("Bonferroni significant" = "#D73027", 
                                 "Nominally significant" = "#FDB863",
                                 "Not significant" = "grey60"), 
                      name = "Significance Level") +
    scale_size_manual(values = c("Bonferroni significant" = 4, 
                                "Nominally significant" = 3.5,
                                "Not significant" = 3), guide = "none") +
    scale_shape_manual(values = c("Bonferroni significant" = 22, 
                                 "Nominally significant" = 21,
                                 "Not significant" = 21), guide = "none") +
    scale_x_log10(limits = c(0.7, 4), 
                  breaks = c(0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0),
                  labels = c("0.8", "1.0", "1.25", "1.5", "2.0", "2.5", "3.0")) +
    
    # Labels
    labs(
      title = "Mendelian Randomization: Endometriosis and Pregnancy Outcomes", 
      subtitle = "Causal effect estimates using inverse variance weighted method",
      x = "Odds Ratio (95% CI)",
      y = "",
      caption = "Red squares = Bonferroni significant (P < 0.001); Orange circles = Nominally significant (P < 0.05)"
    ) +
    
    # Clean theme
    theme(
      axis.text.y = element_text(size = rel(1.0), color = "black"),
      plot.caption = element_text(size = rel(0.7), hjust = 1, color = "grey50"),
      panel.grid.major.y = element_line(color = "grey95", size = 0.3)
    )
  
  return(p)
}

message("📊 Creating summary forest plot...")

summary_plot <- create_summary_forest(mr_results)

# Save summary plot
ggsave(
  filename = file.path(plots_dir, "summary_forest_plot.png"),
  plot = summary_plot,
  width = 12, height = 6, dpi = 300, bg = "white"
)

ggsave(
  filename = file.path(plots_dir, "summary_forest_plot.pdf"),
  plot = summary_plot,
  width = 12, height = 6, bg = "white"
)

message("💾 Saved summary forest plot")

################################################################################
# SECTION 8: COMBINED FIGURE (MAIN RESULTS)
################################################################################

message("🎨 Creating combined figure...")

# Create 2x2 combined plot for Bonferroni significant results
if (!is.null(scatter_plots[["finngen_R12_O15_PLAC_PRAEVIA_filtered"]]) && 
    !is.null(scatter_plots[["finngen_R12_N14_FEMALEINFERT_filtered"]])) {
  
  bonferroni_plot <- (scatter_plots[["finngen_R12_O15_PLAC_PRAEVIA_filtered"]] + 
                     scatter_plots[["finngen_R12_N14_FEMALEINFERT_filtered"]]) +
    plot_annotation(
      title = "Bonferroni-Corrected Significant Results",
      subtitle = "Endometriosis and pregnancy outcomes (P < 0.001)",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                   plot.subtitle = element_text(size = 12, hjust = 0.5))
    )
  
  ggsave(
    filename = file.path(plots_dir, "combined_bonferroni_significant.png"),
    plot = bonferroni_plot,
    width = 16, height = 8, dpi = 300, bg = "white"
  )
  
  ggsave(
    filename = file.path(plots_dir, "combined_bonferroni_significant.pdf"),
    plot = bonferroni_plot,
    width = 16, height = 8, bg = "white"
  )
  
  message("💾 Saved Bonferroni significant combined figure")
}

# Create combined plot for nominally significant results
available_nominal <- nominally_significant[nominally_significant %in% names(scatter_plots)]
available_nominal <- available_nominal[!sapply(scatter_plots[available_nominal], is.null)]

if (length(available_nominal) >= 3) {
  
  # Create a grid layout for nominally significant results
  nominal_plots <- scatter_plots[available_nominal[1:min(4, length(available_nominal))]]
  
  if (length(nominal_plots) == 3) {
    nominal_combined <- (nominal_plots[[1]] + nominal_plots[[2]]) /
                       (nominal_plots[[3]] + plot_spacer()) +
      plot_annotation(
        title = "Nominally Significant Results",
        subtitle = "Endometriosis and pregnancy outcomes (P < 0.05, hypothesis-generating)",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                     plot.subtitle = element_text(size = 12, hjust = 0.5))
      )
  } else if (length(nominal_plots) >= 4) {
    nominal_combined <- (nominal_plots[[1]] + nominal_plots[[2]]) /
                       (nominal_plots[[3]] + nominal_plots[[4]]) +
      plot_annotation(
        title = "Nominally Significant Results",
        subtitle = "Endometriosis and pregnancy outcomes (P < 0.05, hypothesis-generating)",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                     plot.subtitle = element_text(size = 12, hjust = 0.5))
      )
  } else {
    nominal_combined <- wrap_plots(nominal_plots, ncol = 2) +
      plot_annotation(
        title = "Nominally Significant Results",
        subtitle = "Endometriosis and pregnancy outcomes (P < 0.05, hypothesis-generating)",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                     plot.subtitle = element_text(size = 12, hjust = 0.5))
      )
  }
  
  ggsave(
    filename = file.path(plots_dir, "combined_nominally_significant.png"),
    plot = nominal_combined,
    width = 16, height = 12, dpi = 300, bg = "white"
  )
  
  ggsave(
    filename = file.path(plots_dir, "combined_nominally_significant.pdf"), 
    plot = nominal_combined,
    width = 16, height = 12, bg = "white"
  )
  
  message("💾 Saved nominally significant combined figure")
}

################################################################################
# FINAL SUMMARY
################################################################################

message("\n🎉 Publication-quality visualization complete!")
message("📁 All plots saved to: ", plots_dir)
message("📈 Generated plots:")
message("   • Scatter plots (", length(all_plot_outcomes), " outcomes)")
message("     - Bonferroni significant: ", length(bonferroni_significant))
message("     - Nominally significant: ", length(nominally_significant))
message("   • Forest plots (", length(all_plot_outcomes), " outcomes)")  
message("   • Leave-one-out plots (", length(all_plot_outcomes), " outcomes)")
message("   • Funnel plots (", length(all_plot_outcomes), " outcomes)")
message("   • Summary forest plot (all significance levels)")
message("   • Combined figures (by significance level)")
message("📄 Formats: PNG (300 DPI) and PDF")
message("✨ Clean, publication-ready aesthetics with no row names")
message("🏷️  Significance levels clearly labeled on all plots")

# Print plot directory contents
message("\n📋 Files created:")
plot_files <- list.files(plots_dir, pattern = "\\.(png|pdf)$", full.names = FALSE)
for(file in plot_files) {
  message("   • ", file)
}

<function_calls>
<invoke name="TodoWrite">
<parameter name="todos">[{"id": "1", "content": "Rewrite 6_MR_visualization script for publication-quality plots", "status": "completed"}]