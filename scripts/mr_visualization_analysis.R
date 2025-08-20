################################################################################
# MENDELIAN RANDOMIZATION VISUALIZATION ANALYSIS
# Endometriosis and Pregnancy Outcomes Study
# 
# Purpose: Generate comprehensive MR plots for significant outcomes
# Author: [Your Name]
# Date: August 2025
# 
# Sections:
# 1. Data Preparation and Methods Configuration
# 2. Scatter Plots (IVW + MR-Egger + Weighted Median)
# 3. Forest Plots for Key Outcomes
# 4. Leave-One-Out Analysis (MR-Egger)
# 5. Funnel Plots for Publication Bias
# 6. One-to-Many Forest Plot (All Outcomes)
################################################################################

# Load required libraries
library(TwoSampleMR)
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)

# Optional: SVG export capability
if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite")
}

################################################################################
# SECTION 1: DATA PREPARATION AND CONFIGURATION
################################################################################

# Load harmonized data if not already in memory
# Uncomment and modify path as needed:
# dat <- read.csv("path/to/harmonized_data.csv")
# dat <- data.table::fread("path/to/harmonized_data.csv") %>% as.data.frame()

# Define MR methods for analysis
mr_methods <- c(
  "mr_ivw",
  "mr_egger_regression", 
  "mr_weighted_median",
  "mr_weighted_mode",
  "mr_simple_mode"
)

# Key significant outcomes for analysis
significant_outcomes <- c(
  "finngen_R12_N14_FEMALEINFERT_filtered",     # Female infertility
  "5sBPKR",                                    # Premature rupture of membranes
  "finngen_R12_O15_PLAC_PRAEVIA_filtered",     # Placenta praevia
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered" # Placental abruption
)

# Human-readable outcome labels
outcome_labels <- c(
  "finngen_R12_N14_FEMALEINFERT_filtered"     = "Female infertility (FinnGen)",
  "5sBPKR"                                    = "Premature rupture of membranes",
  "finngen_R12_O15_PLAC_PRAEVIA_filtered"     = "Placenta praevia (FinnGen)",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered" = "Placental abruption (FinnGen)"
)

# Run MR analysis with all methods
if (!exists("res_all")) {
  message("Running MR analysis with multiple methods...")
  res_all <- mr(dat, method_list = mr_methods)
}

# Generate single SNP results for forest and funnel plots
if (!exists("res_single")) {
  message("Generating single SNP MR results...")
  res_single <- mr_singlesnp(dat)
}

################################################################################
# SECTION 2: SCATTER PLOTS
################################################################################

message("Generating scatter plots...")

# Create output directory
scatter_dir <- "mr_scatter_plots"
if (!dir.exists(scatter_dir)) dir.create(scatter_dir, recursive = TRUE)

# Generate scatter plots for each significant outcome
scatter_plots <- map(significant_outcomes, function(outcome_id) {
  
  # Filter data for current outcome
  dat_outcome <- dat %>% filter(id.outcome == outcome_id)
  
  # Filter MR results for current outcome
  res_outcome <- res_all %>%
    filter(
      id.outcome == outcome_id,
      method %in% c(
        "Inverse variance weighted",
        "MR Egger", 
        "Weighted median",
        "Simple mode",
        "Weighted mode"
      )
    ) %>%
    arrange(method)
  
  # Skip if insufficient data
  if (nrow(res_outcome) == 0 || nrow(dat_outcome) == 0) {
    message("Insufficient data for scatter plot: ", outcome_id)
    return(NULL)
  }
  
  # Generate scatter plot
  plot_list <- mr_scatter_plot(res_outcome, dat_outcome)
  
  if (length(plot_list) == 0) {
    message("No scatter plot generated for: ", outcome_id)
    return(NULL)
  }
  
  # Customize plot
  plot_final <- plot_list[[1]] +
    ggtitle(paste0("MR Analysis: ", outcome_labels[[outcome_id]] %||% outcome_id)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 11)
    )
  
  # Save plots
  png_file <- file.path(scatter_dir, paste0(outcome_id, "_scatter.png"))
  pdf_file <- file.path(scatter_dir, paste0(outcome_id, "_scatter.pdf"))
  
  ggsave(png_file, plot_final, width = 8, height = 6, dpi = 300)
  ggsave(pdf_file, plot_final, width = 8, height = 6)
  
  message("Saved scatter plot: ", basename(png_file))
  return(plot_final)
})

names(scatter_plots) <- significant_outcomes

################################################################################
# SECTION 3: FOREST PLOTS FOR KEY OUTCOMES
################################################################################

message("Generating forest plots...")

# Create output directory
forest_dir <- "mr_forest_plots"
if (!dir.exists(forest_dir)) dir.create(forest_dir, recursive = TRUE)

# Function to generate and save forest plot
generate_forest_plot <- function(dat, outcome_id, title_text) {
  
  # Filter data for outcome
  dat_outcome <- dat %>% filter(id.outcome == outcome_id)
  
  if (nrow(dat_outcome) == 0) {
    message("No data for forest plot: ", outcome_id)
    return(NULL)
  }
  
  # Generate single SNP results for this outcome
  res_single_outcome <- mr_singlesnp(dat_outcome)
  
  if (nrow(res_single_outcome) == 0) {
    message("No single SNP results for: ", outcome_id)
    return(NULL)
  }
  
  # Create forest plot
  plot_list <- mr_forest_plot(res_single_outcome)
  
  if (length(plot_list) == 0) {
    message("No forest plot generated for: ", outcome_id)
    return(NULL)
  }
  
  # Customize plot
  plot_final <- plot_list[[1]] +
    ggtitle(title_text) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12)
    )
  
  # Save plots
  png_file <- file.path(forest_dir, paste0(outcome_id, "_forest.png"))
  pdf_file <- file.path(forest_dir, paste0(outcome_id, "_forest.pdf"))
  
  ggsave(png_file, plot_final, width = 8, height = 10, dpi = 300)
  ggsave(pdf_file, plot_final, width = 8, height = 10)
  
  message("Saved forest plot: ", basename(png_file))
  return(plot_final)
}

# Generate forest plots for all significant outcomes
forest_plots <- map(significant_outcomes, 
                   ~generate_forest_plot(dat, .x, outcome_labels[[.x]]))
names(forest_plots) <- significant_outcomes

################################################################################
# SECTION 4: LEAVE-ONE-OUT ANALYSIS (MR-EGGER)
################################################################################

message("Generating leave-one-out plots...")

# Create output directories
loo_dir <- "mr_leaveoneout_plots_egger"
loo_tables_dir <- file.path(loo_dir, "tables")
if (!dir.exists(loo_dir)) dir.create(loo_dir, recursive = TRUE)
if (!dir.exists(loo_tables_dir)) dir.create(loo_tables_dir, recursive = TRUE)

# Generate leave-one-out plots
loo_plots <- list()

for (outcome_id in significant_outcomes) {
  
  # Filter data for outcome
  dat_outcome <- dat %>% filter(id.outcome == outcome_id)
  
  # Check minimum SNP requirement for MR-Egger
  if (nrow(dat_outcome) < 3) {
    message("Insufficient SNPs for leave-one-out analysis: ", outcome_id, 
            " (", nrow(dat_outcome), " SNPs)")
    next
  }
  
  message("Running MR-Egger leave-one-out for: ", outcome_id)
  
  # Perform leave-one-out analysis
  loo_results <- tryCatch({
    mr_leaveoneout(dat_outcome, method = mr_egger_regression)
  }, error = function(e) {
    warning("Leave-one-out analysis failed for ", outcome_id, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(loo_results) || nrow(loo_results) == 0) {
    message("No leave-one-out results for: ", outcome_id)
    next
  }
  
  # Save numerical results
  csv_file <- file.path(loo_tables_dir, paste0(outcome_id, "_leaveoneout_egger.csv"))
  write.csv(loo_results, csv_file, row.names = FALSE)
  
  # Generate plot
  plot_list <- mr_leaveoneout_plot(loo_results)
  
  if (length(plot_list) == 0) {
    message("No leave-one-out plot generated for: ", outcome_id)
    next
  }
  
  # Customize plot
  plot_final <- plot_list[[1]] +
    ggtitle(paste0("Leave-One-Out Analysis (MR-Egger): ", outcome_labels[[outcome_id]])) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11)
    )
  
  loo_plots[[outcome_id]] <- plot_final
  
  # Save plots
  png_file <- file.path(loo_dir, paste0(outcome_id, "_leaveoneout_egger.png"))
  pdf_file <- file.path(loo_dir, paste0(outcome_id, "_leaveoneout_egger.pdf"))
  
  ggsave(png_file, plot_final, width = 8, height = 7, dpi = 300)
  ggsave(pdf_file, plot_final, width = 8, height = 7)
  
  message("Saved leave-one-out plot: ", basename(png_file))
}

################################################################################
# SECTION 5: FUNNEL PLOTS FOR PUBLICATION BIAS
################################################################################

message("Generating funnel plots...")

# Create output directory
funnel_dir <- "mr_funnel_plots"
if (!dir.exists(funnel_dir)) dir.create(funnel_dir, recursive = TRUE)

# Generate funnel plots
funnel_plots <- list()

for (outcome_id in significant_outcomes) {
  
  # Filter single SNP results for outcome
  res_outcome <- res_single %>% filter(id.outcome == outcome_id)
  
  if (nrow(res_outcome) == 0) {
    message("No single SNP data for funnel plot: ", outcome_id)
    next
  }
  
  # Generate funnel plot
  plot_list <- mr_funnel_plot(res_outcome)
  
  if (length(plot_list) == 0) {
    message("No funnel plot generated for: ", outcome_id)
    next
  }
  
  # Customize plot
  plot_final <- plot_list[[1]] +
    ggtitle(paste0("Funnel Plot: ", outcome_labels[[outcome_id]])) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12)
    )
  
  funnel_plots[[outcome_id]] <- plot_final
  
  # Save plots
  png_file <- file.path(funnel_dir, paste0(outcome_id, "_funnel.png"))
  pdf_file <- file.path(funnel_dir, paste0(outcome_id, "_funnel.pdf"))
  
  ggsave(png_file, plot_final, width = 7, height = 6, dpi = 300)
  ggsave(pdf_file, plot_final, width = 7, height = 6)
  
  message("Saved funnel plot: ", basename(png_file))
}

################################################################################
# SECTION 6: ONE-TO-MANY FOREST PLOT (COMPREHENSIVE ANALYSIS)
################################################################################

message("Generating comprehensive one-to-many forest plot...")

# Extended list of outcomes for comprehensive analysis
all_outcomes <- c(
  # Core significant outcomes
  "finngen_R12_N14_FEMALEINFERT_filtered",
  "5sBPKR", 
  "finngen_R12_O15_PLAC_PRAEVIA_filtered",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered",
  # Additional outcomes
  "oAOiBt", "Rwbyma", "zMUPe9", "vywWPW",  # Hypertensive/metabolic
  "O2sj1s", "9QjNDJ", "i3cwH7",             # Preterm/PROM
  "TTnKdM", "mxO2D8", "OYxJPL",             # Caesarean
  # Bleeding/PPH
  "Postpartum_hemorrhage_filtered",
  "Postpartum_hemorrhage_due_to_atony_filtered", 
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered",
  "Antepartum_bleeding_filtered",
  "Early_bleeding_with_any_outcome_filtered",
  "Early_bleeding_ending_in_live_birth_filtered",
  # Growth & size
  "zKkCNG", "ycgY62", "yoKaH0", "iS2FMp", "dhBaTu",
  # Neonatal outcomes
  "z2noFK", "dBcqKG", "wvATQb", "RJauGo", "6g03jM"
)

# Comprehensive outcome labels
comprehensive_labels <- c(
  # Core outcomes
  "finngen_R12_N14_FEMALEINFERT_filtered"     = "Female infertility (FinnGen)",
  "5sBPKR"                                    = "Premature rupture of membranes",
  "finngen_R12_O15_PLAC_PRAEVIA_filtered"     = "Placenta praevia (FinnGen)",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered" = "Placental abruption (FinnGen)",
  # Hypertensive/metabolic
  "oAOiBt" = "Hypertensive disorders of pregnancy",
  "Rwbyma" = "Preeclampsia",
  "zMUPe9" = "Gestational hypertension", 
  "vywWPW" = "Gestational diabetes",
  # Preterm/PROM
  "O2sj1s" = "Preterm birth (all)",
  "9QjNDJ" = "Preterm birth (subsample)",
  "i3cwH7" = "Very preterm birth",
  # Delivery mode
  "TTnKdM" = "Caesarean section (all)",
  "mxO2D8" = "Elective caesarean section",
  "OYxJPL" = "Emergency caesarean section",
  # Bleeding/PPH
  "Postpartum_hemorrhage_filtered"                          = "Postpartum haemorrhage (overall)",
  "Postpartum_hemorrhage_due_to_atony_filtered"             = "PPH - uterine atony",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered" = "PPH - retained placenta",
  "Antepartum_bleeding_filtered"                            = "Antepartum bleeding",
  "Early_bleeding_with_any_outcome_filtered"                = "Early bleeding (any outcome)",
  "Early_bleeding_ending_in_live_birth_filtered"            = "Early bleeding (live birth)",
  # Growth & size
  "zKkCNG" = "Low birthweight (<2500g)",
  "ycgY62" = "Large for gestational age",
  "yoKaH0" = "Small for gestational age", 
  "iS2FMp" = "Birthweight Z-score",
  "dhBaTu" = "High birthweight (>4000g)",
  # Neonatal
  "z2noFK" = "Low Apgar at 1 minute",
  "dBcqKG" = "Low Apgar at 5 minutes",
  "wvATQb" = "Apgar score at 1 minute",
  "RJauGo" = "Apgar score at 5 minutes", 
  "6g03jM" = "NICU admission"
)

# Clinical groupings for faceting
clinical_groups <- c(
  "finngen_R12_N14_FEMALEINFERT_filtered"     = "Fertility",
  "5sBPKR"                                    = "Preterm & PROM",
  "finngen_R12_O15_PLAC_PRAEVIA_filtered"     = "Placental disorders",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered" = "Placental disorders",
  "oAOiBt" = "Hypertensive disorders", "Rwbyma" = "Hypertensive disorders",
  "zMUPe9" = "Hypertensive disorders", "vywWPW" = "Metabolic",
  "O2sj1s" = "Preterm & PROM", "9QjNDJ" = "Preterm & PROM", "i3cwH7" = "Preterm & PROM",
  "TTnKdM" = "Delivery mode", "mxO2D8" = "Delivery mode", "OYxJPL" = "Delivery mode",
  "Postpartum_hemorrhage_filtered" = "Bleeding/PPH",
  "Postpartum_hemorrhage_due_to_atony_filtered" = "Bleeding/PPH",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered" = "Bleeding/PPH",
  "Antepartum_bleeding_filtered" = "Bleeding/PPH",
  "Early_bleeding_with_any_outcome_filtered" = "Bleeding/PPH",
  "Early_bleeding_ending_in_live_birth_filtered" = "Bleeding/PPH",
  "zKkCNG" = "Growth & size", "ycgY62" = "Growth & size", "yoKaH0" = "Growth & size",
  "iS2FMp" = "Growth & size", "dhBaTu" = "Growth & size",
  "z2noFK" = "Apgar & NICU", "dBcqKG" = "Apgar & NICU", "wvATQb" = "Apgar & NICU",
  "RJauGo" = "Apgar & NICU", "6g03jM" = "Apgar & NICU"
)

# Ensure MR results exist
if (!exists("res")) {
  if (!exists("dat")) stop("Need either 'res' or 'dat' in memory.")
  res <- mr(dat, method_list = mr_methods)
}

# Prepare results for forest plot
results_formatted <- res %>%
  filter(id.outcome %in% all_outcomes) %>%
  subset_on_method() %>%  # Selects IVW or Wald ratio automatically
  mutate(
    outcome_label = dplyr::recode(outcome, !!!comprehensive_labels),
    group = dplyr::recode(id.outcome, !!!clinical_groups),
    OR = exp(b),
    LCL = exp(b - 1.96 * se),
    UCL = exp(b + 1.96 * se),
    Effect_CI = sprintf("%.2f (%.2f-%.2f)", OR, LCL, UCL),
    P_value = formatC(pval, format = "e", digits = 2),
    Method = dplyr::recode(method,
                          "Inverse variance weighted" = "IVW",
                          "MR Egger" = "MR-Egger", 
                          "Wald ratio" = "Wald",
                          .default = method)
  )

# Format labels and set group order
results_formatted$outcome_label <- str_wrap(results_formatted$outcome_label, width = 40)

group_order <- c("Fertility", "Placental disorders", "Preterm & PROM", "Delivery mode",
                "Bleeding/PPH", "Hypertensive disorders", "Metabolic", 
                "Growth & size", "Apgar & NICU")
results_formatted$group <- factor(results_formatted$group, levels = group_order)

# Sort by group and effect size
results_formatted <- results_formatted %>% 
  arrange(group, desc(abs(b)))

# Normalize point sizes for visualization
results_formatted <- results_formatted %>% 
  mutate(weight_plot = 0.6 * (1/se) / max(1/se, na.rm = TRUE))

# Set axis limits (OR scale)
or_lower <- max(0.30, min(results_formatted$LCL, na.rm = TRUE) * 0.9)
or_upper <- min(3.00, max(results_formatted$UCL, na.rm = TRUE) * 1.1)

# Create output directory
forest_comprehensive_dir <- "mr_1to_many_forest_pretty"
if (!dir.exists(forest_comprehensive_dir)) dir.create(forest_comprehensive_dir, recursive = TRUE)

# Generate comprehensive forest plot
comprehensive_plot <- forest_plot_1_to_many(
  results_formatted,
  b = "b", se = "se",
  exponentiate = TRUE,
  ao_slc = FALSE,
  lo = or_lower, up = or_upper,
  TraitM = "outcome_label",
  by = "group",
  trans = "log2",
  xlab = "Odds ratio per unit increase in genetic liability to endometriosis (95% CI)",
  weight = "weight_plot",
  subheading_size = 10,
  col1_title = "Outcome",
  col1_width = 4.8,
  col_text_size = 3.2,
  addcols = c("nsnp", "Method", "Effect_CI", "P_value"),
  addcol_widths = c(0.9, 1.0, 2.1, 1.2),
  addcol_titles = c("No. SNPs", "Method", "Effect (95% CI)", "P-value")
) +
  theme(
    plot.margin = margin(8, 14, 8, 10),
    strip.text.y = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  ggtitle("Mendelian Randomization: Endometriosis and Pregnancy Outcomes")

# Calculate figure dimensions
n_rows <- nrow(results_formatted)
fig_height <- max(10, min(0.35 * n_rows, 18))
fig_width <- 11

# Save comprehensive forest plot
ggsave(file.path(forest_comprehensive_dir, "mr_comprehensive_forest.png"), 
       comprehensive_plot, width = fig_width, height = fig_height, dpi = 300)
ggsave(file.path(forest_comprehensive_dir, "mr_comprehensive_forest.pdf"), 
       comprehensive_plot, width = fig_width, height = fig_height)
ggsave(file.path(forest_comprehensive_dir, "mr_comprehensive_forest.svg"), 
       comprehensive_plot, width = fig_width, height = fig_height)

message("Saved comprehensive forest plot")

################################################################################
# ANALYSIS SUMMARY
################################################################################

message("\n=== MR VISUALIZATION ANALYSIS COMPLETE ===")
message("Generated plots:")
message("- Scatter plots: ", length(scatter_plots[!sapply(scatter_plots, is.null)]), " files")
message("- Forest plots: ", length(forest_plots[!sapply(forest_plots, is.null)]), " files")
message("- Leave-one-out plots: ", length(loo_plots), " files")
message("- Funnel plots: ", length(funnel_plots), " files")
message("- Comprehensive forest plot: 1 file")
message("\nAll plots saved in high resolution (300 DPI) for publication.")

# Display comprehensive forest plot
print(comprehensive_plot)

<function_calls>
<invoke name="TodoWrite">
<parameter name="todos">[{"id": "1", "content": "Restructure and clean R code for MR visualization", "status": "completed"}]