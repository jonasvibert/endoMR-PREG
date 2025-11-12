#!/usr/bin/env Rscript
###############################################################################
# Create Supplementary Table S2: Sensitivity Analysis Results
# 
# This script generates a comprehensive sensitivity analysis table for the 
# main findings of the endoMR-PREG study, including:
# - IVW (primary analysis)
# - MR-Egger regression  
# - Weighted median
# - Heterogeneity assessment (Cochran's Q)
# - Pleiotropy testing (MR-Egger intercept)
#
# Focus on statistically significant and nominally significant associations
###############################################################################

# Load required packages
library(dplyr)
library(readr)
library(data.table)
library(knitr)
library(here)

# Install and load table formatting packages
if (!requireNamespace("flextable", quietly = TRUE)) {
  install.packages("flextable")
}
if (!requireNamespace("officer", quietly = TRUE)) {
  install.packages("officer")
}
library(flextable)
library(officer)

# Set working directory and paths
results_dir <- here::here("results")

###############################################################################
# 1. LOAD ALL SENSITIVITY ANALYSIS RESULTS
###############################################################################

message("Loading sensitivity analysis results...")

# Load primary IVW results
ivw_results <- fread(file.path(results_dir, "ivw_results.csv"))

# Load MR-Egger results
egger_results <- fread(file.path(results_dir, "egger_results.csv"))

# Load weighted median results
weighted_median_results <- fread(file.path(results_dir, "weighted_median_results.csv"))

# Load heterogeneity results
heterogeneity_results <- fread(file.path(results_dir, "heterogeneity_results.csv"))

# Load pleiotropy results
pleiotropy_results <- fread(file.path(results_dir, "pleiotropy_results.csv"))

message("All sensitivity analysis files loaded successfully.")

###############################################################################
# 2. SELECT OUTCOMES OF INTEREST
###############################################################################

# Focus on statistically significant and nominally significant associations
# Based on Bonferroni threshold (p < 0.001) and nominal significance (p < 0.05)

outcomes_of_interest <- c(
  "finngen_R12_N14_FEMALEINFERT",           # Significant (p = 6.07e-22)
  "finngen_R12_O15_PLAC_PRAEVIA",           # Significant (p = 1.46e-06) 
  "rup_memb",                               # Nominal (p = 0.025)
  "el_cs",                                  # Nominal (p = 0.023)
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",     # Nominal (p = 0.031)
  "lowapgar1"                               # Nominal (p = 0.073)
)

# Outcome labels for the table
outcome_labels <- c(
  "finngen_R12_N14_FEMALEINFERT" = "Female infertility",
  "finngen_R12_O15_PLAC_PRAEVIA" = "Placenta praevia",
  "rup_memb" = "Premature rupture of membranes",
  "el_cs" = "Elective caesarean section",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR" = "Placental abruption",
  "lowapgar1" = "Low Apgar score at 1 minute"
)

###############################################################################
# 3. FILTER AND FORMAT RESULTS
###############################################################################

# Helper function to calculate OR and 95% CI from beta and SE
calculate_or_ci <- function(beta, se) {
  or <- exp(beta)
  ci_lower <- exp(beta - 1.96 * se)
  ci_upper <- exp(beta + 1.96 * se)
  
  # Format OR and CI
  or_formatted <- sprintf("%.2f", or)
  ci_formatted <- sprintf("%.2f–%.2f", ci_lower, ci_upper)
  or_ci <- paste0(or_formatted, " (", ci_formatted, ")")
  
  return(or_ci)
}

# Helper function to format p-values
format_pvalue <- function(pval) {
  if (is.na(pval)) return("NA")
  if (pval < 0.001) {
    return(formatC(pval, format = "e", digits = 2))
  } else {
    return(sprintf("%.3f", pval))
  }
}

# Filter IVW results for outcomes of interest
ivw_filtered <- ivw_results %>%
  filter(outcome %in% outcomes_of_interest) %>%
  mutate(
    outcome_label = outcome_labels[outcome],
    or_ci_ivw = calculate_or_ci(b, se),
    pval_ivw = sapply(pval, format_pvalue)
  ) %>%
  select(outcome, outcome_label, nsnp, or_ci_ivw, pval_ivw)

# Filter MR-Egger results
egger_filtered <- egger_results %>%
  filter(outcome %in% outcomes_of_interest) %>%
  mutate(
    or_ci_egger = calculate_or_ci(b, se),
    pval_egger = sapply(pval, format_pvalue)
  ) %>%
  select(outcome, or_ci_egger, pval_egger)

# Filter weighted median results
wm_filtered <- weighted_median_results %>%
  filter(outcome %in% outcomes_of_interest) %>%
  mutate(
    # Convert string numbers to numeric
    b_numeric = as.numeric(b),
    se_numeric = as.numeric(se),
    pval_numeric = as.numeric(pval),
    or_ci_wm = calculate_or_ci(b_numeric, se_numeric),
    pval_wm = sapply(pval_numeric, format_pvalue)
  ) %>%
  select(outcome, or_ci_wm, pval_wm)

# Filter heterogeneity results (IVW method only)
heterogeneity_filtered <- heterogeneity_results %>%
  filter(outcome %in% outcomes_of_interest, method == "Inverse variance weighted") %>%
  mutate(
    q_pval = sapply(Q_pval, format_pvalue)
  ) %>%
  select(outcome, Q, Q_df, q_pval)

# Filter pleiotropy results (MR-Egger intercept)
pleiotropy_filtered <- pleiotropy_results %>%
  filter(outcome %in% outcomes_of_interest) %>%
  mutate(
    intercept_pval = sapply(pval, format_pvalue)
  ) %>%
  select(outcome, egger_intercept, intercept_pval)

###############################################################################
# 4. MERGE ALL RESULTS INTO ONE TABLE
###############################################################################

message("Merging sensitivity analysis results...")

# Merge all results
sensitivity_table <- ivw_filtered %>%
  left_join(egger_filtered, by = "outcome") %>%
  left_join(wm_filtered, by = "outcome") %>%
  left_join(heterogeneity_filtered, by = "outcome") %>%
  left_join(pleiotropy_filtered, by = "outcome")

# Reorder by significance (significant first, then nominal)
significance_order <- c(
  "finngen_R12_N14_FEMALEINFERT",
  "finngen_R12_O15_PLAC_PRAEVIA", 
  "el_cs",
  "rup_memb",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",
  "lowapgar1"
)

sensitivity_table <- sensitivity_table %>%
  arrange(factor(outcome, levels = significance_order))

# Format final table with proper column names
final_table <- sensitivity_table %>%
  select(
    `Outcome` = outcome_label,
    `No. SNPs` = nsnp,
    `IVW OR (95% CI)` = or_ci_ivw,
    `IVW P-value` = pval_ivw,
    `MR-Egger OR (95% CI)` = or_ci_egger,
    `MR-Egger P-value` = pval_egger,
    `Weighted Median OR (95% CI)` = or_ci_wm,
    `Weighted Median P-value` = pval_wm,
    `Cochran's Q` = Q,
    `Q P-value` = q_pval,
    `MR-Egger Intercept` = egger_intercept,
    `Intercept P-value` = intercept_pval
  )

###############################################################################
# 5. CREATE FORMATTED TABLE FOR MANUSCRIPT
###############################################################################

message("Creating formatted table...")

# Create flextable
ft <- flextable(final_table)

# Apply formatting
ft <- ft %>%
  # Set table width and column widths
  width(width = c(2.0, 0.6, 1.2, 0.8, 1.2, 0.8, 1.2, 0.8, 0.8, 0.8, 1.0, 0.8)) %>%
  
  # Header formatting
  bg(bg = "#F0F0F0", part = "header") %>%
  bold(part = "header") %>%
  align(align = "center", part = "header") %>%
  
  # Body formatting
  align(j = 1, align = "left", part = "body") %>%   # Outcome names left-aligned
  align(j = 2:12, align = "center", part = "body") %>%  # Numbers centered
  
  # Font and size
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  fontsize(size = 11, part = "header") %>%
  
  # Borders
  border_outer(border = fp_border(color = "black", width = 1), part = "all") %>%
  border_inner_h(border = fp_border(color = "gray", width = 0.5), part = "all") %>%
  border_inner_v(border = fp_border(color = "gray", width = 0.5), part = "all") %>%
  
  # Add padding
  padding(padding = 3, part = "all")

# Add footnote
footnote_text <- paste0(
  "Abbreviations: IVW, inverse variance weighted; OR, odds ratio; CI, confidence interval; SNP, single nucleotide polymorphism.\n",
  "The MR-Egger intercept test assesses directional pleiotropy (P < 0.05 indicates potential pleiotropy). ",
  "Cochran's Q test evaluates heterogeneity across genetic instruments (P < 0.05 indicates significant heterogeneity). ",
  "Results are shown for outcomes with P < 0.05 in the primary IVW analysis. ",
  "Bonferroni-corrected significance threshold: P < 0.001."
)

###############################################################################
# 6. EXPORT RESULTS
###############################################################################

message("Exporting Supplementary Table S2...")

# Export as CSV for easy viewing
write_csv(final_table, file.path(results_dir, "Supplementary_Table_S2_Sensitivity_Analysis.csv"))

# Export as Word document
doc <- read_docx()
doc <- doc %>%
  body_add_par("Supplementary Table S2. Sensitivity Analysis Results for Mendelian Randomization of Endometriosis and Pregnancy Outcomes", 
               style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft) %>%
  body_add_par("") %>%
  body_add_par(footnote_text, style = "Normal")

print(doc, target = file.path(results_dir, "Supplementary_Table_S2_Sensitivity_Analysis.docx"))

# Export as HTML for easy viewing (skip if pandoc not available)
tryCatch({
  save_as_html(ft, path = file.path(results_dir, "Supplementary_Table_S2_Sensitivity_Analysis.html"))
  message("HTML file created successfully")
}, error = function(e) {
  message("HTML export skipped (pandoc not available)")
})

###############################################################################
# 7. SUMMARY STATISTICS
###############################################################################

message("Generating summary of sensitivity analysis results...")

# Summary of key findings
summary_stats <- list(
  outcomes_analyzed = nrow(final_table),
  significant_outcomes = sum(as.numeric(gsub("[^0-9.]", "", sensitivity_table$pval_ivw)) < 0.001, na.rm = TRUE),
  nominal_outcomes = sum(as.numeric(gsub("[^0-9.]", "", sensitivity_table$pval_ivw)) < 0.05 & 
                        as.numeric(gsub("[^0-9.]", "", sensitivity_table$pval_ivw)) >= 0.001, na.rm = TRUE),
  heterogeneity_detected = sum(as.numeric(sensitivity_table$q_pval) < 0.05, na.rm = TRUE),
  pleiotropy_detected = sum(as.numeric(sensitivity_table$intercept_pval) < 0.05, na.rm = TRUE)
)

# Print summary
cat("\n=== SUPPLEMENTARY TABLE S2 SUMMARY ===\n")
cat("Outcomes analyzed:", summary_stats$outcomes_analyzed, "\n")
cat("Bonferroni-significant outcomes (p < 0.001):", summary_stats$significant_outcomes, "\n") 
cat("Nominally significant outcomes (p < 0.05):", summary_stats$nominal_outcomes, "\n")
cat("Outcomes with significant heterogeneity:", summary_stats$heterogeneity_detected, "\n")
cat("Outcomes with evidence of pleiotropy:", summary_stats$pleiotropy_detected, "\n")

# Save summary
saveRDS(summary_stats, file.path(results_dir, "supplementary_table_s2_summary.rds"))

message("Supplementary Table S2 creation completed successfully!")
message("Files created:")
message("- ", file.path(results_dir, "Supplementary_Table_S2_Sensitivity_Analysis.csv"))
message("- ", file.path(results_dir, "Supplementary_Table_S2_Sensitivity_Analysis.docx"))
message("- ", file.path(results_dir, "Supplementary_Table_S2_Sensitivity_Analysis.html"))

###############################################################################
# END OF SCRIPT
###############################################################################