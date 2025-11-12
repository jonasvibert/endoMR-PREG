#!/usr/bin/env Rscript
###############################################################################
# Filter to 30 pregnancy outcome variables for endoMR-PREG analysis
# Performs MR analysis with FDR correction and creates forest plot
# Author: Modified by Claude Code
# Date: 2025-10-12
###############################################################################

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(readr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(grid)
  library(purrr)
})

# ---------- Setup ----------
results_dir <- here("results")
plots_dir   <- here("plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

message("=== 30 PREGNANCY OUTCOMES MR ANALYSIS WITH FDR CORRECTION ===")

# ---------- Load Data ----------
# Load harmonised data
harm_file <- file.path(results_dir, "harmonised_rahmioglu_bpo.csv")
if (!file.exists(harm_file)) {
  stop("Harmonised data not found: ", harm_file, "\nRun harmonisation script first.")
}

dat <- readr::read_csv(harm_file, show_col_types = FALSE)
message("Loaded harmonised data: ", nrow(dat), " SNP-outcome pairs")

# ---------- Define 30 optimized pregnancy outcomes ----------
outcomes_optimized <- c(
  "pretb_all", "el_cs", "rup_memb", "pretb_subsamp", "em_cs", 
  "lowapgar1", "pe_subsamp", "lbw_all", "ga_all", "gh_subsamp", 
  "lga", "sb_subsamp", "hdp_subsamp", "depr_subsamp", "hbw_all", 
  "nicu", "lowapgar5", "posttb_all", "cs", "anaemia_preg_all", 
  "sga", "induction", "finngen_R12_O15_PLAC_PRAEVIA", 
  "finngen_R12_O15_PLAC_PREMAT_SEPAR", "finngen_R12_O15_PLAC_DISORD", 
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered", 
  "Postpartum_hemorrhage_filtered", "Postpartum_hemorrhage_due_to_atony_filtered", 
  "Antepartum_bleeding_filtered"
)

# Create outcome labels mapping
outcome_labels <- c(
  "pretb_all" = "Preterm birth (all)",
  "el_cs" = "Elective caesarean section",
  "rup_memb" = "Premature rupture of membranes",
  "pretb_subsamp" = "Preterm birth (subsample)",
  "em_cs" = "Emergency caesarean section",
  "lowapgar1" = "Low Apgar score at 1 min",
  "pe_subsamp" = "Preeclampsia (subsample)",
  "lbw_all" = "Low birthweight (<2500 g)",
  "ga_all" = "Gestational age (all)",
  "gh_subsamp" = "Gestational hypertension (subsample)",
  "lga" = "Large for gestational age",
  "sb_subsamp" = "Stillbirth (subsample)",
  "hdp_subsamp" = "Hypertensive disorders of pregnancy (subsample)",
  "depr_subsamp" = "Postpartum depression (subsample)",
  "hbw_all" = "High birthweight (>4000 g)",
  "nicu" = "NICU admission",
  "lowapgar5" = "Low Apgar score at 5 min",
  "posttb_all" = "Post-term birth",
  "cs" = "Caesarean section (all)",
  "anaemia_preg_all" = "Anaemia in pregnancy (all)",
  "sga" = "Small for gestational age",
  "induction" = "Induction of labour",
  "finngen_R12_O15_PLAC_PRAEVIA" = "Placenta praevia",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR" = "Placental abruption",
  "finngen_R12_O15_PLAC_DISORD" = "Other placental disorders",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered" = "PPH – retained placenta",
  "Postpartum_hemorrhage_filtered" = "Postpartum haemorrhage",
  "Postpartum_hemorrhage_due_to_atony_filtered" = "PPH – atony",
  "Antepartum_bleeding_filtered" = "Antepartum bleeding"
)

# ---------- Filter data to 30 outcomes ----------
dat_filtered <- dat %>%
  filter(outcome %in% outcomes_optimized)

available_outcomes <- unique(dat_filtered$outcome)
missing_outcomes <- setdiff(outcomes_optimized, available_outcomes)

message("Available outcomes: ", length(available_outcomes), "/", length(outcomes_optimized))
if (length(missing_outcomes) > 0) {
  message("Missing outcomes: ", paste(missing_outcomes, collapse = ", "))
}

if (nrow(dat_filtered) == 0) {
  stop("No data available for specified outcomes")
}

# ---------- Perform MR Analysis ----------
message("Running MR analysis for ", length(available_outcomes), " outcomes...")

# Function to run MR for a single outcome
run_mr_single <- function(outcome_name) {
  tryCatch({
    outcome_dat <- dat_filtered %>% filter(outcome == outcome_name)
    
    if (nrow(outcome_dat) < 3) {
      return(NULL)
    }
    
    # Run IVW
    ivw_res <- mr(outcome_dat, method_list = c("mr_ivw"))
    
    if (nrow(ivw_res) == 0) {
      return(NULL)
    }
    
    # Add outcome label
    ivw_res$outcome_clean <- outcome_labels[outcome_name]
    if (is.na(ivw_res$outcome_clean)) {
      ivw_res$outcome_clean <- outcome_name
    }
    
    return(ivw_res)
  }, error = function(e) {
    message("Error with outcome ", outcome_name, ": ", e$message)
    return(NULL)
  })
}

# Run MR for all outcomes
mr_results_list <- map(available_outcomes, run_mr_single)
mr_results_list <- compact(mr_results_list)

if (length(mr_results_list) == 0) {
  stop("No MR results obtained")
}

# Combine results
all_mr_results <- bind_rows(mr_results_list)

message("MR analysis completed for ", nrow(all_mr_results), " outcomes")

# ---------- Apply FDR Correction ----------
message("Applying FDR correction...")

all_mr_results <- all_mr_results %>%
  mutate(
    # Calculate FDR-adjusted p-values
    pval_fdr = p.adjust(pval, method = "fdr"),
    # Create significance categories
    sig_nominal = pval < 0.05,
    sig_fdr = pval_fdr < 0.05,
    # Calculate OR and CI
    OR = exp(b),
    OR_lower = exp(b - 1.96 * se),
    OR_upper = exp(b + 1.96 * se),
    # Format for display
    OR_fmt = sprintf("%.2f", OR),
    CI_fmt = sprintf("(%.2f, %.2f)", OR_lower, OR_upper),
    p_fmt = case_when(
      pval < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", pval)
    ),
    p_fdr_fmt = case_when(
      pval_fdr < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", pval_fdr)
    ),
    # Color coding for plots
    color_group = case_when(
      sig_fdr & OR > 1 ~ "FDR significant (risk)",
      sig_fdr & OR < 1 ~ "FDR significant (protective)",
      sig_nominal & OR > 1 ~ "Nominally significant (risk)",
      sig_nominal & OR < 1 ~ "Nominally significant (protective)",
      TRUE ~ "Non-significant"
    )
  ) %>%
  arrange(pval_fdr, pval)

# ---------- Save Results ----------
write_csv(all_mr_results, file.path(results_dir, "mr_results_30_outcomes_fdr.csv"))

message("Results summary:")
message("- Total outcomes analyzed: ", nrow(all_mr_results))
message("- Nominally significant (p < 0.05): ", sum(all_mr_results$sig_nominal))
message("- FDR significant (q < 0.05): ", sum(all_mr_results$sig_fdr))

# ---------- Create Forest Plot ----------
message("Creating forest plot...")

# Prepare data for forest plot - order by significance (FDR first, then nominal)
forest_data <- all_mr_results %>%
  # Order by FDR significance first, then nominal significance, then effect size
  arrange(pval_fdr, pval, desc(abs(b))) %>%
  mutate(
    # Create factor with levels REVERSED so most significant appears at TOP in ggplot
    outcome_clean = factor(outcome_clean, levels = rev(unique(outcome_clean))),
    # Create abbreviated outcome names for better display
    outcome_short = str_trunc(outcome_clean, 35),
    # Add row numbers for positioning (1 = most significant, should be at top)
    row_order = 1:n()
  )

# Color palette
color_palette <- c(
  "FDR significant (risk)" = "#d73027",
  "FDR significant (protective)" = "#1a9850", 
  "Nominally significant (risk)" = "#fc8d59",
  "Nominally significant (protective)" = "#91bfdb",
  "Non-significant" = "#999999"
)

# Calculate plot dimensions
n_outcomes <- nrow(forest_data)
plot_height <- max(8, min(20, n_outcomes * 0.4 + 4))

# Set up axis limits to include space for text columns
xmin <- 0.1
xmax <- 8.0
text_start <- 15.0  # Start position for text columns (on linear scale)

# Column positions for text
x_or <- text_start
x_ci <- text_start * 1.8
x_p <- text_start * 2.8
x_fdr <- text_start * 3.5

# Header position
y_header <- n_outcomes + 0.8

# Create data frame for text positioning
# Match the y positions exactly to how ggplot will display the factors
text_positions <- forest_data %>%
  mutate(
    # Get the numeric position of each factor level (1 = bottom, n = top)
    y_numeric = as.numeric(outcome_clean),
    or_text = OR_fmt,
    ci_text = CI_fmt,
    p_text = p_fmt,
    fdr_text = p_fdr_fmt
  )

# Create base forest plot - use factor levels as-is (no fct_rev needed)
forest_plot <- ggplot(forest_data, aes(y = outcome_clean)) +
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper, color = color_group),
                 height = 0.3, linewidth = 0.8, alpha = 0.8) +
  geom_point(aes(x = OR, color = color_group), size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6, color = "gray60") +
  
  # Add column headers
  annotate("text", x = x_or, y = y_header, label = "OR", 
           hjust = 0.5, fontface = "bold", size = 4) +
  annotate("text", x = x_ci, y = y_header, label = "95% CI", 
           hjust = 0.5, fontface = "bold", size = 4) +
  annotate("text", x = x_p, y = y_header, label = "P-value", 
           hjust = 0.5, fontface = "bold", size = 4) +
  annotate("text", x = x_fdr, y = y_header, label = "FDR P", 
           hjust = 0.5, fontface = "bold", size = 4) +
  
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4),
    labels = c("0.25", "0.5", "1", "2", "4"),
    limits = c(xmin, x_fdr * 1.2)  # Extended to include text columns
  ) +
  scale_color_manual(values = color_palette, name = "Significance") +
  
  labs(
    title = "Association of Genetically Predicted Endometriosis with 30 Pregnancy Outcomes",
    subtitle = paste("Two-sample Mendelian randomization (IVW method) with FDR correction (n =", 
                     sum(all_mr_results$sig_fdr), "FDR significant outcomes)"),
    x = "Odds Ratio (log scale)",
    y = NULL
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(20, 50, 20, 20)
  )

# Add text annotations for each row
for (i in 1:nrow(text_positions)) {
  row <- text_positions[i, ]
  forest_plot <- forest_plot +
    annotate("text", x = x_or, y = row$y_numeric, label = row$or_text,
             hjust = 0.5, size = 3.2) +
    annotate("text", x = x_ci, y = row$y_numeric, label = row$ci_text,
             hjust = 0.5, size = 3.2) +
    annotate("text", x = x_p, y = row$y_numeric, label = row$p_text,
             hjust = 0.5, size = 3.2) +
    annotate("text", x = x_fdr, y = row$y_numeric, label = row$fdr_text,
             hjust = 0.5, size = 3.2)
}

# Save forest plot
ggsave(file.path(plots_dir, "forest_30_outcomes_fdr.png"),
       forest_plot, width = 16, height = plot_height, dpi = 300, bg = "white")

message("Forest plot saved: ", file.path(plots_dir, "forest_30_outcomes_fdr.png"))

# ---------- Summary ----------
message("\n=== ANALYSIS COMPLETED ===")
message("Results saved to: ", file.path(results_dir, "mr_results_30_outcomes_fdr.csv"))
message("Forest plot saved to: ", file.path(plots_dir, "forest_30_outcomes_fdr.png"))
message("\nFDR-significant outcomes:")

fdr_sig <- all_mr_results %>% filter(sig_fdr)
if (nrow(fdr_sig) > 0) {
  for (i in 1:nrow(fdr_sig)) {
    message(sprintf("  %s: OR = %.2f (%.2f-%.2f), P = %s, FDR P = %s",
                   fdr_sig$outcome_clean[i], fdr_sig$OR[i], 
                   fdr_sig$OR_lower[i], fdr_sig$OR_upper[i],
                   fdr_sig$p_fmt[i], fdr_sig$p_fdr_fmt[i]))
  }
} else {
  message("  No outcomes significant after FDR correction")
}

message("Analysis completed successfully!")