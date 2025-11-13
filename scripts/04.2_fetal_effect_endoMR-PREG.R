+ ################################################################################
+ # FETAL GENETIC EFFECT ANALYSIS
  + # Endometriosis and Pregnancy Outcomes Study
  + # 
  + # Purpose: Analyze fetal genetic effects on pregnancy outcomes using MR
  + # Author: [Your Name]
  + # Date: August 2025
  + # 
  + # Background: 
  + # In pregnancy outcomes, genetic effects can arise from:
  + # 1. Maternal genetics (analyzed in main script)
  + # 2. Fetal genetics (analyzed here)
  + # 3. Maternal-fetal genetic interactions
  + #
  + # This script focuses on fetal genetic contributions to pregnancy outcomes
  + # using available fetal/birth outcome GWAS data
  + ################################################################################
+ 
  + # Load required libraries
library(TwoSampleMR)
+ library(dplyr)
+ library(purrr)
+ library(ggplot2)
+ library(stringr)
+ library(data.table)
+ 
  + ################################################################################
+ # SECTION 1: FETAL GENETIC DATA PREPARATION
  + ################################################################################
+ 
  + message("=== FETAL GENETIC EFFECT ANALYSIS ===")
+ 
  + # Define fetal-relevant outcomes (birth outcomes influenced by fetal genetics)
  + # Using standardized outcome codes from the main analysis (filtered to 29 outcomes)
  + fetal_outcomes <- c(
    +   # Birth size outcomes (strong fetal genetic component)
      +   "lbw_all",    # Low birthweight (<2500g)
    +   "hbw_all",    # High birthweight (>4000g)
    +   "sga",        # Small for gestational age
    +   
      +   # Gestational timing (potential fetal contribution)
      +   "pretb_all",  # Preterm birth (all)
    +   "vpretb_all", # Very preterm birth
    +   "ga_all",     # Gestational age
    +   "rup_memb",   # Premature rupture of membranes
    +   
      +   # Neonatal outcomes (fetal development-related)
      +   "lowapgar1",  # Low Apgar at 1 minute
    +   "lowapgar5",  # Low Apgar at 5 minutes
    +   "nicu"        # NICU admission
    + )
  + 
    + # Fetal outcome labels (updated to match filtered 29 outcomes)
    + fetal_labels <- c(
      +   "lbw_all" = "Low birthweight (<2500g)",
      +   "hbw_all" = "High birthweight (>4000g)",
      +   "sga" = "Small for gestational age",
      +   "pretb_all" = "Preterm birth (all)",
      +   "vpretb_all" = "Very preterm birth",
      +   "ga_all" = "Gestational age",
      +   "rup_memb" = "Premature rupture of membranes",
      +   "lowapgar1" = "Low Apgar at 1 minute",
      +   "lowapgar5" = "Low Apgar at 5 minutes",
      +   "nicu" = "NICU admission"
      + )
    + 
      + # Clinical groupings for fetal effects (updated to match filtered 29 outcomes)
      + fetal_groups <- c(
        +   "lbw_all" = "Birth size", "hbw_all" = "Birth size", "sga" = "Birth size",
        +   "pretb_all" = "Gestational timing", "vpretb_all" = "Gestational timing", 
        +   "ga_all" = "Gestational timing", "rup_memb" = "Gestational timing",
        +   "lowapgar1" = "Neonatal health", "lowapgar5" = "Neonatal health", "nicu" = "Neonatal health"
        + )
      + 
        + ################################################################################
      + # SECTION 2: FETAL GENETIC INSTRUMENT PREPARATION
        + ################################################################################
      + 
        + # Function to prepare fetal genetic instruments
        + prepare_fetal_instruments <- function() {
          +   
            +   message("Preparing fetal genetic instruments...")
          +   
            +   # Note: In practice, you would need fetal GWAS data
            +   # This example assumes you have prepared fetal instruments
            +   # Sources might include:
            +   # - Birth weight GWAS (maternal vs fetal effects separated)
            +   # - Gestational duration GWAS
            +   # - Fetal growth GWAS
            +   
            +   # Example structure for fetal instruments
            +   # Replace with your actual fetal genetic data
            +   
            +   if (!exists("fetal_instruments")) {
              +     message("Loading fetal genetic instruments...")
              +     
                +     # Example: Load fetal birthweight instruments
                +     # fetal_instruments <- fread("data/fetal_birthweight_instruments.txt")
                +     
                +     # For demonstration, create placeholder structure
                +     fetal_instruments <- data.frame(
                  +       SNP = character(0),
                  +       beta = numeric(0),
                  +       se = numeric(0),
                  +       effect_allele = character(0),
                  +       other_allele = character(0),
                  +       eaf = numeric(0),
                  +       pval = numeric(0),
                  +       exposure = character(0),
                  +       stringsAsFactors = FALSE
                  +     )
                +     
                  +     message("Warning: No fetal genetic instruments loaded.")
                +     message("Please prepare fetal GWAS instruments and load into 'fetal_instruments'")
                +     return(NULL)
                +   }
          +   
            +   # Format fetal instruments for TwoSampleMR
            +   fetal_exp <- format_data(
              +     fetal_instruments,
              +     type = "exposure",
              +     snp_col = "SNP",
              +     beta_col = "beta",
              +     se_col = "se",
              +     effect_allele_col = "effect_allele",
              +     other_allele_col = "other_allele",
              +     eaf_col = "eaf",
              +     pval_col = "pval"
              +   )
            +   
              +   return(fetal_exp)
            + }
        + 
          + ################################################################################
        + # SECTION 3: FETAL GENETIC EFFECT ANALYSIS
          + ################################################################################
        + 
          + analyze_fetal_effects <- function(fetal_exposure_data, outcome_data) {
            +   
              +   if (is.null(fetal_exposure_data)) {
                +     message("Skipping fetal analysis - no fetal instruments available")
                +     return(NULL)
                +   }
            +   
              +   message("Running fetal genetic effect analysis...")
            +   
              +   # Harmonize fetal exposure with outcomes
              +   fetal_harmonised <- harmonise_data(
                +     exposure_dat = fetal_exposure_data,
                +     outcome_dat = outcome_data
                +   )
              +   
                +   # Filter for fetal-relevant outcomes
                +   fetal_harmonised <- fetal_harmonised %>%
                  +     filter(id.outcome %in% fetal_outcomes)
                +   
                  +   if (nrow(fetal_harmonised) == 0) {
                    +     message("No harmonised data for fetal analysis")
                    +     return(NULL)
                    +   }
                +   
                  +   # Run MR analysis
                  +   fetal_mr_results <- mr(fetal_harmonised)
                  +   
                    +   # Add fetal-specific formatting
                    +   fetal_results_formatted <- fetal_mr_results %>%
                      +     filter(method == "Inverse variance weighted") %>%
                      +     mutate(
                        +       outcome_label = dplyr::recode(id.outcome, !!!fetal_labels),
                        +       group = dplyr::recode(id.outcome, !!!fetal_groups),
                        +       OR = exp(b),
                        +       LCL = exp(b - 1.96 * se),
                        +       UCL = exp(b + 1.96 * se),
                        +       Effect_CI = sprintf("%.2f (%.2f-%.2f)", OR, LCL, UCL),
                        +       P_value = formatC(pval, format = "e", digits = 2),
                        +       Analysis_type = "Fetal genetic effect"
                        +     )
                    +   
                      +   return(list(
                        +     harmonised = fetal_harmonised,
                        +     results = fetal_mr_results,
                        +     formatted = fetal_results_formatted
                        +   ))
                    + }
          + 
            + ################################################################################
          + # SECTION 4: MATERNAL VS FETAL EFFECT COMPARISON
            + ################################################################################
          + 
            + compare_maternal_fetal_effects <- function(maternal_results, fetal_results) {
              +   
                +   if (is.null(fetal_results)) {
                  +     message("Cannot compare - no fetal results available")
                  +     return(NULL)
                  +   }
              +   
                +   message("Comparing maternal vs fetal genetic effects...")
              +   
                +   # Prepare maternal results for comparison
                +   maternal_formatted <- maternal_results %>%
                  +     filter(
                    +       method == "Inverse variance weighted",
                    +       id.outcome %in% fetal_outcomes
                    +     ) %>%
                  +     mutate(
                    +       outcome_label = dplyr::recode(id.outcome, !!!fetal_labels),
                    +       group = dplyr::recode(id.outcome, !!!fetal_groups),
                    +       OR = exp(b),
                    +       LCL = exp(b - 1.96 * se),
                    +       UCL = exp(b + 1.96 * se),
                    +       Effect_CI = sprintf("%.2f (%.2f-%.2f)", OR, LCL, UCL),
                    +       P_value = formatC(pval, format = "e", digits = 2),
                    +       Analysis_type = "Maternal genetic effect"
                    +     )
                +   
                  +   # Combine maternal and fetal results
                  +   comparison_data <- bind_rows(
                    +     maternal_formatted %>% select(id.outcome, outcome_label, group, b, se, OR, LCL, UCL, 
                                                        +                                   Effect_CI, P_value, Analysis_type, nsnp),
                    +     fetal_results$formatted %>% select(id.outcome, outcome_label, group, b, se, OR, LCL, UCL,
                                                             +                                        Effect_CI, P_value, Analysis_type, nsnp)
                    +   ) %>%
                    +     arrange(group, id.outcome, Analysis_type)
                  +   
                    +   return(comparison_data)
                  + }
            + 
              + ################################################################################
            + # SECTION 5: FETAL EFFECT VISUALIZATION
              + ################################################################################
            + 
              + create_fetal_comparison_plot <- function(comparison_data) {
                +   
                  +   if (is.null(comparison_data)) {
                    +     message("Cannot create comparison plot - no comparison data")
                    +     return(NULL)
                    +   }
                +   
                  +   message("Creating maternal vs fetal effect comparison plot...")
                +   
                  +   # Create output directory
                  +   fetal_dir <- "fetal_genetic_analysis"
                  +   if (!dir.exists(fetal_dir)) dir.create(fetal_dir, recursive = TRUE)
                  +   
                    +   # Create comparison forest plot
                    +   comparison_plot <- ggplot(comparison_data, aes(x = OR, y = outcome_label)) +
                      +     geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
                      +     geom_errorbarh(aes(xmin = LCL, xmax = UCL, color = Analysis_type), 
                                           +                    height = 0.2, position = position_dodge(width = 0.5)) +
                      +     geom_point(aes(color = Analysis_type), size = 3, 
                                       +                position = position_dodge(width = 0.5)) +
                      +     facet_wrap(~group, scales = "free_y", ncol = 1) +
                      +     scale_color_manual(
                        +       values = c("Maternal genetic effect" = "#2E86AB", 
                                           +                 "Fetal genetic effect" = "#A23B72"),
                        +       name = "Effect Origin"
                        +     ) +
                      +     scale_x_log10(
                        +       breaks = c(0.5, 0.7, 1.0, 1.4, 2.0),
                        +       labels = c("0.5", "0.7", "1.0", "1.4", "2.0")
                        +     ) +
                      +     labs(
                        +       title = "Maternal vs Fetal Genetic Effects on Pregnancy Outcomes",
                        +       subtitle = "Endometriosis Genetic Liability",
                        +       x = "Odds Ratio (95% CI)",
                        +       y = "Pregnancy Outcome"
                        +     ) +
                      +     theme_minimal(base_size = 11) +
                      +     theme(
                        +       plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                        +       plot.subtitle = element_text(size = 12, hjust = 0.5),
                        +       strip.text = element_text(size = 11, face = "bold"),
                        +       legend.position = "bottom",
                        +       panel.grid.minor = element_blank()
                        +     )
                    +   
                      +   # Save comparison plot
                      +   ggsave(file.path(fetal_dir, "maternal_vs_fetal_effects.png"), 
                                 +          comparison_plot, width = 10, height = 8, dpi = 300)
                    +   ggsave(file.path(fetal_dir, "maternal_vs_fetal_effects.pdf"), 
                               +          comparison_plot, width = 10, height = 8)
                    +   
                      +   return(comparison_plot)
                    + }
              + 
                + ################################################################################
              + # SECTION 6: FETAL SENSITIVITY ANALYSIS
                + ################################################################################
              + 
                + fetal_sensitivity_analysis <- function(fetal_harmonised_data) {
                  +   
                    +   if (is.null(fetal_harmonised_data)) {
                      +     message("No fetal data for sensitivity analysis")
                      +     return(NULL)
                      +   }
                  +   
                    +   message("Running fetal genetic effect sensitivity analyses...")
                  +   
                    +   # Create output directory
                    +   fetal_sens_dir <- "fetal_genetic_analysis/sensitivity"
                    +   if (!dir.exists(fetal_sens_dir)) dir.create(fetal_sens_dir, recursive = TRUE)
                    +   
                      +   sensitivity_results <- list()
                      +   
                        +   # For each fetal outcome with sufficient SNPs
                        +   for (outcome_id in unique(fetal_harmonised_data$id.outcome)) {
                          +     
                            +     dat_outcome <- fetal_harmonised_data %>% filter(id.outcome == outcome_id)
                            +     
                              +     if (nrow(dat_outcome) < 3) {
                                +       message("Insufficient SNPs for sensitivity analysis: ", outcome_id)
                                +       next
                                +     }
                            +     
                              +     # Heterogeneity test
                              +     het_test <- mr_heterogeneity(dat_outcome)
                              +     
                                +     # Pleiotropy test (if ≥3 SNPs)
                                +     if (nrow(dat_outcome) >= 3) {
                                  +       pleio_test <- mr_pleiotropy_test(dat_outcome)
                                  +     } else {
                                    +       pleio_test <- NULL
                                    +     }
                              +     
                                +     # Leave-one-out analysis
                                +     loo_results <- mr_leaveoneout(dat_outcome)
                                +     
                                  +     # Store results
                                  +     sensitivity_results[[outcome_id]] <- list(
                                    +       heterogeneity = het_test,
                                    +       pleiotropy = pleio_test,
                                    +       leaveoneout = loo_results
                                    +     )
                                  +     
                                    +     # Create leave-one-out plot
                                    +     if (nrow(loo_results) > 0) {
                                      +       loo_plot <- mr_leaveoneout_plot(loo_results)[[1]] +
                                        +         ggtitle(paste0("Fetal Effect Leave-One-Out: ", fetal_labels[[outcome_id]])) +
                                        +         theme_minimal()
                                      +       
                                        +       ggsave(file.path(fetal_sens_dir, paste0(outcome_id, "_fetal_loo.png")), 
                                                       +              loo_plot, width = 8, height = 6, dpi = 300)
                                      +     }
                                  +   }
                      +   
                        +   return(sensitivity_results)
                      + }
                + 
                  + ################################################################################
                + # SECTION 7: FETAL ANALYSIS SUMMARY TABLE
                  + ################################################################################
                + 
                  + create_fetal_summary_table <- function(comparison_data, sensitivity_results) {
                    +   
                      +   if (is.null(comparison_data)) {
                        +     message("Cannot create summary table - no comparison data")
                        +     return(NULL)
                        +   }
                    +   
                      +   message("Creating fetal genetic effect summary table...")
                    +   
                      +   # Create comprehensive summary
                      +   summary_table <- comparison_data %>%
                        +     select(outcome_label, group, Analysis_type, OR, Effect_CI, P_value, nsnp) %>%
                        +     pivot_wider(
                          +       names_from = Analysis_type,
                          +       values_from = c(OR, Effect_CI, P_value, nsnp),
                          +       names_sep = "_"
                          +     ) %>%
                        +     arrange(group, outcome_label)
                      +   
                        +   # Add heterogeneity information if available
                        +   if (!is.null(sensitivity_results)) {
                          +     het_summary <- map_dfr(names(sensitivity_results), function(outcome_id) {
                            +       het_data <- sensitivity_results[[outcome_id]]$heterogeneity
                            +       if (!is.null(het_data)) {
                              +         data.frame(
                                +           outcome_id = outcome_id,
                                +           Q_pval = het_data$Q_pval[het_data$method == "Inverse variance weighted"],
                                +           stringsAsFactors = FALSE
                                +         )
                              +       }
                            +     })
                          +     
                            +     if (nrow(het_summary) > 0) {
                              +       # Add outcome labels for joining
                                +       het_summary <- het_summary %>%
                                  +         mutate(outcome_label = dplyr::recode(outcome_id, !!!fetal_labels))
                                +       
                                  +       # Join with summary table
                                  +       summary_table <- summary_table %>%
                                    +         left_join(het_summary %>% select(outcome_label, Q_pval), by = "outcome_label")
                                  +     }
                          +   }
                      +   
                        +   # Save summary table
                        +   fetal_dir <- "fetal_genetic_analysis"
                        +   write.csv(summary_table, file.path(fetal_dir, "fetal_genetic_effects_summary.csv"), 
                                      +             row.names = FALSE)
                        +   
                          +   return(summary_table)
                        + }
                  + 
                    + ################################################################################
                  + # SECTION 8: MAIN EXECUTION FUNCTION
                    + ################################################################################
                  + 
                    + run_fetal_genetic_analysis <- function(maternal_mr_results = NULL) {
                      +   
                        +   message("\n=== STARTING FETAL GENETIC EFFECT ANALYSIS ===")
                      +   
                        +   # Step 1: Prepare fetal instruments
                        +   fetal_exposure <- prepare_fetal_instruments()
                        +   
                          +   # Step 2: Load outcome data (assumes same outcomes as maternal analysis)
                          +   if (!exists("outcome_data_formatted")) {
                            +     message("Warning: outcome_data_formatted not found")
                            +     message("Please ensure outcome data is prepared from main analysis")
                            +     return(NULL)
                            +   }
                        +   
                          +   # Step 3: Run fetal analysis
                          +   fetal_analysis <- analyze_fetal_effects(fetal_exposure, outcome_data_formatted)
                          +   
                            +   # Step 4: Compare with maternal effects
                            +   comparison_data <- compare_maternal_fetal_effects(maternal_mr_results, fetal_analysis)
                            +   
                              +   # Step 5: Create visualization
                              +   comparison_plot <- create_fetal_comparison_plot(comparison_data)
                              +   
                                +   # Step 6: Sensitivity analysis
                                +   sensitivity_results <- fetal_sensitivity_analysis(fetal_analysis$harmonised)
                                +   
                                  +   # Step 7: Summary table
                                  +   summary_table <- create_fetal_summary_table(comparison_data, sensitivity_results)
                                  +   
                                    +   message("\n=== FETAL GENETIC ANALYSIS COMPLETE ===")
                                  +   
                                    +   return(list(
                                      +     fetal_results = fetal_analysis,
                                      +     comparison_data = comparison_data,
                                      +     comparison_plot = comparison_plot,
                                      +     sensitivity = sensitivity_results,
                                      +     summary_table = summary_table
                                      +   ))
                                  + }
                    + 
                      + ################################################################################
                    + # USAGE INSTRUCTIONS
                      + ################################################################################
                    + 
                      + # To run the fetal genetic analysis:
                      + # 
                      + # 1. Ensure you have fetal genetic instruments prepared
                      + # 2. Load the main MR results (maternal effects)
                      + # 3. Run: fetal_analysis_results <- run_fetal_genetic_analysis(res)
                      + #
                      + # Note: This script requires fetal GWAS data which may need to be:
                      + # - Downloaded from appropriate consortiums
                      + # - Processed to separate maternal vs fetal effects
                      + # - Formatted for TwoSampleMR package
                      + 
                      + message("Fetal genetic effect analysis script loaded.")
                    + message("Run run_fetal_genetic_analysis() to execute the analysis.")
                    + message("Ensure fetal genetic instruments are prepared first.")