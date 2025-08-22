#!/usr/bin/env Rscript
###############################################################################
# MR Visualization Script for endoMR-PREG Study                               
# Input  : harmonised data and MR results from previous analyses               
# Output : Comprehensive plots for MR analysis visualization                   
###############################################################################

### 1. SETUP and package loading ################################################
pkgs <- c("TwoSampleMR", "dplyr", "ggplot2", "here", "readr", "data.table", 
          "gridExtra", "cowplot", "scales", "RColorBrewer")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
}

# Paths
results_dir <- here::here("results")
plots_dir <- here::here("plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Load data files
harm_file <- here::here("results", "harmonised_rahmioglu_bpo.csv")
dat <- data.table::fread(harm_file) %>% as.data.frame()

# Load MR results
all_res <- readr::read_csv(file.path(results_dir, "all_mr_methods.csv"))
single_res <- readr::read_csv(file.path(results_dir, "singlesnp_results.csv"))
loo_res <- readr::read_csv(file.path(results_dir, "leaveoneout_snp_results.csv"))

### 2. OUTCOME LABELS ###########################################################
outcome_labels <- c(
  pretb_all      = "Preterm birth (all)",
  el_cs          = "Elective caesarean section",
  rup_memb       = "Premature rupture of membranes",
  pretb_subsamp  = "Preterm birth (subsample)",
  apgar1         = "Apgar score at 1 min",
  vpretb_all     = "Very preterm birth (all)",
  em_cs          = "Emergency caesarean section",
  lowapgar1      = "Low Apgar score at 1 min",
  pe_subsamp     = "Preeclampsia (subsample)",
  lbw_all        = "Low birthweight (<2500 g)",
  ga_all         = "Gestational age (all)",
  gh_subsamp     = "Gestational hypertension (subsample)",
  ga_subsamp     = "Gestational age (subsample)",
  lga            = "Large for gestational age",
  zbw_all        = "Z-score birthweight (all)",
  nvp_sev_subsamp= "Severe nausea/vomiting (subsample)",
  nvp_sev_all    = "Severe nausea/vomiting (all)",
  sb_subsamp     = "Stillbirth (subsample)",
  gdm_subsamp    = "Gestational diabetes (subsample)",
  apgar5         = "Apgar score at 5 min",
  hdp_subsamp    = "Hypertensive disorders of pregnancy (subsample)",
  r_misc_subsamp = "Recurrent miscarriage (subsample)",
  bf_dur_4c      = "Breastfeeding ≥4 months",
  depr_subsamp   = "Postpartum depression (subsample)",
  hbw_all        = "High birthweight (>4000 g)",
  nicu           = "NICU admission",
  lowapgar5      = "Low Apgar score at 5 min",
  bf_sus         = "Breastfeeding cessation",
  posttb_all     = "Post-term birth",
  misc_subsamp   = "Miscarriage (subsample)",
  bf_ini         = "Breastfeeding initiation",
  bf_est         = "Exclusive breastfeeding",
  s_misc_subsamp = "Single miscarriage (subsample)",
  cs             = "Caesarean section (all)",
  hyp            = "Hyperemesis gravidarum",
  anaemia_preg_all = "Anaemia in pregnancy (all)",
  sga            = "Small for gestational age",
  induction      = "Induction of labour",
  finngen_R12_N14_FEMALEINFERT = "Female infertility",
  finngen_R12_O15_PLAC_PRAEVIA = "Placenta praevia",
  finngen_R12_O15_PLAC_PREMAT_SEPAR = "Placental abruption",
  finngen_R12_O15_PLAC_DISORD = "Other placental disorders",
  finngen_R12_O15_PREG_ECTOP   = "Ectopic pregnancy",
  Early_bleeding_with_any_outcome_filtered = "Early bleeding (any outcome)",
  Postpartum_hemorrhage_due_to_retained_placenta_filtered = "PPH – retained placenta",
  Early_bleeding_ending_in_live_birth_filtered = "Early bleeding – live birth",
  Postpartum_hemorrhage_filtered = "Postpartum haemorrhage",
  Postpartum_hemorrhage_due_to_atony_filtered = "PPH – atony",
  Antepartum_bleeding_filtered = "Antepartum bleeding"
)

# Apply outcome labels to data
labels_df <- data.frame(
  outcome = names(outcome_labels),
  outcome_full = unname(outcome_labels),
  stringsAsFactors = FALSE
)

# Merge labels with data
dat <- merge(dat, labels_df, by = "outcome", all.x = TRUE)
dat$outcome_full[is.na(dat$outcome_full)] <- dat$outcome[is.na(dat$outcome_full)]

# Update MR results with readable names
all_res <- merge(all_res, labels_df, by = "outcome", all.x = TRUE)
all_res$outcome_full[is.na(all_res$outcome_full)] <- all_res$outcome[is.na(all_res$outcome_full)]

single_res <- merge(single_res, labels_df, by = "outcome", all.x = TRUE)
single_res$outcome_full[is.na(single_res$outcome_full)] <- single_res$outcome[is.na(single_res$outcome_full)]

loo_res <- merge(loo_res, labels_df, by = "outcome", all.x = TRUE)
loo_res$outcome_full[is.na(loo_res$outcome_full)] <- loo_res$outcome[is.na(loo_res$outcome_full)]

### 3. SCATTER PLOTS ############################################################
message("Creating scatter plots...")

# Function to create scatter plots for significant outcomes
create_scatter_plots <- function(dat, res, p_threshold = 0.05) {
  # Filter for significant IVW results
  sig_outcomes <- res %>%
    filter(method == "Inverse variance weighted", pval < p_threshold) %>%
    pull(outcome) %>%
    unique()
  
  scatter_plots <- list()
  
  for (outcome_name in sig_outcomes) {
    # Filter data for this outcome
    dat_outcome <- dat %>% filter(outcome == outcome_name)
    res_outcome <- res %>% filter(outcome == outcome_name)
    
    if (nrow(dat_outcome) > 0 && nrow(res_outcome) > 0) {
      # Create scatter plot
      p <- tryCatch({
        mr_scatter_plot(res_outcome, dat_outcome)[[1]] +
          ggtitle(paste("Endometriosis →", unique(dat_outcome$outcome_full)[1])) +
          theme_minimal() +
          theme(plot.title = element_text(size = 12, hjust = 0.5))
      }, error = function(e) {
        message("Error creating scatter plot for ", outcome_name, ": ", e$message)
        NULL
      })
      
      if (!is.null(p)) {
        scatter_plots[[outcome_name]] <- p
        
        # Save individual plot
        ggsave(
          filename = file.path(plots_dir, paste0("scatter_", gsub("[^A-Za-z0-9]", "_", outcome_name), ".png")),
          plot = p,
          width = 8, height = 6, dpi = 300
        )
      }
    }
  }
  
  return(scatter_plots)
}

scatter_plots <- create_scatter_plots(dat, all_res)

### 4. FOREST PLOTS #############################################################
message("Creating forest plots...")

# Function to create forest plots
create_forest_plots <- function(single_res, p_threshold = 0.05) {
  # Get significant outcomes from IVW results
  sig_outcomes <- single_res %>%
    filter(SNP == "All - Inverse variance weighted", pval < p_threshold) %>%
    pull(outcome) %>%
    unique()
  
  forest_plots <- list()
  
  for (outcome_name in sig_outcomes) {
    single_outcome <- single_res %>% filter(outcome == outcome_name)
    
    if (nrow(single_outcome) > 0) {
      p <- tryCatch({
        mr_forest_plot(single_outcome)[[1]] +
          ggtitle(paste("Forest Plot: Endometriosis →", unique(single_outcome$outcome_full)[1])) +
          theme_minimal() +
          theme(plot.title = element_text(size = 12, hjust = 0.5))
      }, error = function(e) {
        message("Error creating forest plot for ", outcome_name, ": ", e$message)
        NULL
      })
      
      if (!is.null(p)) {
        forest_plots[[outcome_name]] <- p
        
        # Save individual plot
        ggsave(
          filename = file.path(plots_dir, paste0("forest_", gsub("[^A-Za-z0-9]", "_", outcome_name), ".png")),
          plot = p,
          width = 10, height = 8, dpi = 300
        )
      }
    }
  }
  
  return(forest_plots)
}

forest_plots <- create_forest_plots(single_res)

### 5. LEAVE-ONE-OUT PLOTS ######################################################
message("Creating leave-one-out plots...")

# Function to create leave-one-out plots
create_loo_plots <- function(loo_res, p_threshold = 0.05) {
  # Get significant outcomes
  sig_outcomes <- loo_res %>%
    filter(SNP == "All", pval < p_threshold) %>%
    pull(outcome) %>%
    unique()
  
  loo_plots <- list()
  
  for (outcome_name in sig_outcomes) {
    loo_outcome <- loo_res %>% filter(outcome == outcome_name)
    
    if (nrow(loo_outcome) > 0) {
      p <- tryCatch({
        mr_leaveoneout_plot(loo_outcome)[[1]] +
          ggtitle(paste("Leave-One-Out: Endometriosis →", unique(loo_outcome$outcome_full)[1])) +
          theme_minimal() +
          theme(plot.title = element_text(size = 12, hjust = 0.5))
      }, error = function(e) {
        message("Error creating leave-one-out plot for ", outcome_name, ": ", e$message)
        NULL
      })
      
      if (!is.null(p)) {
        loo_plots[[outcome_name]] <- p
        
        # Save individual plot
        ggsave(
          filename = file.path(plots_dir, paste0("loo_", gsub("[^A-Za-z0-9]", "_", outcome_name), ".png")),
          plot = p,
          width = 10, height = 8, dpi = 300
        )
      }
    }
  }
  
  return(loo_plots)
}

loo_plots <- create_loo_plots(loo_res)

### 6. FUNNEL PLOTS #############################################################
message("Creating funnel plots...")

# Function to create funnel plots
create_funnel_plots <- function(single_res, p_threshold = 0.05) {
  # Get significant outcomes
  sig_outcomes <- single_res %>%
    filter(SNP == "All - Inverse variance weighted", pval < p_threshold) %>%
    pull(outcome) %>%
    unique()
  
  funnel_plots <- list()
  
  for (outcome_name in sig_outcomes) {
    single_outcome <- single_res %>% filter(outcome == outcome_name)
    
    if (nrow(single_outcome) > 0) {
      p <- tryCatch({
        mr_funnel_plot(single_outcome)[[1]] +
          ggtitle(paste("Funnel Plot: Endometriosis →", unique(single_outcome$outcome_full)[1])) +
          theme_minimal() +
          theme(plot.title = element_text(size = 12, hjust = 0.5))
      }, error = function(e) {
        message("Error creating funnel plot for ", outcome_name, ": ", e$message)
        NULL
      })
      
      if (!is.null(p)) {
        funnel_plots[[outcome_name]] <- p
        
        # Save individual plot
        ggsave(
          filename = file.path(plots_dir, paste0("funnel_", gsub("[^A-Za-z0-9]", "_", outcome_name), ".png")),
          plot = p,
          width = 8, height = 6, dpi = 300
        )
      }
    }
  }
  
  return(funnel_plots)
}

funnel_plots <- create_funnel_plots(single_res)

### 7. COMPREHENSIVE FOREST PLOT ###############################################
message("Creating comprehensive forest plot...")

# Prepare data for comprehensive forest plot
prep_forest_data <- function(all_res) {
  # Filter for IVW results and add confidence intervals
  forest_data <- all_res %>%
    filter(method == "Inverse variance weighted") %>%
    mutate(
      or = exp(b),
      ci_lower = exp(b - 1.96 * se),
      ci_upper = exp(b + 1.96 * se),
      outcome_clean = outcome_full,
      significant = pval < 0.05
    ) %>%
    arrange(desc(abs(b))) %>%
    slice_head(n = 30)  # Top 30 results for visualization
  
  return(forest_data)
}

forest_data <- prep_forest_data(all_res)

# Create comprehensive forest plot
comprehensive_forest <- ggplot(forest_data, aes(x = or, y = reorder(outcome_clean, or))) +
  geom_point(aes(color = significant), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = significant), height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red"), 
                     labels = c("FALSE" = "P ≥ 0.05", "TRUE" = "P < 0.05")) +
  scale_x_log10(breaks = c(0.5, 0.8, 1.0, 1.25, 1.5, 2.0),
                labels = c("0.5", "0.8", "1.0", "1.25", "1.5", "2.0")) +
  labs(
    title = "Mendelian Randomization Results: Endometriosis → Pregnancy Outcomes",
    subtitle = "Odds ratios (95% CI) from Inverse Variance Weighted analysis",
    x = "Odds Ratio (log scale)",
    y = "Pregnancy Outcome",
    color = "Statistical Significance"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom"
  )

# Save comprehensive forest plot
ggsave(
  filename = file.path(plots_dir, "comprehensive_forest_plot.png"),
  plot = comprehensive_forest,
  width = 12, height = 10, dpi = 300
)

### 8. EFFECT SIZE COMPARISON PLOT ##############################################
message("Creating method comparison plot...")

# Compare different MR methods
method_comparison <- all_res %>%
  filter(method %in% c("Inverse variance weighted", "Weighted median", "MR Egger")) %>%
  mutate(
    or = exp(b),
    ci_lower = exp(b - 1.96 * se),
    ci_upper = exp(b + 1.96 * se),
    method_clean = case_when(
      method == "Inverse variance weighted" ~ "IVW",
      method == "Weighted median" ~ "WM",
      method == "MR Egger" ~ "Egger",
      TRUE ~ method
    )
  ) %>%
  # Focus on top outcomes by effect size
  group_by(outcome_full) %>%
  summarise(max_effect = max(abs(b), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_effect)) %>%
  slice_head(n = 15) %>%
  pull(outcome_full) %>%
  {filter(all_res, outcome_full %in% .)} %>%
  filter(method %in% c("Inverse variance weighted", "Weighted median", "MR Egger")) %>%
  mutate(
    or = exp(b),
    ci_lower = exp(b - 1.96 * se),
    ci_upper = exp(b + 1.96 * se),
    method_clean = case_when(
      method == "Inverse variance weighted" ~ "IVW",
      method == "Weighted median" ~ "WM",
      method == "MR Egger" ~ "Egger",
      TRUE ~ method
    )
  )

comparison_plot <- ggplot(method_comparison, aes(x = or, y = reorder(outcome_full, or), color = method_clean)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), 
                 position = position_dodge(width = 0.6), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  scale_x_log10(breaks = c(0.5, 0.8, 1.0, 1.25, 1.5, 2.0),
                labels = c("0.5", "0.8", "1.0", "1.25", "1.5", "2.0")) +
  labs(
    title = "MR Method Comparison: Endometriosis → Pregnancy Outcomes",
    subtitle = "Odds ratios (95% CI) comparing different MR methods",
    x = "Odds Ratio (log scale)",
    y = "Pregnancy Outcome",
    color = "MR Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(plots_dir, "method_comparison_plot.png"),
  plot = comparison_plot,
  width = 12, height = 10, dpi = 300
)

### 9. COMBINED PLOTS FOR KEY OUTCOMES ##########################################
message("Creating combined plots for key outcomes...")

# Identify top significant outcomes
top_outcomes <- all_res %>%
  filter(method == "Inverse variance weighted", pval < 0.05) %>%
  arrange(pval) %>%
  slice_head(n = 6) %>%
  pull(outcome)

# Create combined plots for each top outcome
for (outcome_name in top_outcomes) {
  outcome_clean <- gsub("[^A-Za-z0-9]", "_", outcome_name)
  
  # Check if plots exist for this outcome
  scatter_plot <- scatter_plots[[outcome_name]]
  forest_plot <- forest_plots[[outcome_name]]
  loo_plot <- loo_plots[[outcome_name]]
  funnel_plot <- funnel_plots[[outcome_name]]
  
  # Create list of available plots
  available_plots <- list()
  if (!is.null(scatter_plot)) available_plots[["scatter"]] <- scatter_plot
  if (!is.null(forest_plot)) available_plots[["forest"]] <- forest_plot
  if (!is.null(loo_plot)) available_plots[["loo"]] <- loo_plot
  if (!is.null(funnel_plot)) available_plots[["funnel"]] <- funnel_plot
  
  # Create combined plot if we have at least 2 plots
  if (length(available_plots) >= 2) {
    combined_plot <- plot_grid(plotlist = available_plots, ncol = 2, align = "hv")
    
    # Add title
    title <- ggdraw() + 
      draw_label(paste("MR Analysis: Endometriosis →", unique(dat$outcome_full[dat$outcome == outcome_name])[1]),
                 fontface = 'bold', size = 16)
    
    final_combined <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
    
    ggsave(
      filename = file.path(plots_dir, paste0("combined_", outcome_clean, ".png")),
      plot = final_combined,
      width = 16, height = 12, dpi = 300
    )
  }
}

### 10. SUMMARY STATISTICS ######################################################
message("Generating plot summary...")

# Create summary of plots created
plot_summary <- data.frame(
  Plot_Type = c("Scatter plots", "Forest plots", "Leave-one-out plots", "Funnel plots", "Combined plots"),
  Count = c(length(scatter_plots), length(forest_plots), length(loo_plots), 
            length(funnel_plots), length(top_outcomes)),
  Description = c(
    "SNP effects scatter plots for significant outcomes",
    "Forest plots showing individual SNP and overall effects",
    "Leave-one-out sensitivity analysis plots",
    "Funnel plots for assessing publication bias",
    "Combined 4-panel plots for top significant outcomes"
  )
)

write.csv(plot_summary, file.path(plots_dir, "plot_summary.csv"), row.names = FALSE)

# Print summary
cat("\n========== PLOT GENERATION SUMMARY ==========\n")
print(plot_summary)
cat("\nAll plots saved to:", plots_dir, "\n")
cat("Key files:\n")
cat("- comprehensive_forest_plot.png: Overview of all results\n")
cat("- method_comparison_plot.png: Comparison of MR methods\n")
cat("- combined_[outcome].png: Multi-panel plots for top outcomes\n")
cat("- Individual plots: scatter_, forest_, loo_, funnel_ prefixes\n")

message("Plot generation completed successfully!")