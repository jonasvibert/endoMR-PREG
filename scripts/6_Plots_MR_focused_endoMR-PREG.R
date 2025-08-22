#!/usr/bin/env Rscript
###############################################################################
# Focused MR Plots for endoMR-PREG Study                                      
# Focus: IVW significant results (p < 0.05) and birthweight outcomes          
###############################################################################

library(TwoSampleMR)
library(dplyr)
library(ggplot2)
library(here)

# Paths
results_dir <- here::here("results")
plots_dir <- here::here("plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Load harmonised data and results
harm_file <- here::here("results", "harmonised_rahmioglu_bpo.csv")
dat <- read.csv(harm_file)
all_res <- read.csv(file.path(results_dir, "all_mr_methods.csv"))
single_res <- read.csv(file.path(results_dir, "singlesnp_results.csv"))
loo_res <- read.csv(file.path(results_dir, "leaveoneout_snp_results.csv"))

# IVW significant outcomes (p < 0.05)
ivw_significant <- c(
  "apgar1",           # Apgar score at 1 min (p = 0.029)
  "el_cs",            # Elective caesarean section (p = 0.023)
  "finngen_R12_N14_FEMALEINFERT",  # Female infertility (p = 6.07e-22)
  "finngen_R12_O15_PLAC_PRAEVIA",  # Placenta praevia (p = 1.46e-06)
  "finngen_R12_O15_PLAC_PREMAT_SEPAR", # Placental abruption (p = 0.031)
  "pretb_all",        # Preterm birth (all) (p = 0.007)
  "pretb_subsamp",    # Preterm birth (subsample) (p = 0.029)
  "vpretb_all",       # Very preterm birth (all) (p = 0.048)
  "rup_memb",         # Premature rupture of membranes (p = 0.025)
  "Postpartum_hemorrhage_filtered",    # Postpartum hemorrhage (p = 0.017)
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered"  # PPH - retained placenta (p = 0.003)
)

# Birthweight outcomes
birthweight_outcomes <- c(
  "hbw_all",    # High birthweight (>4000 g)
  "lbw_all",    # Low birthweight (<2500 g)
  "zbw_all"     # Z-score birthweight (all)
)

# Combined list for plotting
outcomes_to_plot <- c(ivw_significant, birthweight_outcomes)

message("Creating plots for ", length(outcomes_to_plot), " outcomes...")

### SCATTER PLOTS ###
message("Creating scatter plots...")
for (outcome_name in outcomes_to_plot) {
  # Filter data
  dat_outcome <- dat[dat$outcome == outcome_name, ]
  res_outcome <- all_res[all_res$outcome == outcome_name, ]
  
  if (nrow(dat_outcome) > 0 && nrow(res_outcome) > 0) {
    # Create scatter plot
    p1 <- mr_scatter_plot(res_outcome, dat_outcome)
    if (length(p1) > 0) {
      ggsave(p1[[1]], 
             file = file.path(plots_dir, paste0("scatter_", outcome_name, ".png")), 
             width = 7, height = 7)
      message("Saved scatter plot for ", outcome_name)
    }
  }
}

### FOREST PLOTS ###
message("Creating forest plots...")
for (outcome_name in outcomes_to_plot) {
  dat_outcome <- dat[dat$outcome == outcome_name, ]
  
  if (nrow(dat_outcome) > 0) {
    single_res <- mr_singlesnp(dat_outcome)
    p2 <- mr_forest_plot(single_res)
    
    if (length(p2) > 0 && inherits(p2[[1]], "gg")) {
      ggsave(
        filename = file.path(plots_dir, paste0("forest_", outcome_name, ".png")),
        plot = p2[[1]],
        width = 10, height = 8, dpi = 300
      )
      message("Saved forest plot for ", outcome_name)
    } else {
      message("No valid forest plot for ", outcome_name)
    }
  } else {
    message("No data for forest plot: ", outcome_name)
  }
}

### LEAVE-ONE-OUT PLOTS ###
message("Creating leave-one-out plots...")
for (outcome_name in outcomes_to_plot) {
  dat_outcome <- dat[dat$outcome == outcome_name, ]
  
  if (nrow(dat_outcome) > 0) {
    # Générer les résultats leave-one-out pour cet outcome
    loo_res <- tryCatch({
      mr_leaveoneout(dat_outcome)
    }, error = function(e) {
      message("Error in leave-one-out for ", outcome_name, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(loo_res) && nrow(loo_res) > 0) {
      p3 <- mr_leaveoneout_plot(loo_res)
      
      if (length(p3) > 0 && inherits(p3[[1]], "gg")) {
        ggsave(
          filename = file.path(plots_dir, paste0("loo_", outcome_name, ".png")),
          plot = p3[[1]],
          width = 10, height = 8,
          dpi = 300
        )
        message("Saved leave-one-out plot for ", outcome_name)
      } else {
        message("No valid leave-one-out plot for ", outcome_name)
      }
    } else {
      message("No leave-one-out results for ", outcome_name)
    }
  } else {
    message("No data for leave-one-out plot: ", outcome_name)
  }
}


### FUNNEL PLOTS ###
message("Creating funnel plots...")
for (outcome_name in outcomes_to_plot) {
  dat_outcome <- dat[dat$outcome == outcome_name, ]
  
  if (nrow(dat_outcome) > 0) {
    # Générer les résultats single-SNP pour cet outcome
    single_res <- tryCatch({
      mr_singlesnp(dat_outcome)
    }, error = function(e) {
      message("Error in single SNP results for ", outcome_name, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(single_res) && nrow(single_res) > 0) {
      # Créer le funnel plot
      p4 <- mr_funnel_plot(single_res)
      
      if (length(p4) > 0 && inherits(p4[[1]], "gg")) {
        ggsave(
          filename = file.path(plots_dir, paste0("funnel_", outcome_name, ".png")),
          plot = p4[[1]],
          width = 7, height = 7,
          dpi = 300
        )
        message("Saved funnel plot for ", outcome_name)
      } else {
        message("No valid funnel plot for ", outcome_name)
      }
    } else {
      message("No single SNP results for funnel plot: ", outcome_name)
    }
  } else {
    message("No data for funnel plot: ", outcome_name)
  }
}


### SUMMARY FOREST PLOT FOR ALL IVW OUTCOMES ###
message("Creating summary forest plot for all IVW outcomes...")

# Get all IVW results and prepare data
ivw_summary <- all_res %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    or = exp(b),
    ci_lower = exp(b - 1.96 * se),
    ci_upper = exp(b + 1.96 * se),
    significant = pval < 0.05,
    pval_text = case_when(
      pval < 0.001 ~ paste0("P = ", formatC(pval, format = "e", digits = 1)),
      pval < 0.01 ~ paste0("P = ", format(round(pval, 3), nsmall = 3)),
      TRUE ~ paste0("P = ", format(round(pval, 2), nsmall = 2))
    ),
    or_ci_text = paste0(format(round(or, 2), nsmall = 2), 
                       " (", format(round(ci_lower, 2), nsmall = 2), 
                       "–", format(round(ci_upper, 2), nsmall = 2), ")"),
    outcome_clean = case_when(
      outcome == "apgar1" ~ "Apgar score at 1 min",
      outcome == "el_cs" ~ "Elective caesarean section", 
      outcome == "finngen_R12_N14_FEMALEINFERT" ~ "Female infertility",
      outcome == "finngen_R12_O15_PLAC_PRAEVIA" ~ "Placenta praevia",
      outcome == "finngen_R12_O15_PLAC_PREMAT_SEPAR" ~ "Placental abruption",
      outcome == "pretb_all" ~ "Preterm birth (all)",
      outcome == "pretb_subsamp" ~ "Preterm birth (subsample)",
      outcome == "vpretb_all" ~ "Very preterm birth (all)",
      outcome == "rup_memb" ~ "Premature rupture of membranes",
      outcome == "Postpartum_hemorrhage_filtered" ~ "Postpartum hemorrhage",
      outcome == "Postpartum_hemorrhage_due_to_retained_placenta_filtered" ~ "PPH - retained placenta",
      outcome == "Postpartum_hemorrhage_due_to_atony" ~ "PPH - atony",
      outcome == "Antepartum_bleeding" ~ "Antepartum bleeding",
      outcome == "Early_bleeding_with_any_outcome" ~ "Early bleeding (any)",
      outcome == "Early_bleeding_ending_in_live_birth" ~ "Early bleeding (live birth)",
      outcome == "hbw_all" ~ "High birthweight (>4000 g)",
      outcome == "lbw_all" ~ "Low birthweight (<2500 g)",
      outcome == "zbw_all" ~ "Z-score birthweight",
      outcome == "Anaemia in pregnancy (all)" ~ "Anaemia in pregnancy",
      outcome == "Gestational age (all)" ~ "Gestational age",
      outcome == "Apgar score at 5 min" ~ "Apgar score at 5 min",
      outcome == "Low birthweight (<2500 g)" ~ "Low birthweight",
      outcome == "High birthweight (>4000 g)" ~ "High birthweight",
      outcome == "Z-score birthweight (all)" ~ "Z-score birthweight",
      TRUE ~ outcome
    ),
    effect_direction = case_when(
      or > 1 & significant ~ "Risk increase",
      or < 1 & significant ~ "Risk decrease", 
      TRUE ~ "Non-significant"
    )
  ) %>%
  arrange(pval) %>%
  slice_head(n = 25)  # Top 25 results

# Create enhanced summary forest plot
summary_forest <- ggplot(ivw_summary, aes(x = or, y = reorder(outcome_clean, pval))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#CC79A7", alpha = 0.8, size = 0.8) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = effect_direction), 
                 height = 0.4, size = 0.8, alpha = 0.8) +
  geom_point(aes(color = effect_direction, shape = significant), size = 3.5) +
  scale_color_manual(
    values = c("Risk increase" = "#D55E00", "Risk decrease" = "#0072B2", "Non-significant" = "gray50"),
    name = "Effect"
  ) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    name = "Significant"
  ) +
  scale_x_log10(
    breaks = c(0.1, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10, 20, 50, 100),
    labels = c("0.1", "0.3", "0.5", "0.7", "1.0", "1.5", "2.0", "3.0", "5.0", "10", "20", "50", "100")
  ) +
  labs(
    title = "Mendelian Randomization Results: Endometriosis → Pregnancy Outcomes",
    subtitle = "Inverse Variance Weighted (IVW) odds ratios with 95% confidence intervals",
    x = "Odds Ratio (log scale)",
    y = NULL,
    caption = "Ordered by statistical significance (most significant at top)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray95", size = 0.3),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_legend(title = "Effect Direction", override.aes = list(size = 3)),
    shape = guide_legend(title = "P < 0.05", override.aes = list(size = 3))
  )

# Save summary forest plot
ggsave(summary_forest,
       file = file.path(plots_dir, "summary_forest_all_IVW.png"),
       width = 12, height = 10, dpi = 300)

ggsave(summary_forest,
       file = file.path(plots_dir, "summary_forest_all_IVW.pdf"),
       width = 12, height = 10)

message("Saved summary forest plot: summary_forest_all_IVW.png")

message("Plot generation completed!")
message("Plots saved to: ", plots_dir)