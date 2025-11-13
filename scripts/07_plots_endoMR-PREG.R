#!/usr/bin/env Rscript
###############################################################################
# Focused Mendelian Randomization (MR) Plots for EndoMR-PREG Study           
# Focus: IVW significant results (p < 0.05) and birthweight outcomes          
###############################################################################

library(TwoSampleMR)
library(dplyr)
library(ggplot2)
library(here)
library(readr)
library(forcats)
library(scales)

# Define folders
results_dir <- here("results")
plots_dir <- here("plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Load harmonised data and MR results
harmonised_file <- file.path(results_dir, "harmonised_rahmioglu_bpo.csv")
dat <- read_csv(harmonised_file)
all_results <- read_csv(file.path(results_dir, "all_mr_methods.csv"))
single_snp_results <- read_csv(file.path(results_dir, "singlesnp_results.csv"))
leaveoneout_results <- read_csv(file.path(results_dir, "leaveoneout_snp_results.csv"))

# Outcomes with significant IVW p-values (< 0.05) - filtered to match 29 outcomes
significant_ivw_outcomes <- c(
  "lowapgar1",                               # Low Apgar score at 1 minute
  "el_cs",                                   # Elective caesarean section
  "finngen_R12_O15_PLAC_PRAEVIA",            # Placenta praevia
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",       # Placental abruption
  "pretb_all",                               # Preterm birth (all)
  "vpretb_all",                              # Very preterm birth
  "rup_memb",                                # Prelabour rupture of membranes
  "Postpartum_hemorrhage_filtered",          # Postpartum hemorrhage
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered" # PPH due to retained placenta
)

# Additional outcomes of interest: birthweight-related (filtered to 29 outcomes)
birthweight_outcomes <- c(
  "hbw_all",    # High birthweight (>4000 g)
  "lbw_all",    # Low birthweight (<2500 g)
  "sga"         # Small for gestational age
)

# Combined list of outcomes to plot
outcomes_to_plot <- c(significant_ivw_outcomes, birthweight_outcomes)

###############################################################################
# ENHANCED IVW FOREST PLOT FUNCTION
###############################################################################

plot_ivw_forest <- function(res, exposure_name = NULL, dark = TRUE,
                           title = "Endometriosis → Pregnancy Outcomes (IVW Method)") {
  
  # Load configuration if not already loaded
  if (!exists("OUTCOME_LABELS_FILTERED")) {
    source(here("config", "config.R"))
    source(here("config", "utils.R"))
  }
  
  df <- res %>%
    { if (!is.null(exposure_name)) filter(., .data$exposure == exposure_name) else . } %>%
    # Apply clean outcome labels using our utility function
    apply_outcome_labels() %>%
    add_outcome_metadata() %>%
    filter(grepl("^Inverse variance weighted", .data$method, ignore.case = TRUE)) %>%
    filter(!is.na(.data$b), !is.na(.data$se)) %>%
    mutate(
      OR = exp(.data$b),
      LCL = exp(.data$b - 1.96*.data$se),
      UCL = exp(.data$b + 1.96*.data$se),
      fdr = p.adjust(.data$pval, "fdr"),
      bonferroni = p.adjust(.data$pval, "bonferroni"),
      sig = case_when(
        bonferroni < 0.05 ~ "Bonferroni < 0.05",
        fdr < 0.05 ~ "FDR < 0.05", 
        .data$pval < 0.05 ~ "p < 0.05",
        TRUE ~ "NS"
      ),
      # Use the clean outcome labels (ensure character vector)
      outcome_lab = as.character(outcome_full)
    ) %>%
    arrange(.data$OR) %>%
    mutate(outcome_lab = fct_reorder(.data$outcome_lab, .data$OR))
  
  # Ensure nsnp column exists
  if (!"nsnp" %in% names(df) || all(is.na(df$nsnp))) df$nsnp <- 1L
  
  # Create base plot
  p <- ggplot(df, aes(.data$OR, .data$outcome_lab)) +
    geom_errorbarh(aes(xmin = .data$LCL, xmax = .data$UCL, color = .data$sig), 
                   height = 0.3, alpha = 0.8, size = 0.8) +
    geom_point(aes(size = .data$nsnp, fill = .data$sig), 
               shape = 21, stroke = 0.5, color = "black", alpha = 0.9) +
    geom_vline(xintercept = 1, linetype = "dashed", size = 0.6, color = "gray50") +
    scale_x_log10(
      breaks = c(0.50, 0.70, 1, 1.25, 1.50, 2.00, 3.00),
      labels = number_format(accuracy = 0.01),
      limits = c(0.4, max(df$UCL, na.rm = TRUE) * 1.1)
    ) +
    scale_size_continuous(
      range = c(2.5, 7), 
      name = "n SNPs",
      breaks = c(min(df$nsnp, na.rm = TRUE), median(df$nsnp, na.rm = TRUE), max(df$nsnp, na.rm = TRUE)),
      labels = function(x) round(x)
    ) +
    scale_fill_manual(
      values = c("NS" = "grey70", "p < 0.05" = "#2b8cbe", 
                "FDR < 0.05" = "#e6550d", "Bonferroni < 0.05" = "#d62728"),
      breaks = c("NS", "p < 0.05", "FDR < 0.05", "Bonferroni < 0.05"), 
      name = "Significance"
    ) +
    scale_color_manual(
      values = c("NS" = "grey70", "p < 0.05" = "#2b8cbe", 
                "FDR < 0.05" = "#e6550d", "Bonferroni < 0.05" = "#d62728"),
      breaks = c("NS", "p < 0.05", "FDR < 0.05", "Bonferroni < 0.05"), 
      name = "Significance"
    ) +
    labs(
      x = "Odds Ratio (95% CI)",
      y = "",
      title = title,
      subtitle = paste0("n = ", nrow(df), " pregnancy outcomes | ", 
                       sum(df$sig != "NS", na.rm = TRUE), " significant associations"),
      caption = paste0("Point size = number of SNPs | ",
                      "Red: Bonferroni significant | Orange: FDR significant | Blue: Nominally significant")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linetype = "solid", size = 0.3, color = "grey85"),
      panel.grid.major.x = element_line(linetype = "solid", size = 0.2, color = "grey90"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
      plot.caption = element_text(hjust = 0, size = 9, color = "gray50"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11, face = "bold"),
      plot.margin = margin(20, 20, 20, 20),
      strip.text = element_text(face = "bold")
    ) +
    guides(
      fill = guide_legend(override.aes = list(size = 4), title.position = "top"),
      color = "none",
      size = guide_legend(title.position = "top")
    )
  
  # Apply dark theme if requested
  if (dark) {
    p <- p + theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      text = element_text(color = "white"),
      axis.text = element_text(color = "white"),
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.major.y = element_line(color = "grey60"),
      plot.title = element_text(color = "white"),
      plot.subtitle = element_text(color = "grey80"),
      plot.caption = element_text(color = "grey70")
    )
  }
  
  return(p)
}

###############################################################################
# GENERATE ENHANCED IVW FOREST PLOT
###############################################################################

message("Creating enhanced IVW forest plot for all outcomes...")

# Create the enhanced forest plot
p_enhanced <- plot_ivw_forest(
  all_results, 
  exposure_name = "Endometriosis (Rahmioglu)", 
  dark = FALSE
)

# Save both light and dark versions
ggsave(
  filename = file.path(plots_dir, "ivw_forest_all_outcomes.png"),
  plot = p_enhanced,
  width = 10, height = 12, dpi = 300, bg = "white"
)

# Dark version
p_enhanced_dark <- plot_ivw_forest(
  all_results, 
  exposure_name = "Endometriosis (Rahmioglu)", 
  dark = TRUE
)

ggsave(
  filename = file.path(plots_dir, "ivw_forest_all_outcomes_dark.png"),
  plot = p_enhanced_dark,
  width = 10, height = 12, dpi = 300, bg = "black"
)

message("Saved enhanced IVW forest plots: ivw_forest_all_outcomes.png (light & dark versions)")

message("Generating scatter, forest, leave-one-out, and funnel plots for ", length(outcomes_to_plot), " outcomes...")

### SCATTER PLOTS ###
message("Creating scatter plots...")
for (outcome_name in outcomes_to_plot) {
  
  # Filter harmonised and results data for this outcome
  dat_outcome <- dat %>% filter(outcome == outcome_name)
  res_outcome <- all_results %>% filter(outcome == outcome_name)
  
  # Only proceed if multiple MR methods are available
  methods_present <- unique(res_outcome$method)
  if (nrow(dat_outcome) > 0 && length(methods_present) > 1) {
    
    # Keep only methods supported by scatter plot
    allowed_methods <- c("Inverse variance weighted", "MR Egger", "Weighted median", "Simple mode", "Weighted mode")
    res_outcome <- res_outcome %>% filter(method %in% allowed_methods)
    
    # Generate scatter plot
    p1 <- tryCatch({
      mr_scatter_plot(res_outcome, dat_outcome)
    }, error = function(e) {
      message("Scatter plot error for ", outcome_name, ": ", e$message)
      NULL
    })
    
    # Save if plot successful
    if (!is.null(p1) && length(p1) > 0) {
      ggsave(p1[[1]],
             filename = file.path(plots_dir, paste0("scatter_", outcome_name, ".png")),
             width = 7, height = 7, dpi = 300)
      message("Saved scatter plot for: ", outcome_name)
    }
  } else {
    message("Skipping scatter plot for ", outcome_name, ": not enough MR methods available.")
  }
}


### FOREST PLOTS ###
for (outcome_name in outcomes_to_plot) {
  dat_filtered <- dat %>% filter(outcome == outcome_name)
  
  if (nrow(dat_filtered) > 0) {
    p <- tryCatch({
      mr_forest_plot(mr_singlesnp(dat_filtered))
    }, error = function(e) {
      message("Forest plot error for ", outcome_name, ": ", e$message)
      NULL
    })
    
    if (!is.null(p) && inherits(p[[1]], "gg")) {
      ggsave(file.path(plots_dir, paste0("forest_", outcome_name, ".png")),
             p[[1]], width = 10, height = 8, dpi = 300)
      message("Saved forest plot for: ", outcome_name)
    }
  }
}

### LEAVE-ONE-OUT PLOTS ###
for (outcome_name in outcomes_to_plot) {
  dat_filtered <- dat %>% filter(outcome == outcome_name)
  
  if (nrow(dat_filtered) > 0) {
    loo <- tryCatch(mr_leaveoneout(dat_filtered), error = function(e) NULL)
    
    if (!is.null(loo) && nrow(loo) > 0) {
      p <- mr_leaveoneout_plot(loo)
      if (!is.null(p) && inherits(p[[1]], "gg")) {
        ggsave(file.path(plots_dir, paste0("loo_", outcome_name, ".png")),
               p[[1]], width = 10, height = 8, dpi = 300)
        message("Saved leave-one-out plot for: ", outcome_name)
      }
    }
  }
}

### FUNNEL PLOTS ###
for (outcome_name in outcomes_to_plot) {
  dat_filtered <- dat %>% filter(outcome == outcome_name)
  
  if (nrow(dat_filtered) > 0) {
    single_snp <- tryCatch(mr_singlesnp(dat_filtered), error = function(e) NULL)
    
    if (!is.null(single_snp) && nrow(single_snp) > 0) {
      p <- mr_funnel_plot(single_snp)
      if (!is.null(p) && inherits(p[[1]], "gg")) {
        ggsave(file.path(plots_dir, paste0("funnel_", outcome_name, ".png")),
               p[[1]], width = 7, height = 7, dpi = 300)
        message("Saved funnel plot for: ", outcome_name)
      }
    }
  }
}


###############################################################################
# v2 Forest plot summary for all IVW results from all_results
###############################################################################

library(metafor)
library(dplyr)

# Keep only IVW results and calculate OR and 95% CI
df <- all_results %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    OR = exp(b),
    CI_low = exp(b - 1.96 * se),
    CI_high = exp(b + 1.96 * se)
  )

# Clean outcome labels (forced renaming, even if outcome_clean already exists)
label_map <- c(
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered" = "PPH - retained placenta",
  "Postpartum_hemorrhage_filtered" = "Postpartum hemorrhage",
  "Early_bleeding_ending_in_live_birth_filtered" = "Early bleeding (live birth)",
  "Antepartum_bleeding_filtered" = "Antepartum bleeding",
  "Postpartum_hemorrhage_due_to_atony_filtered" = "PPH - uterine atony",
  "finngen_R12_N14_FEMALEINFERT" = "Female infertility",
  "finngen_R12_O15_PLAC_PRAEVIA" = "Placenta praevia",
  "finngen_R12_O99_PLACENTA_DISORDER" = "Other placental disorders",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR" = "Placental abruption",
  "el_cs" = "Elective caesarean section",
  "apgar1" = "Low Apgar score at 1 min",
  "em_cs" = "Emergency caesarean section",
  "apgar5" = "Low Apgar score at 5 min",
  "rup_memb" = "Premature rupture of membranes",
  "gest_all" = "Gestational age",
  "lga" = "Large for gestational age",
  "gest_subsamp" = "Gestational age (subset)",
  "rm_subsamp" = "Recurrent miscarriage",
  "hbw_all" = "High birthweight (>4000 g)",
  "postterm" = "Post-term birth",
  "ppd" = "Postpartum depression",
  "ectopic" = "Ectopic pregnancy",
  "nvp_subsamp" = "Severe nausea/vomiting (subset)",
  "nvp_all" = "Severe nausea/vomiting",
  "cs_all" = "Caesarean section",
  "miscarriage_subsamp" = "Miscarriage",
  "sga" = "Small for gestational age",
  "sm_subsamp" = "Single miscarriage",
  "induction" = "Induction of labour",
  "hg" = "Hyperemesis gravidarum",
  "apgar_score5" = "Apgar score at 5 min",
  "breastfeed_24" = "Breastfeeding ≥4 months",
  "exclusive_breastfeed" = "Exclusive breastfeeding",
  "anaemia" = "Anaemia in pregnancy",
  "breastfeed_init" = "Breastfeeding initiation",
  "breastfeed_stop" = "Breastfeeding cessation",
  "hdp_subsamp" = "Hypertensive disorders",
  "gdm_subsamp" = "Gestational diabetes",
  "nicu" = "NICU admission",
  "apgar_score1" = "Apgar score at 1 min",
  "gh_subsamp" = "Gestational hypertension",
  "preeclampsia_subsamp" = "Preeclampsia",
  "stillbirth_subsamp" = "Stillbirth",
  "lbw_all" = "Low birthweight (<2500 g)",
  "pretb_subsamp" = "Preterm birth (spontaneous)",
  "pretb_all" = "Preterm birth (any)",
  "vpretb_all" = "Very preterm birth (any)",
  "Early_bleeding_with_any_outcome_filtered" = "Early bleeding (any outcome)",
  "zbw_all" = "Z-score birthweight"
)

# Force named character vector (important !)
label_map <- setNames(as.character(label_map), names(label_map))

# Convert outcome to character
df$outcome <- as.character(df$outcome)

# Map to clean names using named vector
df$outcome_clean <- label_map[df$outcome]

# Replace NAs (unmatched) by original outcome
df$outcome_clean[is.na(df$outcome_clean)] <- df$outcome[is.na(df$outcome_clean)]

print(head(df$outcome_clean))
print(head(df$outcome))

# Arrange by descending OR
df <- df %>% arrange(desc(OR))

# Plot dimensions
n <- nrow(df)
plot_height <- 100 + 25 * n  # adjust as needed

# Output file
png("plots/forest_metafor_ivw_v2.png", width = 1400, height = plot_height, res = 120)

# Forest plot
forest(
  x = df$b,
  sei = df$se,
  slab = sprintf("  %s", df$outcome_clean),
  xlab = "Odds Ratio",
  annotate = FALSE,
  header = FALSE, 
  ilab = data.frame(
    sprintf("%.2f", df$OR),
    sprintf("(%.2f, %.2f)", df$CI_low, df$CI_high),
    sprintf("%.2e", df$pval),
    df$nsnp
  ),
  ilab.xpos = c(3.1, 3.9, 4.7, 5.4),
  pch = 16,
  atransf = exp,
  at = log(c(0.5, 1, 2)),
  xlim = c(-6, 7.5),
  rows = seq(n, 1),
  ylim = c(0, n + 4)
)
# Add column headers
text(-6, n + 2.5, "Outcome", font = 2, pos = 4)
text(3.5, n + 2.5, "OR", font = 2)
text(4.5, n + 2.5, "95% CI", font = 2)
text(5.5, n + 2.5, "p-value", font = 2)
text(6.2, n + 2.5, "SNPs", font = 2)

dev.off()




# Dimensions du plot
png("plots/forest_metafor_ivw_v2.png", width = 1600, height = plot_height, res = 120)

# Forest plot ajusté
forest(
  x = df$b,
  sei = df$se,
  slab = sprintf("  %s", df$outcome_clean),
  xlab = "Odds Ratio",
  annotate = FALSE,
  header = FALSE,
  ilab = data.frame(
    sprintf("%.2f", df$OR),
    sprintf("(%.2f, %.2f)", df$CI_low, df$CI_high),
    sprintf("%.2e", df$pval),
    df$nsnp
  ),
  ilab.xpos = c(3.8, 4.7, 5.6, 6.3),  
  pch = 16,
  atransf = exp,
  at = log(c(0.5, 1, 2, 3)),  
  xlim = c(-6, 7.5),         
  rows = seq(n, 1),
  ylim = c(0, n + 4),
  cex = 0.9
)

# En-têtes de colonnes
text(-6, n + 2.5, "Outcome", font = 2, pos = 4, cex = 0.95)
text(3.8, n + 2.5, "OR", font = 2, cex = 0.95)
text(4.7, n + 2.5, "95% CI", font = 2, cex = 0.95)
text(5.6, n + 2.5, "p-value", font = 2, cex = 0.95)
text(6.3, n + 2.5, "SNPs", font = 2, cex = 0.95)

dev.off()


###############################################################################
# Multi-Method Forest Plots for EndoMR-PREG Study
###############################################################################

library(ggplot2)
library(dplyr)
library(gridExtra)

# Define clinical outcome names for focused plot
top_clinical_outcomes <- c(
  "finngen_R12_O15_PLAC_PRAEVIA",
  "finngen_R12_N14_FEMALEINFERT", 
  "rup_memb",
  "el_cs",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",
  "apgar1",
  "pretb_all",
  "finngen_R12_O99_PLACENTA_DISORDER",
  "miscarriage_subsamp",
  "Postpartum_hemorrhage_filtered",
  "sga"
)

# Enhanced outcome name mapping for clinical interpretation
clinical_outcome_names <- c(
  "finngen_R12_O15_PLAC_PRAEVIA" = "Placenta praevia",
  "finngen_R12_N14_FEMALEINFERT" = "Female infertility", 
  "rup_memb" = "Premature rupture of membranes",
  "el_cs" = "Elective caesarean section",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR" = "Placental abruption",
  "apgar1" = "Low Apgar score (1 min)",
  "pretb_all" = "Preterm birth",
  "finngen_R12_O99_PLACENTA_DISORDER" = "Placental disorders",
  "miscarriage_subsamp" = "Miscarriage",
  "Postpartum_hemorrhage_filtered" = "Postpartum hemorrhage",
  "sga" = "Small for gestational age",
  "anaemia_preg_all" = "Anaemia in pregnancy",
  "Antepartum_bleeding" = "Antepartum bleeding",
  "apgar5" = "Low Apgar score (5 min)",
  "apgar_score1" = "Apgar score (1 min)",
  "apgar_score5" = "Apgar score (5 min)",
  "breastfeed_24" = "Breastfeeding ≥4 months",
  "breastfeed_init" = "Breastfeeding initiation",
  "breastfeed_stop" = "Breastfeeding cessation",
  "cs_all" = "Caesarean section",
  "Early_bleeding_ending_in_live_birth" = "Early bleeding (live birth)",
  "Early_bleeding_with_any_outcome" = "Early bleeding (any outcome)",
  "ectopic" = "Ectopic pregnancy",
  "em_cs" = "Emergency caesarean section",
  "exclusive_breastfeed" = "Exclusive breastfeeding",
  "gdm_subsamp" = "Gestational diabetes",
  "gest_all" = "Gestational age",
  "gest_subsamp" = "Gestational age (subset)",
  "gh_subsamp" = "Gestational hypertension",
  "hbw_all" = "High birthweight",
  "hdp_subsamp" = "Hypertensive disorders",
  "hg" = "Hyperemesis gravidarum",
  "induction" = "Induction of labour",
  "lbw_all" = "Low birthweight",
  "lga" = "Large for gestational age",
  "nicu" = "NICU admission",
  "nvp_all" = "Nausea/vomiting",
  "nvp_subsamp" = "Nausea/vomiting (subset)",
  "postterm" = "Post-term birth",
  "Postpartum_hemorrhage_due_to_atony" = "PPH (uterine atony)",
  "Postpartum_hemorrhage_due_to_retained_placenta" = "PPH (retained placenta)",
  "ppd" = "Postpartum depression",
  "preeclampsia_subsamp" = "Preeclampsia",
  "pretb_subsamp" = "Preterm birth (spontaneous)",
  "rm_subsamp" = "Recurrent miscarriage",
  "sm_subsamp" = "Single miscarriage",
  "stillbirth_subsamp" = "Stillbirth",
  "vpretb_all" = "Very preterm birth",
  "zbw_all" = "Birthweight z-score"
)

# Function to create multi-method forest plot
create_multimethod_forest_plot <- function(mr_results, outcomes_to_include = NULL, title_suffix = "") {
  
  # Filter data if specific outcomes requested
  if (!is.null(outcomes_to_include)) {
    plot_data <- mr_results %>% filter(outcome %in% outcomes_to_include)
  } else {
    plot_data <- mr_results
  }
  
  # Keep only main MR methods
  plot_data <- plot_data %>%
    filter(method %in% c("Inverse variance weighted", "MR Egger", "Weighted median")) %>%
    mutate(
      OR = exp(b),
      CI_lower = exp(b - 1.96 * se),
      CI_upper = exp(b + 1.96 * se),
      # Map to clinical names
      outcome_clean = ifelse(outcome %in% names(clinical_outcome_names),
                           clinical_outcome_names[outcome],
                           outcome),
      # Significance levels
      significance = case_when(
        pval < 0.001 ~ "p < 0.001 (Bonferroni)",
        pval < 0.05 ~ "p < 0.05",
        TRUE ~ "Not significant"
      )
    )
  
  # Order by IVW OR magnitude for readability
  ivw_order <- plot_data %>%
    filter(method == "Inverse variance weighted") %>%
    arrange(desc(OR)) %>%
    pull(outcome_clean)
  
  plot_data$outcome_clean <- factor(plot_data$outcome_clean, levels = ivw_order)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = OR, y = outcome_clean, color = method)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), 
                   height = 0.3, position = position_dodge(width = 0.5)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    scale_x_log10(
      breaks = c(0.5, 1, 1.5, 2, 3),
      labels = c("0.5", "1.0", "1.5", "2.0", "3.0")
    ) +
    scale_color_manual(
      values = c("Inverse variance weighted" = "#E31A1C", 
                 "MR Egger" = "#1F78B4", 
                 "Weighted median" = "#33A02C"),
      name = "MR Method"
    ) +
    labs(
      title = paste0("Endometriosis → Pregnancy Outcomes", title_suffix),
      subtitle = "Odds ratios with 95% confidence intervals",
      x = "Odds Ratio (log scale)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 11),
      legend.title = element_text(face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  return(p)
}

# Create forest plot for ALL outcomes
message("Creating comprehensive forest plot for all outcomes...")
forest_all <- create_multimethod_forest_plot(
  all_results, 
  title_suffix = " (All 49 Outcomes)"
)

ggsave(
  filename = file.path(plots_dir, "forest_multimethod_all_outcomes.png"),
  plot = forest_all,
  width = 12, height = 16, dpi = 300
)

message("Saved comprehensive forest plot: forest_multimethod_all_outcomes.png")

# Create forest plot for TOP 11 clinical outcomes
message("Creating focused forest plot for top 11 clinical outcomes...")
forest_clinical <- create_multimethod_forest_plot(
  all_results,
  outcomes_to_include = top_clinical_outcomes,
  title_suffix = " (Top 11 Clinical Outcomes)"
)

ggsave(
  filename = file.path(plots_dir, "forest_multimethod_clinical_top11.png"),
  plot = forest_clinical,
  width = 12, height = 8, dpi = 300
)

message("Saved clinical forest plot: forest_multimethod_clinical_top11.png")

message("Multi-method forest plots completed successfully!")

###############################################################################
# Selected Variables Forest Plot - EXACT SAME STYLE as original
# Based on ivw_forest_all_outcomes.png with metafor library
###############################################################################

# Selected variables (as specified)
selected_variables <- c(
  "finngen_R12_O15_PLAC_PRAEVIA",
  "rup_memb",
  "el_cs",
  "em_cs",
  "cs_all",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",
  "apgar1",
  "apgar5",
  "pretb_all",
  "finngen_R12_O99_PLACENTA_DISORDER",
  "miscarriage_subsamp",
  "Postpartum_hemorrhage_filtered",
  "sga",
  "Antepartum_bleeding_filtered",
  "Early_bleeding_ending_in_live_birth_filtered",
  "Early_bleeding_with_any_outcome_filtered",
  "finngen_R12_O15_PREG_ECTOP_filtered",
  "gdm_subsamp",
  "gh_subsamp",
  "hbw_all",
  "hdp_subsamp",
  "hg",
  "induction",
  "lbw_all",
  "lga",
  "nicu",
  "nvp_all",
  "nvp_subsamp",
  "postterm",
  "Postpartum_hemorrhage_due_to_atony_filtered",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered",
  "ppd",
  "preeclampsia_subsamp",
  "pretb_subsamp",
  "rm_subsamp",
  "sm_subsamp",
  "stillbirth_subsamp",
  "vpretb_all"
)

# Create selected variables forest plot using EXACT same style as original
message("Creating selected variables forest plot with exact original style...")

# Filter for IVW results and selected variables only
df_selected <- all_results %>%
  filter(method == "Inverse variance weighted") %>%
  filter(outcome %in% selected_variables) %>%
  mutate(
    OR = exp(b),
    CI_low = exp(b - 1.96 * se),
    CI_high = exp(b + 1.96 * se)
  )

# Labels for selected variables (exact mapping as specified)
selected_label_map <- c(
  "finngen_R12_O15_PLAC_PRAEVIA" = "Placenta praevia",
  "rup_memb" = "Premature rupture of membranes",
  "el_cs" = "Elective caesarean section",
  "em_cs" = "Emergency caesarean section",
  "cs_all" = "Caesarean section",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR" = "Placental abruption",
  "apgar1" = "Low Apgar score (1 min)",
  "apgar5" = "Low Apgar score (5 min)",
  "pretb_all" = "Preterm birth",
  "finngen_R12_O99_PLACENTA_DISORDER" = "Placental disorders",
  "miscarriage_subsamp" = "Miscarriage",
  "Postpartum_hemorrhage_filtered" = "Postpartum hemorrhage",
  "sga" = "Small for gestational age",
  "Antepartum_bleeding_filtered" = "Antepartum bleeding",
  "Early_bleeding_ending_in_live_birth_filtered" = "Early bleeding (live birth)",
  "Early_bleeding_with_any_outcome_filtered" = "Early bleeding (any outcome)",
  "finngen_R12_O15_PREG_ECTOP_filtered" = "Ectopic pregnancy",
  "gdm_subsamp" = "Gestational diabetes",
  "gh_subsamp" = "Gestational hypertension",
  "hbw_all" = "High birthweight",
  "hdp_subsamp" = "Hypertensive disorders",
  "hg" = "Hyperemesis gravidarum",
  "induction" = "Induction of labour",
  "lbw_all" = "Low birthweight",
  "lga" = "Large for gestational age",
  "nicu" = "NICU admission",
  "nvp_all" = "Nausea/vomiting",
  "nvp_subsamp" = "Nausea/vomiting (subset)",
  "postterm" = "Post-term birth",
  "Postpartum_hemorrhage_due_to_atony_filtered" = "PPH (uterine atony)",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered" = "PPH (retained placenta)",
  "ppd" = "Postpartum depression",
  "preeclampsia_subsamp" = "Preeclampsia",
  "pretb_subsamp" = "Preterm birth (spontaneous)",
  "rm_subsamp" = "Recurrent miscarriage",
  "sm_subsamp" = "Single miscarriage",
  "stillbirth_subsamp" = "Stillbirth",
  "vpretb_all" = "Very preterm birth"
)

# Force named character vector (important!)
selected_label_map <- setNames(as.character(selected_label_map), names(selected_label_map))

# Convert outcome to character
df_selected$outcome <- as.character(df_selected$outcome)

# Map to clean names using named vector
df_selected$outcome_clean <- selected_label_map[df_selected$outcome]

# Replace NAs (unmatched) by original outcome
df_selected$outcome_clean[is.na(df_selected$outcome_clean)] <- df_selected$outcome[is.na(df_selected$outcome_clean)]

# Arrange by descending OR (same as original)
df_selected <- df_selected %>% arrange(desc(OR))

# Plot dimensions
n_selected <- nrow(df_selected)
plot_height_selected <- 100 + 25 * n_selected

# Output file - EXACT same style as original
png("plots/ivw_forest_all_outcomes_selected_variables.png", width = 1600, height = plot_height_selected, res = 120)

# Forest plot with EXACT same parameters as original
forest(
  x = df_selected$b,
  sei = df_selected$se,
  slab = sprintf("  %s", df_selected$outcome_clean),
  xlab = "Odds Ratio",
  annotate = FALSE,
  header = FALSE,
  ilab = data.frame(
    sprintf("%.2f", df_selected$OR),
    sprintf("(%.2f, %.2f)", df_selected$CI_low, df_selected$CI_high),
    sprintf("%.2e", df_selected$pval),
    df_selected$nsnp
  ),
  ilab.xpos = c(3.8, 4.7, 5.6, 6.3),  # Same as original
  pch = 16,
  atransf = exp,
  at = log(c(0.5, 1, 2, 3)),  # Same as original
  xlim = c(-2.5, 7.5),         # Same as original
  rows = seq(n_selected, 1),
  ylim = c(0, n_selected + 4),
  cex = 0.9
)

# Column headers - EXACT same as original
text(-6, n_selected + 2.5, "Outcome", font = 2, pos = 4, cex = 0.95)
text(3.8, n_selected + 2.5, "OR", font = 2, cex = 0.95)
text(4.7, n_selected + 2.5, "95% CI", font = 2, cex = 0.95)
text(5.6, n_selected + 2.5, "p-value", font = 2, cex = 0.95)
text(6.3, n_selected + 2.5, "SNPs", font = 2, cex = 0.95)

dev.off()

message(paste0("✓ Created selected variables forest plot: ivw_forest_all_outcomes_selected_variables.png"))
message(paste0("   - Included ", n_selected, " selected variables"))
message("   - Used exact same style as original forest plot")

