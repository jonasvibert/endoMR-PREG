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

# Outcomes with significant IVW p-values (< 0.05)
significant_ivw_outcomes <- c(
  "apgar1",                                  # Apgar score at 1 minute
  "el_cs",                                   # Elective caesarean section
  "finngen_R12_N14_FEMALEINFERT",            # Female infertility
  "finngen_R12_O15_PLAC_PRAEVIA",            # Placenta praevia
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",       # Placental abruption
  "pretb_all",                               # Preterm birth (all)
  "pretb_subsamp",                           # Preterm birth (subsample)
  "vpretb_all",                              # Very preterm birth
  "rup_memb",                                # Prelabour rupture of membranes
  "Postpartum_hemorrhage_filtered",          # Postpartum hemorrhage
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered" # PPH due to retained placenta
)

# Additional outcomes of interest: birthweight-related
birthweight_outcomes <- c(
  "hbw_all",    # High birthweight (>4000 g)
  "lbw_all",    # Low birthweight (<2500 g)
  "zbw_all"     # Birthweight z-score
)

# Combined list of outcomes to plot
outcomes_to_plot <- c(significant_ivw_outcomes, birthweight_outcomes)

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
  "gest_all" = "Gestational age (all)",
  "lga" = "Large for gestational age",
  "gest_subsamp" = "Gestational age (subsample)",
  "rm_subsamp" = "Recurrent miscarriage (subsample)",
  "hbw_all" = "High birthweight (>4000 g)",
  "postterm" = "Post-term birth",
  "ppd" = "Postpartum depression",
  "ectopic" = "Ectopic pregnancy",
  "nvp_subsamp" = "Severe nausea/vomiting (subsample)",
  "nvp_all" = "Severe nausea/vomiting (all)",
  "cs_all" = "Caesarean section (all)",
  "miscarriage_subsamp" = "Miscarriage (subsample)",
  "sga" = "Small for gestational age",
  "sm_subsamp" = "Single miscarriage (subsample)",
  "induction" = "Induction of labour",
  "hg" = "Hyperemesis gravidarum",
  "apgar_score5" = "Apgar score at 5 min",
  "breastfeed_24" = "Breastfeeding ≥4 months",
  "exclusive_breastfeed" = "Exclusive breastfeeding",
  "anaemia" = "Anaemia in pregnancy",
  "breastfeed_init" = "Breastfeeding initiation",
  "breastfeed_stop" = "Breastfeeding cessation",
  "hdp_subsamp" = "Hypertensive disorders (subsample)",
  "gdm_subsamp" = "Gestational diabetes (subsample)",
  "nicu" = "NICU admission",
  "apgar_score1" = "Apgar score at 1 min",
  "gh_subsamp" = "Gestational hypertension (subsample)",
  "preeclampsia_subsamp" = "Preeclampsia (subsample)",
  "stillbirth_subsamp" = "Stillbirth (subsample)",
  "lbw_all" = "Low birthweight (<2500 g)",
  "pretb_subsamp" = "Preterm birth (subsample)",
  "pretb_all" = "Preterm birth (all)",
  "vpretb_all" = "Very preterm birth (all)",
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

