#!/usr/bin/env Rscript
###############################################################################
# Two‐Sample MR analysis                                                       
# Input  : results/harmonised_rahmioglu_bpo.csv                                 
# Output : CSVs and plots for IVW, Egger, WM, heterogeneity, pleiotropy, LOO,  
#          cohort leave-one-out                                               
###############################################################################

### 1. SETUP and package loading ################################################
pkgs <- c("TwoSampleMR", "MRPRESSO", "dplyr", "ggplot2", "here", "readr", "data.table")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "MRPRESSO") devtools::install_github("rondolab/MR-PRESSO")
    else install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
  library(here)
}

# Paths
results_dir <- here::here("results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
harm_file   <- here::here("results", "harmonised_rahmioglu_bpo.csv")
path_outcome <- here::here("data", "OUTCOME_MR-PREG", "stu_out_dat.txt")

### 2. LOAD and sanity checks ####################################################
dat <- data.table::fread(harm_file)
stopifnot(all(c("id.exposure", "beta.exposure", "beta.outcome") %in% colnames(dat)))

phenotypes <- unique(dat$outcome)

### 2b. OUTCOME LABELS ###########################################################
# Dictionary of human-readable outcome names
outcome_labels <- c(
  pretb_all      = "Preterm birth (any)",
  el_cs          = "Elective caesarean section",
  rup_memb       = "Premature rupture of membranes",
  pretb_subsamp  = "Preterm birth (spontaneous)",
  apgar1         = "Apgar score at 1 min",
  vpretb_all     = "Very preterm birth (any)",
  em_cs          = "Emergency caesarean section",
  lowapgar1      = "Low Apgar score at 1 min",
  pe_subsamp     = "Preeclampsia",
  lbw_all        = "Low birthweight (<2500 g)",
  ga_all         = "Gestational age",
  gh_subsamp     = "Gestational hypertension",
  ga_subsamp     = "Gestational age (subset)",
  lga            = "Large for gestational age",
  zbw_all        = "Z-score birthweight",
  nvp_sev_subsamp= "Severe nausea/vomiting (subset)",
  nvp_sev_all    = "Severe nausea/vomiting",
  sb_subsamp     = "Stillbirth",
  gdm_subsamp    = "Gestational diabetes",
  apgar5         = "Apgar score at 5 min",
  hdp_subsamp    = "Hypertensive disorders of pregnancy",
  r_misc_subsamp = "Recurrent miscarriage",
  bf_dur_4c      = "Breastfeeding ≥4 months",
  depr_subsamp   = "Postpartum depression",
  hbw_all        = "High birthweight (>4000 g)",
  nicu           = "NICU admission",
  lowapgar5      = "Low Apgar score at 5 min",
  bf_sus         = "Breastfeeding cessation",
  posttb_all     = "Post-term birth",
  misc_subsamp   = "Miscarriage",
  bf_ini         = "Breastfeeding initiation",
  bf_est         = "Exclusive breastfeeding",
  s_misc_subsamp = "Single miscarriage",
  cs             = "Caesarean section",
  hyp            = "Hyperemesis gravidarum",
  anaemia_preg_all = "Anaemia in pregnancy",
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

# Add a readable outcome column
# Transformer le dictionnaire en data.frame
labels_df <- data.frame(
  outcome = names(outcome_labels),
  outcome_full = unname(outcome_labels),
  stringsAsFactors = FALSE
)

# Fusionner sur outcome
dat <- merge(dat, labels_df, by = "outcome", all.x = TRUE)

# Si pas de correspondance, fallback au brut
dat$outcome_full[is.na(dat$outcome_full)] <- dat$outcome[is.na(dat$outcome_full)]
dat$id.outcome <- dat$outcome

# Nettoyage : exclure les lignes problématiques
dat <- dat %>%
  dplyr::filter(
    !is.na(beta.outcome),
    !is.na(se.outcome),
    !is.infinite(se.outcome),
    !is.na(pval.outcome),
    se.outcome > 0
  )

phenotypes <- unique(dat$outcome_full)
# View phenotypes
print(phenotypes)

### 3. HELPER FUNCTIONS ##########################################################
run_mr_methods <- function(data, methods) {
  res <- mr(data, method_list = methods)
  as.data.frame(res)
}
export_csv <- function(df, name) {
  if ("outcome" %in% colnames(df) && !"outcome_full" %in% colnames(df)) {
    df <- merge(df, labels_df, by = "outcome", all.x = TRUE)
    df$outcome_full[is.na(df$outcome_full)] <- df$outcome[is.na(df$outcome_full)]
  }
  
  # Format numeric columns: round normally, use scientific if very small
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], function(x) {
    sapply(x, function(val) {
      if (is.na(val)) return(NA)
      if (abs(val) < 0.001) format(val, scientific = TRUE, digits = 3)
      else round(val, 3)
    })
  })
  
  # Export
  write.csv(df, file.path(results_dir, paste0(name, ".csv")), row.names = FALSE)
}

### 4. MAIN MR ESTIMATES ########################################################
# 4.1 Inverse-Variance Weighted (IVW)
ivw_res <- run_mr_methods(dat, "mr_ivw")
ivw_res_raw <- ivw_res

# 4.2 MR-Egger
egger_res <- run_mr_methods(dat, "mr_egger_regression")

# 4.3 Weighted median
wm_res   <- run_mr_methods(dat, "mr_weighted_median")

# 4.4 All default MR methods (IVW, Egger, WM, modes…)
all_res  <- mr(dat)

# Export results
export_csv(ivw_res,  "ivw_results")
export_csv(egger_res,"egger_results")
export_csv(wm_res,   "weighted_median_results")
export_csv(all_res,  "all_mr_methods")
