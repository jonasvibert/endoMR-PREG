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

### 2. LOAD and sanity checks ####################################################
dat <- data.table::fread(harm_file)
stopifnot(all(c("id.exposure", "id.outcome", "beta.exposure", "beta.outcome") %in% colnames(dat)))

phenotypes <- unique(dat$outcome)

### 2b. OUTCOME LABELS ###########################################################
# Dictionary of human-readable outcome names
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

phenotypes <- unique(dat$outcome_full)
# View phenotypes
print(phenotypes)

### 3. HELPER FUNCTIONS ##########################################################
run_mr_methods <- function(data, methods) {
  res <- mr(data, method_list = methods)
  as.data.frame(res)
}
export_csv <- function(df, name) {
  # Add labels if outcome_full is missing
  if ("outcome" %in% colnames(df) && !"outcome_full" %in% colnames(df)) {
    df <- merge(df, labels_df, by = "outcome", all.x = TRUE)
    df$outcome_full[is.na(df$outcome_full)] <- df$outcome[is.na(df$outcome_full)]
  }
  
  # Drop unnecessary columns
  drop_cols <- c("id.exposure", "id.outcome", "outcome", "exposure")
  df <- df[, !(names(df) %in% drop_cols), drop = FALSE]
  
  # Rename outcome_full → outcome
  if ("outcome_full" %in% colnames(df)) {
    names(df)[names(df) == "outcome_full"] <- "outcome"
    col_order <- c("outcome", setdiff(names(df), "outcome"))
    df <- df[, col_order, drop = FALSE]
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

### 5. SENSITIVITY ANALYSES ######################################################
# 5.1 Heterogeneity (Cochran’s Q)
het_res <- mr_heterogeneity(dat)

# 5.2 Horizontal pleiotropy (Egger intercept)
plt_res <- mr_pleiotropy_test(dat)

# 5.3 Single-SNP effects
single_res <- mr_singlesnp(dat)

# 5.4 Leave-one-out by SNP
loo_snp_res <- mr_leaveoneout(dat)

# Export
export_csv(het_res, "heterogeneity_results")
export_csv(plt_res, "pleiotropy_results")
export_csv(single_res, "singlesnp_results")
export_csv(loo_snp_res, "leaveoneout_snp_results")

### 6. MR-PRESSO OUTLIER DETECTION

# 1. load your IVW results & harmonised data
ivw_res <- read.csv(file.path(results_dir, "ivw_results.csv"), stringsAsFactors = FALSE)
dat      <- data.table::fread(harm_file) %>% as.data.frame()

# 2. pick only the outcomes with a significant IVW (p < 0.05)
sig_outcomes <- ivw_res %>%
  filter(method == "Inverse variance weighted", pval < 0.05) %>%
  pull(id.outcome) %>%
  unique()

# 3. run MR-PRESSO on each, with fewer permutations
presso_results <- list()
for(o in sig_outcomes) {
  dat_sub <- dat[dat$id.outcome == o, ]
  if(nrow(dat_sub) < 4) {
    message("Skipping ", o, ": fewer than 4 instruments")
    next
  }
  message("Running MR-PRESSO on: ", o)
  presso_results[[o]] <- tryCatch(
    mr_presso(
      BetaOutcome     = "beta.outcome",
      BetaExposure    = "beta.exposure",
      SdOutcome       = "se.outcome",
      SdExposure      = "se.exposure",
      data            = dat_sub,
      NbDistribution  = 1000,      # only 1k permutations
      SignifThreshold = 0.05
    ),
    error = function(e) {
      warning("MR-PRESSO failed for ", o, ": ", e$message)
      NULL
    }
  )
}

# 4. inspect
presso_results
  
# Export a compact summary for each outcome
for (o in names(presso_results)) {
  res <- presso_results[[o]]
  if (!is.null(res)) {
    export_csv(as.data.frame(res$`Main`), paste0("presso_", o, "_main"))
    if (!is.null(res$`OutlierTest`)) {
      export_csv(as.data.frame(res$`OutlierTest`), paste0("presso_", o, "_outliers"))
    }
  }
}
### EXCLUDED SNPs? -----------------------------------------------------
# Compare with your original list of 41 SNPs
setdiff(exposure_dat$SNP, dat$SNP)

