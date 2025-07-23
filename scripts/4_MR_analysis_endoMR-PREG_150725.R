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
message("✔ Loaded ", nrow(dat), " rows across ", length(phenotypes), " phenotypes")
# View phenotypes
print(phenotypes)

### 3. HELPER FUNCTIONS ##########################################################
run_mr_methods <- function(data, methods) {
  res <- mr(data, method_list = methods)
  as.data.frame(res)
}
export_csv <- function(df, name) {
  write.csv(df, file.path(results_dir, paste0(name, ".csv")), row.names = FALSE)
}

### 4. MAIN MR ESTIMATES ########################################################
# 4.1 Inverse-Variance Weighted (IVW)
ivw_res <- run_mr_methods(dat, "mr_ivw")
print(ivw_res)

# 4.2 MR-Egger
egger_res <- run_mr_methods(dat, "mr_egger_regression")

# 4.3 Weighted median
wm_res   <- run_mr_methods(dat, "mr_weighted_median")

# 4.4 All default MR methods (IVW, Egger, WM, modes…)
#all_res  <- mr(dat)

# Export results
export_csv(ivw_res,  "ivw_results")
export_csv(egger_res,"egger_results")
export_csv(wm_res,   "weighted_median_results")
#export_csv(all_res,  "all_mr_methods")

### 5. SENSITIVITY ANALYSES ######################################################
# 5.1 Heterogeneity (Cochran’s Q)
het_res <- mr_heterogeneity(dat)
export_csv(het_res, "heterogeneity_results")

# 5.2 Horizontal pleiotropy (Egger intercept)
plt_res <- mr_pleiotropy_test(dat)
export_csv(plt_res, "pleiotropy_results")

# 5.3 Single-SNP effects
single_res <- mr_singlesnp(dat)
export_csv(single_res, "singlesnp_results")

# 5.4 Leave-one-out by SNP
loo_snp_res <- mr_leaveoneout(dat)
export_csv(loo_snp_res, "leaveoneout_snp_results")

# Quick prints
print(head(het_res))
print(head(plt_res))
print(head(single_res))
print(head(loo_snp_res))

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

### EXCLUDED SNPs? -----------------------------------------------------
# Harmonise exposure and outcome data
harm <- harmonise_data(exposure_dat, outcome_dat)
# Compare with your original list of 41 SNPs
setdiff(exposure_dat$SNP, harm$SNP)

