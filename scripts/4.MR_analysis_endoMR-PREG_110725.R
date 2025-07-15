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
}

# Paths
project_dir <- here::here()
results_dir <- file.path(project_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
harm_file   <- file.path(results_dir, "harmonised_rahmioglu_bpo.csv")

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

### 6. MR-PRESSO OUTLIER DETECTION ################################################
if (nrow(dat) >= 4) {
  presso <- tryCatch(
    mr_presso(
      BetaOutcome     = "beta.outcome",
      BetaExposure    = "beta.exposure",
      SdOutcome       = "se.outcome",
      SdExposure      = "se.exposure",
      data            = dat,
      NbDistribution  = 10000,
      SignifThreshold = 0.05
    ),
    error = function(e) {
      warning("MR-PRESSO failed: ", e$message)
      NULL
    }
  )
  if (!is.null(presso)) print(presso)
} else {
  message("👉 Skipping MR-PRESSO: fewer than 4 instruments")
}

### 7. PLOTTING ###################################################################
# 7.1 Scatter plot (IVW vs Egger vs modes)
p_scatter <- mr_scatter_plot(all_res, dat)
ggsave(plot = p_scatter[[1]], filename = file.path(results_dir, "scatter_ivw_egger.png"),
       width = 7, height = 7)

# 7.2 Forest plot for single-SNP
p_forest <- mr_forest_plot(single_res)
ggsave(plot = p_forest[[1]], filename = file.path(results_dir, "forest_singlesnp.png"),
       width = 7, height = 7)

# 7.3 Funnel plot for single-SNP
p_funnel <- mr_funnel_plot(single_res)
ggsave(plot = p_funnel[[1]], filename = file.path(results_dir, "funnel_singlesnp.png"),
       width = 7, height = 7)

# 7.4 Leave-one-out SNP plot
p_loo <- mr_leaveoneout_plot(loo_snp_res)
ggsave(plot = p_loo[[1]], filename = file.path(results_dir, "leaveoneout_plot.png"),
       width = 7, height = 7)

### 8. LEAVE-ONE-OUT BY COHORT ###################################################
# Requires a file 'stu_out_dat.txt' with columns: SNP, beta, se, effect_allele, other_allele, eaf, Phenotype, study
if (file.exists("stu_out_dat.txt")) {
  mr_data <- readr::read_delim("stu_out_dat.txt", delim = "\t", col_types = cols())
  cohorts  <- unique(mr_data$study)
  
  loo_cohort <- lapply(cohorts, function(coh) {
    df   <- dplyr::filter(mr_data, study != coh)
    exp  <- df %>% rename(beta.exposure = beta, se.exposure = se) %>% mutate(id.exposure = unique(df$Phenotype))
    out  <- df %>% rename(beta.outcome = beta, se.outcome = se, id.outcome = Phenotype)
    hmd  <- harmonise_data(exp, out, action = 2)
    res  <- mr(hmd, method_list = "mr_ivw")
    res$left_out_cohort <- coh
    res
  })
  loo_cohort_df <- dplyr::bind_rows(loo_cohort)
  export_csv(loo_cohort_df, "leaveoneout_by_cohort")
}

message("✅ Script finished. All results are saved in ", results_dir, "/")
