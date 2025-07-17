###############################################################################
# Harmonisation of Rahmioglu Endometriosis SNPs ↔ Bad Pregnancy Outcomes       #
#   (MR-PREG adverse outcomes & Westergaard PPH + FinnGen R12 outcomes)       #
# Output: harmonised_rahmioglu_bpo.csv                                        #
###############################################################################

### 1. SETUP ##################################################################

# Load required packages (install if missing)
pkgs <- c("TwoSampleMR", "dplyr", "data.table", "here")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
}

# Define directories
project_dir <- here::here()         
data_dir    <- here::here("data")
results_dir <- here::here("results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

### 2. EXPOSURE DATA ##########################################################

# Assume `clumped` is your prepared endometriosis exposure data frame
exposure_dat <- clumped
stopifnot("id.exposure" %in% colnames(exposure_dat))

### 3. OUTCOME DATA ###########################################################

# Combine the MR-ready outcomes you created in the previous script:
#   * mr_preg      (MR-PREG adverse pregnancy outcomes)
#   * res_pph list (Westergaard PPH outcomes)
#   * res_fg list  (FinnGen R12 outcomes)

# Reconstruct a single outcome data.frame:
outcomes_pph <- bind_rows(lapply(res_pph, `[[`, "mr"))
outcomes_fg  <- bind_rows(lapply(res_fg,  `[[`, "mr"))

outcome_dat <- bind_rows(
  mr_preg_dat,
  outcomes_pph,
  outcomes_fg
)
stopifnot(all(c("id.outcome", "outcome") %in% colnames(outcome_dat)))

### 4. HARMONISATION ##########################################################

dat <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat  = outcome_dat,
  action       = 2
)

# Quick summary of harmonisation
table(dat$action)        # alignment status
table(dat$remove)        # how many SNPs dropped
table(dat$palindromic)   # palindromic flags
table(dat$outcome)       # SNP counts per outcome

# Optional: inspect first rows
print(head(dat))
print(str(dat))
message("Total harmonised rows: ", nrow(dat))

### 5. SAVE THE HARMONISED DATA ##############################################

out_file <- file.path(results_dir, "harmonised_rahmioglu_bpo.csv")
write.csv(dat, out_file, row.names = FALSE)
message("✔ Harmonised data saved to ", out_file)

### 6. N° PATIENTS PER SELECTED SNP ##############################################
snp_summary_41 <- final_instruments %>%
  select(
    SNP,
    beta = beta.exposure,
    se   = se.exposure,
    pval = pval.exposure,
    eaf  = eaf.exposure,
    n    = samplesize.exposure
  )

View(snp_summary_41)     
# OU
print(snp_summary_41)    # Console
