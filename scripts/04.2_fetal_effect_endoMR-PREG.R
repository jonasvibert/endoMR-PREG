#!/usr/bin/env Rscript
################################################################################
# Script: 04.2_fetal_effect_endoMR-PREG.R
# Project: endoMR-PREG
#
# Purpose:
#   Use trio-based GWAS (maternal / offspring / paternal effects) to estimate
#   MR effects of endometriosis genetic liability on pregnancy outcomes, and
#   compare:
#     - Maternal genetic effects
#     - Fetal (offspring) genetic effects
#     - Paternal genetic effects
#
# Input:
#   - results/Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv (exposure)
#   - data/OUTCOME_MR-PREG/trios_out_dat.txt (or .dat)             (outcomes)
#
# Output:
#   - results/trios_mr_results_maternal.csv
#   - results/trios_mr_results_fetal.csv
#   - results/trios_mr_results_paternal.csv
#   - results/trios_mr_results_comparison.csv
#   - results/plots/Figure_trios_maternal_fetal_paternal.png / .pdf
################################################################################

### 0) SETUP ###################################################################

required_pkgs <- c(
  "TwoSampleMR",
  "dplyr",
  "data.table",
  "here",
  "ggplot2",
  "readr",
  "tidyr"
)

safe_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

invisible(lapply(required_pkgs, safe_install))

project_dir <- here::here()
data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
plots_dir   <- file.path(results_dir, "plots")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)

message("=== TRIOS MR: maternal vs fetal vs paternal genetic effects ===")

################################################################################
# 1) LOAD EXPOSURE: ENDOMETRIOSIS (RAHMIOLU) INSTRUMENTS                        #
################################################################################

# Exposure file saved in Script 01 (clumped instruments)
exposure_file <- file.path(
  results_dir,
  "Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv"
)

stopifnot(file.exists(exposure_file))

exposure_dat <- data.table::fread(exposure_file, data.table = FALSE)

# Minimal sanity checks
stopifnot(all(c("SNP", "beta.exposure", "se.exposure") %in% colnames(exposure_dat)))
stopifnot("id.exposure" %in% colnames(exposure_dat))

message("Loaded exposure instruments: ", nrow(exposure_dat), " SNPs.")

################################################################################
# 2) LOAD TRIOS OUTCOME DATA                                                   #
################################################################################

# Try .txt, then fallback to no extension if needed
trios_path_txt <- file.path(data_dir, "OUTCOME_MR-PREG", "trios_out_dat.txt")
trios_path_raw <- file.path(data_dir, "OUTCOME_MR-PREG", "trios_out_dat")

if (file.exists(trios_path_txt)) {
  trios_path <- trios_path_txt
} else if (file.exists(trios_path_raw)) {
  trios_path <- trios_path_raw
} else {
  stop("Could not find trios_out_dat file in data/OUTCOME_MR-PREG/")
}

message("Reading trios file from: ", trios_path)

trios_raw <- data.table::fread(trios_path, data.table = FALSE)

# Expected columns (from your example)
expected_cols <- c(
  "SNP", "chr", "pos",
  "effect_allele", "other_allele",
  "beta_mat", "se_mat", "p_mat", "n_mat",
  "beta_mat_donuts", "se_mat_donuts", "p_mat_donuts",
  "eaf_mat", "n_studies_mat",
  "beta_off", "se_off", "p_off", "n_off",
  "beta_off_donuts", "se_off_donuts", "p_off_donuts",
  "beta_pat", "se_pat", "p_pat", "n_pat",
  "beta_pat_donuts", "se_pat_donuts", "p_pat_donuts",
  "Phenotype", "type", "analyses", "study"
)

missing_cols <- setdiff(expected_cols, colnames(trios_raw))
if (length(missing_cols)) {
  warning("The following expected columns are missing in trios_out_dat: ",
          paste(missing_cols, collapse = ", "))
}

message("Trios data: ", nrow(trios_raw), " rows, ",
        ncol(trios_raw), " columns, ",
        length(unique(trios_raw$Phenotype)), " distinct phenotypes.")

################################################################################
# 3) LIMIT TO THE 29 MAIN OUTCOMES                                             #
################################################################################

# Must match the 29 outcomes kept in your main script (04_main_analyses_endoMR-PREG.R)
outcome_labels_29 <- c(
  # Placental outcomes
  finngen_R12_O15_PLAC_PRAEVIA                      = "Placenta praevia",
  finngen_R12_O15_PLAC_DISORD                       = "Placental disorders",
  finngen_R12_O15_PLAC_PREMAT_SEPAR                 = "Premature placental separation",
  
  # Caesarean delivery
  el_cs                                            = "Elective caesarean section",
  em_cs                                            = "Emergency caesarean section",
  cs                                               = "Caesarean section",
  
  # Apgar scores
  lowapgar1                                        = "Low Apgar score at 1 min",
  lowapgar5                                        = "Low Apgar score at 5 min",
  
  # Labor and delivery complications
  rup_memb                                         = "Premature rupture of membranes",
  induction                                        = "Labour induction",
  
  # Gestational age and timing
  ga_all                                           = "Gestational age",
  pretb_all                                        = "Preterm birth (any)",
  vpretb_all                                       = "Very preterm birth",
  posttb_all                                       = "Post-term birth",
  
  # Birth weight outcomes
  hbw_all                                          = "High birthweight (>4000g)",
  lbw_all                                          = "Low birthweight (<2500g)",
  sga                                              = "Small for gestational age",
  
  # Maternal health
  depr_subsamp                                     = "Postpartum Depression",
  anaemia_preg_all                                 = "Pregnancy anemia",
  
  # Pregnancy complications
  gdm_subsamp                                      = "Gestational diabetes mellitus",
  hdp_subsamp                                      = "Hypertensive disorders of pregnancy",
  gh_subsamp                                       = "Gestational hypertension",
  pe_subsamp                                       = "Preeclampsia",
  
  # Neonatal outcomes
  nicu                                             = "NICU admission",
  sb_subsamp                                       = "Stillbirth",
  
  # Hemorrhage and bleeding (filtered versions in main analysis)
  Antepartum_bleeding_filtered                     = "Antepartum bleeding",
  Postpartum_hemorrhage_filtered                   = "Postpartum hemorrhage",
  Postpartum_hemorrhage_due_to_atony_filtered      = "PPH due to atony",
  Postpartum_hemorrhage_due_to_retained_placenta_filtered = "PPH due to retained placenta"
)

# Restrict trios data to phenotypes present in the 29-set
trios <- trios_raw %>%
  dplyr::filter(Phenotype %in% names(outcome_labels_29))

message("Trios rows after restricting to 29 outcomes: ", nrow(trios))
message("Phenotypes in trios after restriction: ",
        paste(sort(unique(trios$Phenotype)), collapse = ", "))

################################################################################
# 4) FORMAT TRIOS AS OUTCOME DATA (MATERNAL / FETAL / PATERNAL)                #
################################################################################

format_trios_component <- function(trios_df, component = c("maternal", "fetal", "paternal")) {
  component <- match.arg(component)
  
  if (component == "maternal") {
    beta_col <- "beta_mat"
    se_col   <- "se_mat"
    p_col    <- "p_mat"
    n_col    <- "n_mat"
    eaf_col  <- "eaf_mat"
    origin   <- "Maternal"
  } else if (component == "fetal") {
    beta_col <- "beta_off"
    se_col   <- "se_off"
    p_col    <- "p_off"
    n_col    <- "n_off"
    eaf_col  <- "eaf_mat"   # use maternal EAF as approximation
    origin   <- "Fetal"
  } else {
    beta_col <- "beta_pat"
    se_col   <- "se_pat"
    p_col    <- "p_pat"
    n_col    <- "n_pat"
    eaf_col  <- "eaf_mat"   # use maternal EAF as approximation
    origin   <- "Paternal"
  }
  
  dat_comp <- trios_df %>%
    dplyr::transmute(
      SNP            = SNP,
      chr            = chr,
      pos            = pos,
      effect_allele  = effect_allele,
      other_allele   = other_allele,
      beta           = !!rlang::sym(beta_col),
      se             = !!rlang::sym(se_col),
      pval           = !!rlang::sym(p_col),
      n              = !!rlang::sym(n_col),
      eaf            = !!rlang::sym(eaf_col),
      outcome        = Phenotype,
      id.outcome     = Phenotype,
      study          = study,
      component      = origin
    )
  
  # Drop rows with missing beta or se
  dat_comp <- dat_comp %>%
    dplyr::filter(
      !is.na(beta),
      !is.na(se),
      se > 0
    )
  
  # Format as outcome data for TwoSampleMR
  out_formatted <- TwoSampleMR::format_data(
    dat               = dat_comp,
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    samplesize_col    = "n",
    phenotype_col     = "outcome",
    id_col            = "id.outcome",
    chr_col           = "chr",
    pos_col           = "pos"
  )
  
  out_formatted$component <- origin
  out_formatted
}

outcome_mat <- format_trios_component(trios, "maternal")
outcome_fet <- format_trios_component(trios, "fetal")
outcome_pat <- format_trios_component(trios, "paternal")

message("Maternal outcome rows: ", nrow(outcome_mat))
message("Fetal outcome rows   : ", nrow(outcome_fet))
message("Paternal outcome rows: ", nrow(outcome_pat))

################################################################################
# 5) MR ANALYSIS FUNCTION                                                      #
################################################################################

run_trios_mr <- function(exposure_dat, outcome_dat, label_origin) {
  message("Harmonising and running MR for ", label_origin, " outcomes...")
  
  dat_harm <- TwoSampleMR::harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat  = outcome_dat,
    action       = 2
  )
  
  # MR (all default methods, but we will mainly use IVW)
  mr_res <- TwoSampleMR::mr(dat_harm)
  
  # Add human-readable labels
  mr_res <- mr_res %>%
    dplyr::mutate(
      outcome_full = dplyr::recode(id.outcome, !!!outcome_labels_29),
      origin       = label_origin
    )
  
  list(harmonised = dat_harm, results = mr_res)
}

################################################################################
# 6) RUN MR FOR MATERNAL, FETAL, PATERNAL                                     #
################################################################################

res_mat <- run_trios_mr(exposure_dat, outcome_mat, "Maternal")
res_fet <- run_trios_mr(exposure_dat, outcome_fet, "Fetal")
res_pat <- run_trios_mr(exposure_dat, outcome_pat, "Paternal")

mr_mat <- res_mat$results
mr_fet <- res_fet$results
mr_pat <- res_pat$results

# Keep IVW only for comparison
ivw_mat <- mr_mat %>% dplyr::filter(method == "Inverse variance weighted")
ivw_fet <- mr_fet %>% dplyr::filter(method == "Inverse variance weighted")
ivw_pat <- mr_pat %>% dplyr::filter(method == "Inverse variance weighted")

################################################################################
# 7) COMBINE AND FORMAT RESULTS                                               #
################################################################################

ivw_all <- dplyr::bind_rows(ivw_mat, ivw_fet, ivw_pat) %>%
  dplyr::mutate(
    OR   = exp(b),
    LCL  = exp(b - 1.96 * se),
    UCL  = exp(b + 1.96 * se),
    OR_CI = sprintf("%.2f (%.2f–%.2f)", OR, LCL, UCL),
    pval_fmt = formatC(pval, format = "e", digits = 2)
  )

# FDR correction across all (maternal+fetal+paternal) IVW tests
ivw_all$qval <- p.adjust(ivw_all$pval, method = "fdr")

# Export separate CSVs
readr::write_csv(ivw_mat, file.path(results_dir, "trios_mr_results_maternal.csv"))
readr::write_csv(ivw_fet, file.path(results_dir, "trios_mr_results_fetal.csv"))
readr::write_csv(ivw_pat, file.path(results_dir, "trios_mr_results_paternal.csv"))
readr::write_csv(ivw_all, file.path(results_dir, "trios_mr_results_comparison.csv"))

message("Saved MR results (maternal / fetal / paternal) in results/")

################################################################################
# 8) FOREST-PLOT COMPARISON: MATERNAL vs FETAL vs PATERNAL                   #
################################################################################

# Order outcomes by maternal OR magnitude (or any criterion you prefer)
ivw_all$order_outcome <- dplyr::dense_rank(dplyr::desc(abs(ivw_all$b)))

p_trios <- ggplot(ivw_all, aes(x = OR, y = reorder(outcome_full, order_outcome))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(xmin = LCL, xmax = UCL, colour = origin),
    height = 0.25,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    aes(colour = origin),
    size = 2.7,
    position = position_dodge(width = 0.6)
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.7, 1.0, 1.4, 2.0),
    labels = c("0.5", "0.7", "1.0", "1.4", "2.0")
  ) +
  labs(
    title    = "Maternal vs fetal vs paternal genetic effects",
    subtitle = "Endometriosis genetic liability (trios-based MR)",
    x        = "Odds ratio (95% CI, log scale)",
    y        = "Pregnancy outcome"
  ) +
  scale_colour_manual(
    values = c(
      "Maternal" = "#1f78b4",
      "Fetal"    = "#33a02c",
      "Paternal" = "#e31a1c"
    ),
    name = "Genetic effect"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(plots_dir, "Figure_trios_maternal_fetal_paternal.png"),
  p_trios, width = 9, height = 7, dpi = 300
)
ggsave(
  file.path(plots_dir, "Figure_trios_maternal_fetal_paternal.pdf"),
  p_trios, width = 9, height = 7
)

message("Saved comparison plot to results/plots/Figure_trios_maternal_fetal_paternal.*")

################################################################################
# 9) SUMMARY MESSAGE                                                          #
################################################################################

message("=== TRIOS MR ANALYSIS COMPLETE ===")
message("  - CSVs: trios_mr_results_*.csv in results/")
message("  - Figure: Figure_trios_maternal_fetal_paternal.(png/pdf) in results/plots/")
