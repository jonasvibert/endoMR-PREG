#!/usr/bin/env Rscript
################################################################################
# Script: 04.2_fetal_effect_endoMR-PREG.R
# Project: endoMR-PREG
#
# Purpose:
#   Use trio-based GWAS (maternal / offspring / paternal effects) to estimate
#   MR effects of endometriosis genetic liability on pregnancy outcomes.
#
#   Uses the MUTUALLY ADJUSTED (DONUTS) estimates from the trio-based GWAS,
#   following the Lawlor group standard (Warrington et al. 2019, Nat Genet):
#
#     (1) Maternal effect adjusted for fetal genotype  → beta_mat_donuts
#         Model: Outcome = β_m·G_m + β_f·G_f
#         Interpretation: pure maternal genetic effect, not confounded by
#         maternally inherited alleles acting via fetal biology.
#
#     (2) Fetal effect adjusted for maternal genotype  → beta_off_donuts
#         Sensitivity analysis: is the signal primarily fetal?
#
#     (3) Paternal effect (conditioned)               → beta_pat_donuts
#         Comparison arm: negative control for fetal effects;
#         paternal genotype can influence offspring via transmitted alleles
#         but does not share the maternal intrauterine environment.
#
# Input:
#   - results/Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv  (exposure)
#   - data/OUTCOME_MR-PREG/trios_out_dat.txt                     (outcomes)
#
# Output:
#   - results/trios_adj_mr_results_maternal.csv
#   - results/trios_adj_mr_results_fetal.csv
#   - results/trios_adj_mr_results_paternal.csv
#   - results/trios_adj_mr_results_comparison.csv
#   - results/trios_adj_mr_results_by_outcome_long.csv
#   - results/plots/Figure_trios_adj_maternal_fetal_paternal.png / .pdf
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

set.seed(42)

message("=== TRIOS MR (ADJUSTED): maternal vs fetal vs paternal — conditioned estimates ===")

################################################################################
# 1) LOAD EXPOSURE: ENDOMETRIOSIS (RAHMIOGLU) INSTRUMENTS                      #
################################################################################

exposure_file <- file.path(
  results_dir,
  "Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv"
)

stopifnot(file.exists(exposure_file))

exposure_dat <- data.table::fread(exposure_file, data.table = FALSE)

stopifnot(all(c("SNP", "beta.exposure", "se.exposure") %in% colnames(exposure_dat)))
stopifnot("id.exposure" %in% colnames(exposure_dat))

message("Loaded exposure instruments: ", nrow(exposure_dat), " SNPs.")

################################################################################
# 2) LOAD TRIOS OUTCOME DATA                                                   #
################################################################################

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

# Columns required for the ADJUSTED analysis
required_donuts_cols <- c(
  "SNP", "chr", "pos",
  "effect_allele", "other_allele", "eaf_mat",
  "beta_mat_donuts", "se_mat_donuts", "p_mat_donuts", "n_mat",
  "beta_off_donuts", "se_off_donuts", "p_off_donuts", "n_off",
  "beta_pat_donuts", "se_pat_donuts", "p_pat_donuts", "n_pat",
  "Phenotype", "study"
)

missing_cols <- setdiff(required_donuts_cols, colnames(trios_raw))
if (length(missing_cols)) {
  stop(
    "The following required DONUTS columns are missing in trios_out_dat:\n  ",
    paste(missing_cols, collapse = ", "),
    "\nThis script requires the conditioned (donuts) estimates."
  )
}

message(
  "Trios data: ", nrow(trios_raw), " rows, ",
  ncol(trios_raw), " columns, ",
  length(unique(trios_raw$Phenotype)), " distinct phenotypes."
)

# Diagnostic: NA counts in donuts columns before any filtering
message("\n--- NA counts in conditioned (donuts) columns (full dataset) ---")
for (col in c("beta_mat_donuts", "beta_off_donuts", "beta_pat_donuts")) {
  n_na    <- sum(is.na(trios_raw[[col]]))
  n_total <- nrow(trios_raw)
  message(sprintf("  %-22s : %d NA / %d rows (%.1f%%)",
                  col, n_na, n_total, 100 * n_na / n_total))
}

################################################################################
# 3) LIMIT TO THE 30 MAIN OUTCOMES                                             #
################################################################################

# Must match the 30 outcomes kept in 04_main_analyses_endoMR-PREG.R

vars_keep_trios <- c(
  # Placenta & bleeding (7)
  "Antepartum_bleeding",
  "Postpartum_hemorrhage",
  "Postpartum_hemorrhage_due_to_atony",
  "Postpartum_hemorrhage_due_to_retained_placenta",
  "finngen_R12_O15_PLAC_PRAEVIA",
  "finngen_R12_O15_PLAC_DISORD",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",

  # Membranes (1)
  "rup_memb",

  # Preterm birth (2)
  "pretb_all",
  "vpretb_all",

  # Growth / GA / weight (6)
  "ga_all",
  "sga",
  "lbw_all",
  "hbw_all",
  "lga",
  "zbw_all",

  # Neonatal / Apgar / perinatal death (4)
  "lowapgar1",
  "lowapgar5",
  "nicu",
  "sb_subsamp",

  # Maternal complications (6)
  "anaemia_preg_all",
  "gdm_subsamp",
  "gh_subsamp",
  "hdp_subsamp",
  "pe_subsamp",
  "depr_subsamp",

  # Other obstetric (4)
  "induction",
  "posttb_all",
  "el_cs",
  "em_cs"
)

outcome_labels_30 <- c(
  # Bleeding (4)
  Antepartum_bleeding                            = "Antepartum bleeding",
  Postpartum_hemorrhage                          = "Postpartum hemorrhage (any)",
  Postpartum_hemorrhage_due_to_atony             = "PPH due to atony",
  Postpartum_hemorrhage_due_to_retained_placenta = "PPH due to retained placenta",

  # Placenta (3)
  finngen_R12_O15_PLAC_PRAEVIA                   = "Placenta praevia",
  finngen_R12_O15_PLAC_DISORD                    = "Placental disorders",
  finngen_R12_O15_PLAC_PREMAT_SEPAR              = "Premature placental separation",

  # Membranes (1)
  rup_memb                                       = "Premature rupture of membranes",

  # Birth timing (4)
  pretb_all                                      = "Preterm birth <37 weeks (any)",
  vpretb_all                                     = "Very preterm birth <34 weeks",
  posttb_all                                     = "Post-term birth",
  ga_all                                         = "Gestational age",

  # Fetal growth (5)
  sga                                            = "Small for gestational age",
  lbw_all                                        = "Low birthweight <2500g",
  hbw_all                                        = "High birthweight >4000g",
  lga                                            = "Large for gestational age",
  zbw_all                                        = "Z-score birthweight",

  # Neonatal adaptation / perinatal death (4)
  lowapgar1                                      = "Low Apgar score at 1 min",
  lowapgar5                                      = "Low Apgar score at 5 min",
  nicu                                           = "NICU admission",
  sb_subsamp                                     = "Stillbirth",

  # Maternal complications (6)
  anaemia_preg_all                               = "Pregnancy anemia",
  gdm_subsamp                                    = "Gestational diabetes",
  gh_subsamp                                     = "Gestational hypertension",
  hdp_subsamp                                    = "Hypertensive disorders of pregnancy",
  pe_subsamp                                     = "Preeclampsia",
  depr_subsamp                                   = "Postpartum depression",

  # Caesarean section (2)
  el_cs                                          = "Elective caesarean section",
  em_cs                                          = "Emergency caesarean section",

  # Other obstetric timing (1)
  induction                                      = "Labour induction"
)

trios <- trios_raw %>%
  dplyr::filter(Phenotype %in% vars_keep_trios)

message("\nTrios rows after restricting to 30 outcomes: ", nrow(trios))
message(
  "Phenotypes present after restriction: ",
  paste(sort(unique(trios$Phenotype)), collapse = ", ")
)

# Diagnostic: NA counts after restriction to 30 outcomes
message("\n--- NA counts in conditioned (donuts) columns (30 outcomes) ---")
for (col in c("beta_mat_donuts", "beta_off_donuts", "beta_pat_donuts")) {
  n_na    <- sum(is.na(trios[[col]]))
  n_total <- nrow(trios)
  message(sprintf("  %-22s : %d NA / %d rows (%.1f%%)",
                  col, n_na, n_total, 100 * n_na / n_total))
}

################################################################################
# 4) FORMAT TRIOS AS OUTCOME DATA — USING CONDITIONED (DONUTS) ESTIMATES       #
#                                                                              #
#   beta_mat_donuts : maternal effect conditioned on offspring genotype        #
#                     → removes fetal genetic confounding of maternal estimate #
#   beta_off_donuts : offspring (fetal) effect conditioned on maternal genotype#
#                     → isolates fetal genetic signal                          #
#   beta_pat_donuts : paternal effect conditioned on offspring genotype        #
#                     → negative control; no shared intrauterine environment   #
################################################################################

format_trios_component <- function(trios_df, component = c("maternal_unadj", "maternal", "fetal", "paternal")) {
  component <- match.arg(component)

  if (component == "maternal_unadj") {
    beta_col <- "beta_mat"          # maternal NOT adjusted for fetal genotype
    se_col   <- "se_mat"
    p_col    <- "p_mat"
    n_col    <- "n_mat"
    eaf_col  <- "eaf_mat"
    origin   <- "Maternal (unadj.)"
  } else if (component == "maternal") {
    beta_col <- "beta_mat_donuts"   # maternal adjusted for fetal genotype
    se_col   <- "se_mat_donuts"
    p_col    <- "p_mat_donuts"
    n_col    <- "n_mat"
    eaf_col  <- "eaf_mat"
    origin   <- "Maternal (adj. fetal)"
  } else if (component == "fetal") {
    beta_col <- "beta_off_donuts"   # fetal adjusted for maternal genotype
    se_col   <- "se_off_donuts"
    p_col    <- "p_off_donuts"
    n_col    <- "n_off"
    eaf_col  <- "eaf_mat"           # trios file provides one EAF (maternal); used for all components
    origin   <- "Fetal (adj.)"
  } else {
    beta_col <- "beta_pat_donuts"   # paternal conditioned
    se_col   <- "se_pat_donuts"
    p_col    <- "p_pat_donuts"
    n_col    <- "n_pat"
    eaf_col  <- "eaf_mat"           # trios file provides one EAF (maternal); used for all components
    origin   <- "Paternal (adj.)"
  }

  n_before <- nrow(trios_df)

  dat_comp <- trios_df %>%
    dplyr::transmute(
      SNP           = SNP,
      chr           = chr,
      pos           = pos,
      effect_allele = effect_allele,
      other_allele  = other_allele,
      beta          = !!rlang::sym(beta_col),
      se            = !!rlang::sym(se_col),
      pval          = !!rlang::sym(p_col),
      n             = !!rlang::sym(n_col),
      eaf           = !!rlang::sym(eaf_col),
      outcome       = Phenotype,
      id.outcome    = Phenotype,
      study         = study,
      component     = origin
    ) %>%
    dplyr::filter(
      !is.na(beta),
      !is.na(se),
      se > 0
    )

  n_after <- nrow(dat_comp)
  message(sprintf(
    "  [%s] %d rows before NA filter → %d retained (dropped %d)",
    origin, n_before, n_after, n_before - n_after
  ))

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

message("\n--- Formatting outcome data ---")
outcome_mat_unadj <- format_trios_component(trios, "maternal_unadj")
outcome_mat       <- format_trios_component(trios, "maternal")
outcome_fet       <- format_trios_component(trios, "fetal")
outcome_pat       <- format_trios_component(trios, "paternal")

message("\nFormatted outcome rows:")
message("  Maternal (unadj.)    : ", nrow(outcome_mat_unadj))
message("  Maternal (adj. fetal): ", nrow(outcome_mat))
message("  Fetal (adj.)         : ", nrow(outcome_fet))
message("  Paternal (adj.)      : ", nrow(outcome_pat))

################################################################################
# 5) MR ANALYSIS FUNCTION                                                      #
################################################################################

run_trios_mr <- function(exposure_dat, outcome_dat, label_origin) {
  message("\nHarmonising and running MR for: ", label_origin)

  dat_harm <- TwoSampleMR::harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat  = outcome_dat,
    action       = 2
  )
  dat_harm <- dat_harm[dat_harm$mr_keep == TRUE, ]

  mr_res <- TwoSampleMR::mr(dat_harm)

  mr_res <- mr_res %>%
    dplyr::mutate(
      outcome_full = dplyr::recode(
        id.outcome,
        !!!outcome_labels_30,
        .default = id.outcome
      ),
      origin = label_origin
    )

  list(harmonised = dat_harm, results = mr_res)
}

################################################################################
# 6) RUN MR FOR ALL FOUR COMPONENTS                                            #
#                                                                              #
#   Main plot (Figure 3): Maternal (unadj.) / Maternal (adj. fetal) /         #
#                         Paternal (adj.)                                      #
#   Supplementary Table S4: all four components including Fetal (adj.)        #
################################################################################

message("\n=== Running MR analyses ===")

res_mat_unadj <- run_trios_mr(exposure_dat, outcome_mat_unadj, "Maternal (unadj.)")
res_mat       <- run_trios_mr(exposure_dat, outcome_mat,       "Maternal (adj. fetal)")
res_fet       <- run_trios_mr(exposure_dat, outcome_fet,       "Fetal (adj.)")
res_pat       <- run_trios_mr(exposure_dat, outcome_pat,       "Paternal (adj.)")

mr_mat_unadj <- res_mat_unadj$results
mr_mat       <- res_mat$results
mr_fet       <- res_fet$results
mr_pat       <- res_pat$results

ivw_mat_unadj <- mr_mat_unadj %>% dplyr::filter(method == "Inverse variance weighted")
ivw_mat       <- mr_mat       %>% dplyr::filter(method == "Inverse variance weighted")
ivw_fet       <- mr_fet       %>% dplyr::filter(method == "Inverse variance weighted")
ivw_pat       <- mr_pat       %>% dplyr::filter(method == "Inverse variance weighted")

message("\nIVW results per component:")
message("  Maternal (unadj.)    : ", nrow(ivw_mat_unadj))
message("  Maternal (adj. fetal): ", nrow(ivw_mat))
message("  Fetal (adj.)         : ", nrow(ivw_fet))
message("  Paternal (adj.)      : ", nrow(ivw_pat))

################################################################################
# 7) COMBINE AND FORMAT RESULTS                                                #
#                                                                              #
#   Binary outcomes  : OR = exp(b), LCL/UCL exponentiated                     #
#   Continuous (zbw_all, ga_all) : OR column = β (NOT exponentiated)          #
#   FDR correction applied per component (30 tests each), not across all 90   #
################################################################################

# Continuous outcomes: report β (do NOT exponentiate)
continuous_outcomes <- c("zbw_all", "ga_all")

ivw_all <- dplyr::bind_rows(ivw_mat_unadj, ivw_mat, ivw_fet, ivw_pat) %>%
  dplyr::mutate(
    scale    = if_else(id.outcome %in% continuous_outcomes, "continuous", "binary"),
    # OR/LCL/UCL only meaningful for binary outcomes
    OR       = if_else(scale == "binary", exp(b),             b),
    LCL      = if_else(scale == "binary", exp(b - 1.96 * se), b - 1.96 * se),
    UCL      = if_else(scale == "binary", exp(b + 1.96 * se), b + 1.96 * se),
    # Label column: "OR (95% CI)" for binary, "β (95% CI)" for continuous
    est_CI   = if_else(
      scale == "binary",
      sprintf("%.2f (%.2f\u2013%.2f)", OR, LCL, UCL),
      sprintf("%.3f (%.3f\u2013%.3f)", OR, LCL, UCL)
    ),
    pval_fmt = formatC(pval, format = "e", digits = 2)
  )

# FDR applied per component (30 tests each), not across all 90 combined
ivw_all <- ivw_all %>%
  dplyr::group_by(origin) %>%
  dplyr::mutate(qval = p.adjust(pval, method = "fdr")) %>%
  dplyr::ungroup()

readr::write_csv(ivw_mat_unadj, file.path(results_dir, "trios_adj_mr_results_maternal_unadj.csv"))
readr::write_csv(ivw_mat,       file.path(results_dir, "trios_adj_mr_results_maternal_adj.csv"))
readr::write_csv(ivw_fet,       file.path(results_dir, "trios_adj_mr_results_fetal.csv"))
readr::write_csv(ivw_pat,       file.path(results_dir, "trios_adj_mr_results_paternal.csv"))
readr::write_csv(ivw_all,       file.path(results_dir, "trios_adj_mr_results_comparison.csv"))

message("\nSaved IVW results (conditioned) to results/trios_adj_mr_results_*.csv")

# Long table: one row per outcome × component, rounded to 3 decimals ----------

ivw_by_outcome <- ivw_all %>%
  dplyr::select(
    outcome_id   = id.outcome,
    outcome_full,
    scale,
    origin,
    nsnp,
    b, se,
    OR, LCL, UCL,
    pval, qval
  ) %>%
  dplyr::mutate(
    origin = factor(origin, levels = c("Maternal (unadj.)", "Maternal (adj. fetal)", "Fetal (adj.)", "Paternal (adj.)")),
    b      = round(b,   3),
    se     = round(se,  3),
    OR     = round(OR,  3),   # = β for continuous outcomes
    LCL    = round(LCL, 3),
    UCL    = round(UCL, 3),
    pval   = signif(pval, 3),
    qval   = signif(qval, 3),
    est_CI = if_else(
      scale == "binary",
      sprintf("%.3f (%.3f\u2013%.3f)", OR, LCL, UCL),
      sprintf("%.3f (%.3f\u2013%.3f)", OR, LCL, UCL)
    )
  ) %>%
  dplyr::arrange(outcome_full, origin)

readr::write_csv(
  ivw_by_outcome,
  file.path(results_dir, "trios_adj_mr_results_by_outcome_long.csv")
)

message("Saved long table to results/trios_adj_mr_results_by_outcome_long.csv")

################################################################################
# 8) FOREST PLOT: MATERNAL (adj.) vs FETAL (adj.) vs PATERNAL (adj.)          #
################################################################################

ivw_all$order_outcome <- dplyr::dense_rank(dplyr::desc(abs(ivw_all$b)))

p_trios <- ggplot(
  ivw_all,
  aes(x = OR, y = reorder(outcome_full, order_outcome))
) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(xmin = LCL, xmax = UCL, colour = origin),
    height   = 0.25,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    aes(colour = origin),
    size     = 2.7,
    position = position_dodge(width = 0.6)
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.7, 1.0, 1.4, 2.0),
    labels = c("0.5", "0.7", "1.0", "1.4", "2.0")
  ) +
  scale_colour_manual(
    values = c(
      "Maternal (adj.)" = "#1f78b4",
      "Fetal (adj.)"    = "#33a02c",
      "Paternal (adj.)" = "#e31a1c"
    ),
    name = "Genetic effect"
  ) +
  labs(
    title    = "Maternal vs fetal vs paternal genetic effects (conditioned estimates)",
    subtitle = "Endometriosis genetic liability \u2014 trio-based MR (beta_*_donuts)",
    x        = "Odds ratio (95% CI, log scale)",
    y        = "Pregnancy outcome"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(plots_dir, "Figure_trios_adj_maternal_fetal_paternal.png"),
  p_trios, width = 9, height = 7, dpi = 300
)
ggsave(
  file.path(plots_dir, "Figure_trios_adj_maternal_fetal_paternal.pdf"),
  p_trios, width = 9, height = 7
)

message("Saved forest plot to results/plots/Figure_trios_adj_maternal_fetal_paternal.*")

################################################################################
# 9) SUMMARY                                                                   #
################################################################################

message("\n=== TRIOS MR (ADJUSTED) COMPLETE ===")
message("  Estimates used:")
message("    Maternal : beta_mat_donuts (adjusted for fetal genotype)")
message("    Fetal    : beta_off_donuts (adjusted for maternal genotype)")
message("    Paternal : beta_pat_donuts (conditioned)")
message("  CSVs    : results/trios_adj_mr_results_*.csv")
message("  Figure  : results/plots/Figure_trios_adj_maternal_fetal_paternal.*")
