### 5. SENSITIVITY ANALYSES ######################################################

# -------------------------------------------------------------------------------
# Helpers (safe fallbacks if not already defined elsewhere)
# -------------------------------------------------------------------------------

if (!exists("export_csv")) {
  export_csv <- function(df, stem, dir_out = results_dir) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
      message("Nothing to export for ", stem)
      return(invisible(NULL))
    }
    if (missing(dir_out) || is.null(dir_out)) dir_out <- "."
    if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    fn <- file.path(dir_out, paste0(stem, "_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv"))
    readr::write_csv(df, fn)
    message("Written: ", fn)
    invisible(fn)
  }
}

if (!exists("labels_df")) {
  labels_df <- data.frame(outcome = character(), outcome_full = character(), stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------------------
# Standard TwoSampleMR sensitivity outputs from harmonised dataset
# -------------------------------------------------------------------------------

het_res     <- mr_heterogeneity(dat)      # 5.1 Cochran’s Q
plt_res     <- mr_pleiotropy_test(dat)    # 5.2 Egger intercept
single_res  <- mr_singlesnp(dat)          # 5.3 Single-SNP effects
loo_snp_res <- mr_leaveoneout(dat)        # 5.4 Leave-one-out by SNP


###############################################################################
# 5.5 LEAVE-ONE-COHORT-OUT MENDELIAN RANDOMIZATION (MR) ANALYSIS
###############################################################################

message("=== LEAVE-ONE-COHORT-OUT MR ANALYSIS ===")

# -------------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------------

message("Loading packages...")
required_pkgs <- c(
  "TwoSampleMR", "dplyr", "ggplot2", "here", "readr",
  "data.table", "gridExtra", "scales", "stringr", "tidyr", "purrr"
)
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(lapply(required_pkgs, function(p) library(p, character.only = TRUE)))
set.seed(1)

PLOT_ONLY_SIGNIFICANT <- TRUE
MAX_IN_COMBINED_PANEL <- 6

results_dir <- here::here("results")
plots_dir   <- here::here("plots")
dir.create(results_dir, FALSE, TRUE)
dir.create(plots_dir, FALSE, TRUE)

outcome_file <- here::here("data", "OUTCOME_MR-PREG", "stu_out_dat.txt")
harm_file    <- here::here("results", "harmonised_rahmioglu_bpo.csv")

# -------------------------------------------------------------------------------
# Load and validate data
# -------------------------------------------------------------------------------

message("Loading outcome data...")
if (!file.exists(outcome_file)) stop("Outcome file not found: ", outcome_file)
if (!file.exists(harm_file))    stop("Harmonized exposure file not found: ", harm_file)

mr_data <- readr::read_table(outcome_file, col_types = cols()) |> as.data.frame()
required_cols_out <- c("SNP", "beta", "se", "pval", "samplesize", "study",
                       "Phenotype", "effect_allele", "other_allele", "eaf",
                       "ncase", "ncontrol")
miss_out <- setdiff(required_cols_out, names(mr_data))
if (length(miss_out)) stop("Missing columns in outcome data: ", paste(miss_out, collapse = ", "))

harm_data <- data.table::fread(harm_file) |> as.data.frame()
required_cols_exp <- c("SNP", "beta.exposure", "se.exposure", "pval.exposure",
                       "effect_allele.exposure", "other_allele.exposure",
                       "eaf.exposure", "samplesize.exposure",
                       "id.exposure", "exposure")
miss_exp <- setdiff(required_cols_exp, names(harm_data))
if (length(miss_exp)) stop("Missing columns in exposure data: ", paste(miss_exp, collapse = ", "))

message("Outcome rows: ", nrow(mr_data), " | Exposure rows: ", nrow(harm_data))

# -------------------------------------------------------------------------------
# Cohorts and labels
# -------------------------------------------------------------------------------

cohort_summary <- mr_data |>
  group_by(study) |>
  summarise(
    n_snps = n_distinct(SNP),
    n_phenotypes = n_distinct(Phenotype),
    mean_samplesize = mean(samplesize, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_snps))

message("Available cohorts: ", paste(cohort_summary$study, collapse = ", "))
if (nrow(cohort_summary) < 2) stop("Need ≥2 cohorts for leave-one-out analysis")

outcome_labels <- c(
  pretb_all = "Preterm birth (any)",
  el_cs = "Elective caesarean section",
  rup_memb = "Premature rupture of membranes",
  pretb_subsamp = "Preterm birth (spontaneous)",
  apgar1 = "Apgar score at 1 min",
  vpretb_all = "Very preterm birth (any)",
  em_cs = "Emergency caesarean section",
  lowapgar1 = "Low Apgar score at 1 min",
  pe_subsamp = "Preeclampsia",
  lbw_all = "Low birthweight (<2500 g)",
  ga_all = "Gestational age",
  gh_subsamp = "Gestational hypertension",
  ga_subsamp = "Gestational age (subset)",
  lga = "Large for gestational age",
  zbw_all = "Z-score birthweight",
  nvp_sev_subsamp = "Severe nausea/vomiting (subset)",
  nvp_sev_all = "Severe nausea/vomiting",
  hbw_all = "High birthweight (>4000 g)",
  anaemia_preg_all = "Pregnancy anaemia",
  finngen_R12_N14_FEMALEINFERT = "Female infertility",
  finngen_R12_O15_PLAC_DISORD = "Placental disorders",
  finngen_R12_O15_PLAC_PRAEVIA = "Placenta praevia",
  finngen_R12_O15_PLAC_PREMAT_SEPAR = "Premature placental separation",
  finngen_R12_O15_PREG_ECTOP = "Ectopic pregnancy"
)

outcome_types <- c(
  pretb_all = "binary", el_cs = "binary", rup_memb = "binary",
  pretb_subsamp = "binary", apgar1 = "continuous", vpretb_all = "binary",
  em_cs = "binary", lowapgar1 = "binary", pe_subsamp = "binary",
  lbw_all = "binary", ga_all = "continuous", gh_subsamp = "binary",
  ga_subsamp = "continuous", lga = "binary", zbw_all = "continuous",
  nvp_sev_subsamp = "binary", nvp_sev_all = "binary", hbw_all = "binary",
  anaemia_preg_all = "binary", finngen_R12_N14_FEMALEINFERT = "binary",
  finngen_R12_O15_PLAC_DISORD = "binary", finngen_R12_O15_PLAC_PRAEVIA = "binary",
  finngen_R12_O15_PLAC_PREMAT_SEPAR = "binary", finngen_R12_O15_PREG_ECTOP = "binary"
)

mr_data <- mr_data |>
  mutate(
    id.outcome = Phenotype,
    outcome_label = dplyr::coalesce(outcome_labels[Phenotype], Phenotype),
    outcome_type = dplyr::coalesce(outcome_types[Phenotype], "binary")
  )

# -------------------------------------------------------------------------------
# Exposure preparation
# -------------------------------------------------------------------------------

message("Preparing exposure dataset...")

exposure_df_raw <- harm_data |>
  select(SNP, beta.exposure, se.exposure, pval.exposure,
         effect_allele.exposure, other_allele.exposure, eaf.exposure,
         samplesize.exposure, id.exposure, exposure) |>
  distinct(SNP, .keep_all = TRUE)

exposure_for_format <- exposure_df_raw |>
  transmute(
    SNP,
    beta = beta.exposure, se = se.exposure, pval = pval.exposure,
    effect_allele = effect_allele.exposure, other_allele = other_allele.exposure,
    eaf = eaf.exposure, samplesize = samplesize.exposure,
    ncase = NA_real_, ncontrol = NA_real_,
    id = id.exposure, exposure = exposure, phenotype = exposure
  )

exposure_formatted <- TwoSampleMR::format_data(
  exposure_for_format, type = "exposure", phenotype_col = "phenotype"
)

message("Exposure SNPs: ", dplyr::n_distinct(exposure_formatted$SNP))

# -------------------------------------------------------------------------------
# Helpers: meta-analysis and harmonisation
# -------------------------------------------------------------------------------

meta_outcome_after_drop <- function(dat_one_outcome) {
  dat_one_outcome |>
    filter(!is.na(SNP), is.finite(beta.outcome), is.finite(se.outcome)) |>
    group_by(SNP) |>
    summarise(
      w = sum(1 / (se.outcome^2), na.rm = TRUE),
      beta = ifelse(w > 0, sum(beta.outcome / (se.outcome^2), na.rm = TRUE) / w, NA_real_),
      se = ifelse(w > 0, sqrt(1 / w), NA_real_),
      pval = ifelse(w > 0, 2 * pnorm(-abs(beta / se)), NA_real_),
      effect_allele = { tmp <- na.omit(effect_allele.outcome); if (length(tmp)) tmp[1] else NA_character_ },
      other_allele  = { tmp <- na.omit(other_allele.outcome); if (length(tmp)) tmp[1] else NA_character_ },
      ss_sum = sum(samplesize.outcome, na.rm = TRUE),
      eaf = ifelse(ss_sum > 0 & any(!is.na(eaf.outcome)),
                   sum(eaf.outcome * samplesize.outcome, na.rm = TRUE) / ss_sum,
                   mean(eaf.outcome, na.rm = TRUE)),
      samplesize = ss_sum,
      ncase = sum(ncase.outcome, na.rm = TRUE),
      ncontrol = sum(ncontrol.outcome, na.rm = TRUE),
      .groups = "drop"
    ) |>
    select(-w, -ss_sum) |>
    filter(is.finite(beta), is.finite(se), !is.na(effect_allele), !is.na(other_allele))
}

build_outcome_formatted <- function(meta_df, outcome_id, outcome_label) {
  if (is.null(meta_df) || nrow(meta_df) == 0) return(NULL)
  tmp <- meta_df |>
    transmute(
      SNP, beta, se, pval, effect_allele, other_allele, eaf,
      samplesize, ncase, ncontrol,
      id = outcome_id, outcome = outcome_label, phenotype = outcome_label
    ) |>
    filter(!is.na(SNP), is.finite(beta), is.finite(se),
           !is.na(effect_allele), !is.na(other_allele))
  if (nrow(tmp) == 0) return(NULL)
  TwoSampleMR::format_data(tmp, type = "outcome", phenotype_col = "phenotype")
}

# -------------------------------------------------------------------------------
# Main leave-one-cohort-out MR
# -------------------------------------------------------------------------------

message("Preparing outcome data frame...")

mr_outcome_raw <- mr_data |>
  transmute(
    SNP,
    beta.outcome = beta, se.outcome = se, pval.outcome = pval,
    effect_allele.outcome = effect_allele, other_allele.outcome = other_allele,
    eaf.outcome = eaf, samplesize.outcome = samplesize,
    ncase.outcome = ncase, ncontrol.outcome = ncontrol,
    id.outcome, outcome = outcome_label, outcome_type, study
  )

outcomes_list <- split(mr_outcome_raw, mr_outcome_raw$id.outcome)

run_loo_for_outcome <- function(one_outcome_df) {
  outcome_id   <- unique(one_outcome_df$id.outcome)[1]
  outcome_name <- unique(one_outcome_df$outcome)[1]
  cohorts      <- sort(unique(one_outcome_df$study))
  message(">> Outcome: ", outcome_id, " | ", outcome_name, " | Cohorts: ", length(cohorts))
  
  # Main analysis
  main_meta <- meta_outcome_after_drop(one_outcome_df)
  main_outcome_fmt <- build_outcome_formatted(main_meta, outcome_id, outcome_name)
  if (is.null(main_outcome_fmt) || nrow(main_outcome_fmt) < 3) return(NULL)
  
  main_harm <- tryCatch(TwoSampleMR::harmonise_data(exposure_formatted, main_outcome_fmt), error = function(e) NULL)
  if (is.null(main_harm) || dplyr::n_distinct(main_harm$SNP) < 3) return(NULL)
  
  main_ivw <- TwoSampleMR::mr(main_harm, method_list = "mr_ivw")
  if (!nrow(main_ivw)) return(NULL)
  
  main_ivw$analysis_type <- "Main analysis"
  main_ivw$nsnp <- dplyr::n_distinct(main_harm$SNP)
  main_ivw$left_out_cohort <- "None"
  main_ivw$n_snps_remaining <- main_ivw$nsnp
  main_ivw$n_cohorts_remaining <- length(cohorts)
  main_ivw$outcome <- outcome_name
  main_ivw$id.outcome <- outcome_id
  
  # Leave-one-cohort-out
  loo_list <- lapply(cohorts, function(left_out) {
    sub_df <- dplyr::filter(one_outcome_df, study != left_out)
    meta <- meta_outcome_after_drop(sub_df)
    out_fmt <- build_outcome_formatted(meta, outcome_id, outcome_name)
    if (is.null(out_fmt) || nrow(out_fmt) < 3) return(NULL)
    harm <- tryCatch(TwoSampleMR::harmonise_data(exposure_formatted, out_fmt), error = function(e) NULL)
    if (is.null(harm) || dplyr::n_distinct(harm$SNP) < 3) return(NULL)
    res <- TwoSampleMR::mr(harm, method_list = "mr_ivw")
    if (!nrow(res)) return(NULL)
    res$analysis_type <- paste0("Without ", left_out)
    res$left_out_cohort <- left_out
    res$n_snps_remaining <- dplyr::n_distinct(harm$SNP)
    res$n_cohorts_remaining <- dplyr::n_distinct(sub_df$study)
    res$outcome <- outcome_name
    res$id.outcome <- outcome_id
    res
  })
  
  dplyr::bind_rows(main_ivw, dplyr::bind_rows(loo_list))
}

message("Running main + leave-one-cohort-out analyses...")

results_list <- lapply(outcomes_list, run_loo_for_outcome)
all_results_raw <- dplyr::bind_rows(results_list)
if (!nrow(all_results_raw)) stop("No MR results produced")

# -------------------------------------------------------------------------------
# Robustness summary and export
# -------------------------------------------------------------------------------

all_results <- all_results_raw |>
  mutate(
    outcome = dplyr::coalesce(outcome, outcome_labels[id.outcome], id.outcome),
    outcome_type = dplyr::coalesce(outcome_types[id.outcome], "binary"),
    effect_label = ifelse(outcome_type == "binary", "OR", "Beta"),
    est_display = ifelse(outcome_type == "binary", exp(b), b),
    ci_low = ifelse(outcome_type == "binary", exp(b - 1.96 * se), b - 1.96 * se),
    ci_high = ifelse(outcome_type == "binary", exp(b + 1.96 * se), b + 1.96 * se),
    est_ci = sprintf("%.3f (%.3f–%.3f)", est_display, ci_low, ci_high),
    p_display = case_when(
      is.na(pval) ~ NA_character_,
      pval < 0.001 ~ "<0.001",
      pval < 0.01 ~ sprintf("%.3f", pval),
      TRUE ~ sprintf("%.3f", pval)
    )
  ) |>
  select(id.outcome, outcome, method, b, se, pval, nsnp, analysis_type,
         left_out_cohort, n_snps_remaining, n_cohorts_remaining,
         outcome_type, effect_label, est_display, ci_low, ci_high, est_ci, p_display)

readr::write_csv(all_results, file.path(results_dir, "leave_one_cohort_results.csv"))
readr::write_csv(cohort_summary, file.path(results_dir, "cohort_composition.csv"))

message("Saved results to: ", results_dir)
