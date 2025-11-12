#!/usr/bin/env Rscript
# 5_Tables_MR_analysis_endoMR-PREG.R
###############################################################################
# Generate summary tables for endometriosis → pregnancy outcomes MR analysis
###############################################################################

# =========================
# Table 1: Sample sizes and SNP counts for pregnancy outcomes
# =========================
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(readr); library(here)
})

source(here::here("config", "config.R"))
source(here::here("config", "utils.R"))

setup_logging("5_tables_MR_analysis")
log_info("=== Starting table generation for endoMR-PREG study ===")

# Load required packages
load_packages(c("dplyr", "stringr", "readr", "openxlsx", "knitr", "data.table", "purrr"))

# ---- Load harmonised data ----
harm_file <- file.path(PATHS$results, "harmonised_rahmioglu_bpo.csv")
if (!file.exists(harm_file)) {
  log_error("Harmonised data not found. Please run harmonisation script first.")
  stop("Harmonised data not found.")
}

dat <- data.table::fread(harm_file)
log_info("Loaded harmonised data with", nrow(dat), "rows")

# ---- Whitelist 41 SNPs ----
# Load from clumped data
snps_file <- file.path(PATHS$results, "endometriosis_clumped_snps.tsv")
if (file.exists(snps_file)) {
  clumped_snps <- data.table::fread(snps_file)
  snp_whitelist <- unique(clumped_snps$SNP)
} else {
  # Use unique SNPs from harmonised data
  snp_whitelist <- unique(dat$SNP)
}

log_info("Using", length(snp_whitelist), "SNP instruments")

# ---- Outcome labels for pregnancy traits ----
outcome_labels <- c(
  # FinnGen
  "finngen_R12_N14_FEMALEINFERT_filtered"     = "Female infertility",
  "finngen_R12_O15_PLAC_DISORD_filtered"      = "Placental disorders (overall)",
  "finngen_R12_O15_PLAC_PRAEVIA_filtered"     = "Placenta praevia",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered"= "Placental abruption",
  "finngen_R12_O15_PREG_ECTOP_filtered"       = "Ectopic pregnancy",
  
  # MR-PREG
  "anaemia_preg_all"   = "Anaemia in pregnancy",
  "apgar1"             = "Apgar score at 1 min",
  "apgar5"             = "Apgar score at 5 min",
  "bf_dur_4c"          = "Breastfeeding ≥4 months",
  "bf_est"             = "Exclusive breastfeeding",
  "bf_ini"             = "Breastfeeding initiation",
  "bf_sus"             = "Breastfeeding cessation",
  "cs"                 = "Caesarean section",
  "el_cs"              = "Elective caesarean section",
  "em_cs"              = "Emergency caesarean section",
  "depr_subsamp"       = "Postpartum depression",
  "ga_all"             = "Gestational age",
  "ga_subsamp"         = "Gestational age (subset)",
  "gdm_subsamp"        = "Gestational diabetes",
  "gh_subsamp"         = "Gestational hypertension",
  "hdp_subsamp"        = "Hypertensive disorders",
  "hbw_all"            = "High birthweight (>4000 g)",
  "lbw_all"            = "Low birthweight (<2500 g)",
  "lga"                = "Large for gestational age",
  "sga"                = "Small for gestational age",
  "hyp"                = "Hyperemesis gravidarum",
  "induction"          = "Induction of labour",
  "lowapgar1"          = "Low Apgar score at 1 min",
  "lowapgar5"          = "Low Apgar score at 5 min",
  "misc_subsamp"       = "Miscarriage",
  "nicu"               = "NICU admission",
  "nvp_sev_all"        = "Severe nausea/vomiting",
  "nvp_sev_subsamp"    = "Severe nausea/vomiting (subset)",
  "pe_subsamp"         = "Preeclampsia",
  "posttb_all"         = "Post-term birth",
  "pretb_all"          = "Preterm birth (any)",
  "pretb_subsamp"      = "Preterm birth (spontaneous)",
  "vpretb_all"         = "Very preterm birth",
  "r_misc_subsamp"     = "Recurrent miscarriage",
  "s_misc_subsamp"     = "Single miscarriage",
  "sb_subsamp"         = "Stillbirth",
  "rup_memb"           = "Premature rupture of membranes",
  "zbw_all"            = "Birthweight Z-score",
  
  # Westergaard PPH
  "Antepartum_bleeding_filtered"                 = "Antepartum bleeding",
  "Early_bleeding_ending_in_live_birth_filtered" = "Early bleeding (ending in live birth)",
  "Early_bleeding_with_any_outcome_filtered"     = "Early bleeding (any outcome)",
  "Postpartum_hemorrhage_due_to_atony_filtered"  = "PPH (uterine atony)",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered" = "PPH (retained placenta)",
  "Postpartum_hemorrhage_filtered"               = "Postpartum hemorrhage"
)

# ---- Pregnancy outcomes summary ----
log_info("Creating sample size summary for pregnancy outcomes...")

# Load outcome data and create summary by data source
preg_summary <- dat %>%
  filter(SNP %in% snp_whitelist) %>%
  group_by(outcome) %>%
  summarise(
    `No. SNPs` = n_distinct(SNP),
    `N total`  = first(samplesize.outcome),
    Cases      = NA_integer_,
    Controls   = NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(
    Source = case_when(
      grepl("^finngen_R12_", outcome) ~ "FinnGen R12",
      grepl("Postpartum|Antepartum|Early_bleeding", outcome) ~ "Westergaard (PPH)",
      TRUE ~ "MR-PREG"
    )
  )

# ---- Fill FinnGen cases/controls from R12 manifest (robust) ----
R12_MANIFEST_URL <- "https://zenodo.org/records/15815805/files/finngen_R12_manifest.csv"

fg_endpoints <- unique(preg_summary$outcome[grepl("^finngen_R12_", preg_summary$outcome)])
# Remove _filtered suffix for manifest lookup
fg_endpoints_clean <- gsub("_filtered$", "", fg_endpoints)

# Read manifest
if (length(fg_endpoints) > 0) {
  log_info("Fetching FinnGen R12 sample sizes...")
  mani <- tryCatch(
    data.table::fread(R12_MANIFEST_URL, showProgress = FALSE),
    error = function(e) {
      log_warn("Could not fetch FinnGen R12 manifest:", e$message)
      NULL
    }
  )
  
  if (!is.null(mani)) {
    # Find phenocode column and URL column (robust)
    phenocode_col <- c("phenocode","phenotype","code")
    phenocode_col <- phenocode_col[phenocode_col %in% names(mani)][1]
    
    if (!is.na(phenocode_col)) {
      url_like_cols <- names(mani)[grepl("url|path|file", names(mani), ignore.case = TRUE)]
      url_col <- if (length(url_like_cols)) url_like_cols[1] else NA_character_
      
      # Standardize phenocode with finngen_R12_ prefix
      mani <- mani %>%
        mutate(phenocode_std = if_else(
          str_starts(.data[[phenocode_col]], "finngen_R12_"),
          .data[[phenocode_col]],
          paste0("finngen_R12_", .data[[phenocode_col]])
        ))
      
      # Build targets with URL from manifest
      if (!is.na(url_col)) {
        targets <- mani %>%
          filter(phenocode_std %in% fg_endpoints_clean) %>%
          transmute(
            phenocode = phenocode_std,
            url = .data[[url_col]]
          ) %>%
          distinct()
      } else {
        # fallback: construct standard R12 path
        targets <- mani %>%
          filter(phenocode_std %in% fg_endpoints_clean) %>%
          transmute(
            phenocode = phenocode_std,
            url = paste0(
              "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/",
              phenocode_std, ".gz"
            )
          ) %>%
          distinct()
      }
      
      # Read each summary stats and extract n_cases / n_controls
      read_remote_gz <- function(url) {
        tryCatch(
          data.table::fread(url, nThread = max(1, parallel::detectCores()-1), showProgress = FALSE),
          error = function(e) {
            log_warn("Failed to read", basename(url))
            NULL
          }
        )
      }
      
      get_sizes <- function(ep, url) {
        dt <- read_remote_gz(url)
        if (is.null(dt) || !nrow(dt)) return(tibble(outcome = ep, Cases = NA_integer_, Controls = NA_integer_))
        
        # be tolerant to naming variants
        ncase_col <- c("n_cases","ncase","N_cases"); ncase_col <- ncase_col[ncase_col %in% names(dt)][1]
        nctrl_col <- c("n_controls","ncontrol","N_controls"); nctrl_col <- nctrl_col[nctrl_col %in% names(dt)][1]
        
        nc <- if (!is.na(ncase_col)) suppressWarnings(max(dt[[ncase_col]], na.rm = TRUE)) else NA
        nt <- if (!is.na(nctrl_col)) suppressWarnings(max(dt[[nctrl_col]], na.rm = TRUE)) else NA
        
        if (!is.finite(nc)) nc <- NA_integer_
        if (!is.finite(nt)) nt <- NA_integer_
        
        tibble(outcome = ep, Cases = as.integer(nc), Controls = as.integer(nt))
      }
      
      # Try to get FinnGen sample sizes (with error handling)
      if (nrow(targets) > 0) {
        fg_sizes <- tryCatch({
          purrr::map2_dfr(targets$phenocode, targets$url, get_sizes)
        }, error = function(e) {
          log_warn("Error downloading FinnGen data, using NA values")
          data.frame(outcome = character(0), Cases = integer(0), Controls = integer(0))
        })
        
        # Join back to preg_summary (add _filtered suffix back)
        if (nrow(fg_sizes) > 0) {
          fg_sizes_filtered <- fg_sizes %>%
            mutate(outcome_filtered = paste0(outcome, "_filtered"))
          
          preg_summary <- preg_summary %>%
            left_join(fg_sizes_filtered, by = c("outcome" = "outcome_filtered")) %>%
            mutate(
              Cases    = dplyr::coalesce(Cases.y, Cases.x),
              Controls = dplyr::coalesce(Controls.y, Controls.x),
              `N total` = ifelse(!is.na(Cases) & !is.na(Controls), Cases + Controls, `N total`)
            ) %>%
            select(-Cases.x, -Controls.x, -Cases.y, -Controls.y)
        }
      }
    }
  }
}

# ---- Add known sample sizes for MR-PREG outcomes ----
# Extract case/control info where available from the harmonised data
case_control_info <- dat %>%
  filter(!is.na(ncase.outcome) | !is.na(ncontrol.outcome)) %>%
  group_by(outcome) %>%
  summarise(
    Cases_extracted = first(ncase.outcome),
    Controls_extracted = first(ncontrol.outcome),
    .groups = "drop"
  )

preg_summary <- preg_summary %>%
  left_join(case_control_info, by = "outcome") %>%
  mutate(
    Cases = dplyr::coalesce(Cases, Cases_extracted),
    Controls = dplyr::coalesce(Controls, Controls_extracted),
    `N total` = ifelse(!is.na(Cases) & !is.na(Controls) & is.na(`N total`), 
                       Cases + Controls, `N total`)
  ) %>%
  select(-Cases_extracted, -Controls_extracted)

# ---- Combine & tidy ----
table1_preg <- preg_summary %>%
  select(outcome, Source, `No. SNPs`, `N total`, Cases, Controls) %>%
  rename(Outcome = outcome) %>%
  mutate(Outcome = dplyr::recode(Outcome, !!!outcome_labels)) %>%
  arrange(Source, Outcome)

# ---- Save & preview ----
table1_file <- file.path(PATHS$results, "Table1_pregnancy_sample_sizes.csv")
write.csv(table1_preg, table1_file, row.names = FALSE)
log_info("Table 1 saved:", table1_file)

if (requireNamespace("knitr", quietly = TRUE)) {
  cat("\n=== TABLE 1: SAMPLE SIZES ===\n")
  print(knitr::kable(table1_preg, align = "lrrrr", caption = "Table 1. Sample sizes and SNP counts for pregnancy outcomes"))
} else {
  print(table1_preg, n = nrow(table1_preg))
}

###############################################
# Table 2. Primary Mendelian Randomization (IVW) Estimates
# for the Association of Genetically Predicted Endometriosis
# With Pregnancy Outcomes
###############################################

# Load IVW results
ivw_file <- file.path(PATHS$results, "ivw_results.csv")
if (!file.exists(ivw_file)) {
  log_warn("IVW results not found. Please run MR analysis first.")
} else {
  log_info("Creating Table 2: Primary IVW estimates...")
  
  # Load results
  ivw_res <- read.csv(ivw_file, stringsAsFactors = FALSE)
  
  # Format and rename for publication
  ivw_tab_pub <- ivw_res %>%
    filter(method == "Inverse variance weighted") %>%
    mutate(
      `Data source` = case_when(
        grepl("^finngen_R12_", outcome) ~ "FinnGen R12",
        grepl("Postpartum|Antepartum|Early_bleeding", outcome) ~ "Westergaard (PPH)",
        TRUE ~ "MR-PREG"
      ),
      OR_val = exp(b),
      CI_low = exp(b - 1.96 * se),
      CI_high = exp(b + 1.96 * se),
      `OR` = sprintf("%.2f", OR_val),
      `95% CI` = sprintf("(%.2f–%.2f)", CI_low, CI_high),
      `P` = case_when(
        pval < 0.001 ~ sprintf("%.2e", pval),
        pval < 0.01 ~ sprintf("%.3f", pval),
        TRUE ~ sprintf("%.2f", pval)
      ),
      `P (Bonferroni)` = ifelse(pval * n() < 0.05, "<0.05", sprintf("%.2e", pval * n())),
      `Q (FDR)` = sprintf("%.2e", p.adjust(pval, method = "fdr")),
      `Significance` = case_when(
        pval < 0.05 & pval * n() < 0.05 ~ "Yes (Bonf)",
        pval < 0.05 ~ "Yes",
        TRUE ~ "No"
      ),
      Outcome = recode(outcome, !!!outcome_labels)
    ) %>%
    select(
      `Data source`, Outcome, `No. SNPs` = nsnp, `OR`, `95% CI`, `P`, `P (Bonferroni)`, `Q (FDR)`, `Significance`
    ) %>%
    arrange(`Data source`, Outcome)
  
  # Preview
  cat("\n=== TABLE 2: PRIMARY IVW ESTIMATES ===\n")
  if (requireNamespace("knitr", quietly = TRUE)) {
    print(knitr::kable(ivw_tab_pub, align = "llrrrrrrr", caption = "Table 2. Primary IVW Estimates"))
  } else {
    print(as.data.frame(ivw_tab_pub), row.names = FALSE)
  }
  
  # Export
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    xlsx_file <- file.path(PATHS$results, "ivw_results_pregnancy_table.xlsx")
    write.xlsx(ivw_tab_pub, file = xlsx_file, overwrite = TRUE)
    log_info("Table 2 (Excel) saved:", xlsx_file)
  }
  
  csv_file <- file.path(PATHS$results, "ivw_results_pregnancy_table.csv")
  write.csv(ivw_tab_pub, csv_file, row.names = FALSE)
  log_info("Table 2 (CSV) saved:", csv_file)

  ###############################################
  # Table 3. Sensitivity Analyses Summary (IVW / WM / Egger)
  ###############################################
  
  # Load other MR results
  wm_file <- file.path(PATHS$results, "weighted_median_results.csv")
  egger_file <- file.path(PATHS$results, "egger_results.csv")
  
  if (file.exists(wm_file) && file.exists(egger_file)) {
    log_info("Creating Table 3: Sensitivity analyses summary...")
    
    wm_res <- read.csv(wm_file, stringsAsFactors = FALSE)
    egger_res <- read.csv(egger_file, stringsAsFactors = FALSE)
    
    # Helper formatters
    fmt_or   <- function(b) sprintf("%.2f", exp(b))
    fmt_ci   <- function(b,se) paste0(sprintf("%.2f", exp(b - 1.96*se)),
                                      "–",
                                      sprintf("%.2f", exp(b + 1.96*se)))
    fmt_p    <- function(p)  {
      case_when(
        p < 0.001 ~ formatC(p, format = "e", digits = 2),
        p < 0.01 ~ sprintf("%.3f", p),
        TRUE ~ sprintf("%.2f", p)
      )
    }
    
    # Prepare method-specific slim tables
    ivw_slim <- ivw_res %>%
      transmute(
        outcome,
        `Data source` = case_when(
          grepl("^finngen_R12_", outcome) ~ "FinnGen R12",
          grepl("Postpartum|Antepartum|Early_bleeding", outcome) ~ "Westergaard (PPH)",
          TRUE ~ "MR-PREG"
        ),
        SNPs = nsnp,
        IVW_OR   = fmt_or(b),
        IVW_CI   = fmt_ci(b, se),
        IVW_P    = fmt_p(pval)
      )
    
    wm_slim <- wm_res %>%
      transmute(
        outcome,
        WM_OR    = fmt_or(b),
        WM_CI    = fmt_ci(b, se),
        WM_P     = fmt_p(pval)
      )
    
    egger_slim <- egger_res %>%
      transmute(
        outcome,
        Egger_OR = fmt_or(b),
        Egger_CI = fmt_ci(b, se),
        Egger_P  = fmt_p(pval)
      )
    
    # Join all and label outcomes
    sens_tab <- ivw_slim %>%
      left_join(wm_slim,    by = "outcome") %>%
      left_join(egger_slim, by = "outcome") %>%
      mutate(
        Outcome = recode(outcome, !!!outcome_labels)
      ) %>%
      select(`Data source`, Outcome, SNPs,
             `IVW OR` = IVW_OR, `IVW 95% CI` = IVW_CI, `IVW P` = IVW_P,
             `WM OR`  = WM_OR,  `WM 95% CI`  = WM_CI,  `WM P`  = WM_P,
             `MR-Egger OR` = Egger_OR, `MR-Egger 95% CI` = Egger_CI, `MR-Egger P` = Egger_P) %>%
      arrange(`Data source`, Outcome)
    
    # Print to console
    sens_tab <- tibble::as_tibble(sens_tab)
    cat("\n=== TABLE 3: SENSITIVITY ANALYSES (IVW, Weighted Median, MR-Egger) ===\n")
    if (requireNamespace("knitr", quietly = TRUE)) {
      print(knitr::kable(sens_tab, align = "lllrrrrrrrr", caption = "Table 3. Sensitivity Analyses Summary"))
    } else {
      print(sens_tab, n = Inf, width = Inf, na.print = "")
    }
    
    # Export to Excel and CSV
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      xlsx_sens <- file.path(PATHS$results, "mr_sensitivity_pregnancy_summary.xlsx")
      write.xlsx(sens_tab, xlsx_sens, asTable = TRUE)
      log_info("Table 3 (Excel) saved:", xlsx_sens)
    }
    
    csv_sens <- file.path(PATHS$results, "mr_sensitivity_pregnancy_summary.csv")
    write.csv(sens_tab, csv_sens, row.names = FALSE)
    log_info("Table 3 (CSV) saved:", csv_sens)
  } else {
    log_warn("Weighted median or Egger results not found. Skipping sensitivity table.")
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Table 4. Comprehensive MR Results Summary
# ─────────────────────────────────────────────────────────────────────────────

if (exists("ivw_res")) {
  log_info("Creating Table 4: Comprehensive results summary...")
  
  # Create comprehensive results summary
  comprehensive_summary <- ivw_res %>%
    filter(method == "Inverse variance weighted") %>%
    mutate(
      `Data source` = case_when(
        grepl("^finngen_R12_", outcome) ~ "FinnGen R12",
        grepl("Postpartum|Antepartum|Early_bleeding", outcome) ~ "Westergaard (PPH)",
        TRUE ~ "MR-PREG"
      ),
      OR = exp(b),
      CI_lower = exp(b - 1.96 * se),
      CI_upper = exp(b + 1.96 * se),
      # Significance levels
      Significance_level = case_when(
        pval < 0.001 ~ "p < 0.001",
        pval < 0.01 ~ "p < 0.01",
        pval < 0.05 ~ "p < 0.05",
        pval < 0.10 ~ "p < 0.10",
        TRUE ~ "p ≥ 0.10"
      ),
      # Effect size classification
      Effect_magnitude = case_when(
        OR > 2.0 ~ "Large increase (OR > 2.0)",
        OR > 1.5 ~ "Moderate increase (1.5 < OR ≤ 2.0)",
        OR > 1.2 ~ "Small increase (1.2 < OR ≤ 1.5)",
        OR > 1.05 ~ "Minimal increase (1.05 < OR ≤ 1.2)",
        OR >= 0.95 ~ "No effect (0.95 ≤ OR ≤ 1.05)",
        OR >= 0.83 ~ "Minimal decrease (0.83 ≤ OR < 0.95)",
        OR >= 0.67 ~ "Small decrease (0.67 ≤ OR < 0.83)",
        OR >= 0.50 ~ "Moderate decrease (0.50 ≤ OR < 0.67)",
        TRUE ~ "Large decrease (OR < 0.50)"
      ),
      # Clinical interpretation
      Clinical_relevance = case_when(
        pval < 0.05 & (OR > 1.2 | OR < 0.83) ~ "Clinically significant",
        pval < 0.05 ~ "Statistically significant",
        pval < 0.10 ~ "Suggestive evidence",
        TRUE ~ "No evidence"
      ),
      Outcome_label = recode(outcome, !!!outcome_labels)
    ) %>%
    select(
      `Data source`,
      Outcome = outcome,
      Outcome_label,
      `N SNPs` = nsnp,
      Beta = b,
      SE = se,
      OR,
      `CI lower` = CI_lower,
      `CI upper` = CI_upper,
      `P value` = pval,
      Significance_level,
      Effect_magnitude,
      Clinical_relevance
    ) %>%
    arrange(`Data source`, `P value`)
  
  # Save comprehensive summary
  comp_file <- file.path(PATHS$results, "comprehensive_mr_pregnancy_summary.csv")
  write.csv(comprehensive_summary, comp_file, row.names = FALSE)
  log_info("Table 4 saved:", comp_file)
  
  # Create significance breakdown
  sig_breakdown <- comprehensive_summary %>%
    group_by(`Data source`, Significance_level) %>%
    summarise(N_outcomes = n(), .groups = "drop") %>%
    arrange(`Data source`, Significance_level)
  
  # Create effect size breakdown
  effect_breakdown <- comprehensive_summary %>%
    group_by(`Data source`, Effect_magnitude) %>%
    summarise(N_outcomes = n(), Mean_OR = round(mean(OR), 2), .groups = "drop") %>%
    arrange(`Data source`, desc(Mean_OR))
  
  # Save breakdowns
  write.csv(sig_breakdown, file.path(PATHS$results, "significance_breakdown.csv"), row.names = FALSE)
  write.csv(effect_breakdown, file.path(PATHS$results, "effect_size_breakdown.csv"), row.names = FALSE)
  
  # Print summary statistics to console
  cat("\n=== MR ANALYSIS OVERVIEW ===\n")
  cat("Total outcomes analyzed:", nrow(comprehensive_summary), "\n")
  
  # Count significant results
  sig_results <- comprehensive_summary %>% filter(`P value` < 0.05)
  cat("Significant associations (p < 0.05):", nrow(sig_results), "\n")
  
  # Results by data source
  source_summary <- comprehensive_summary %>%
    group_by(`Data source`) %>%
    summarise(
      Total = n(),
      Significant = sum(`P value` < 0.05),
      `% Significant` = round(100 * Significant / Total, 1),
      .groups = "drop"
    )
  
  cat("\n--- Results by Data Source ---\n")
  print(source_summary)
  
  # Top 10 most significant results
  top_results <- comprehensive_summary %>%
    arrange(`P value`) %>%
    head(10) %>%
    select(`Data source`, Outcome_label, OR, `P value`, Effect_magnitude)
  
  cat("\n--- Top 10 Most Significant Associations ---\n")
  print(top_results, n = 10)
  
  # Publication-ready significant results table
  publication_significant <- comprehensive_summary %>%
    filter(`P value` < 0.05) %>%
    mutate(
      `OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", OR, `CI lower`, `CI upper`),
      `P-value` = case_when(
        `P value` < 0.001 ~ "<0.001",
        `P value` < 0.01 ~ sprintf("%.3f", `P value`),
        TRUE ~ sprintf("%.2f", `P value`)
      )
    ) %>%
    select(
      `Data Source` = `Data source`,
      Outcome = Outcome_label,
      `SNPs` = `N SNPs`,
      `OR (95% CI)`,
      `P-value`,
      `Effect Classification` = Effect_magnitude,
      `Clinical Interpretation` = Clinical_relevance
    ) %>%
    arrange(`Data Source`, `P-value`)
  
  # Save publication table
  pub_file <- file.path(PATHS$results, "significant_associations_publication_ready.csv")
  write.csv(publication_significant, pub_file, row.names = FALSE)
  
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    pub_xlsx <- file.path(PATHS$results, "significant_associations_publication_ready.xlsx")
    write.xlsx(publication_significant, pub_xlsx, overwrite = TRUE)
    log_info("Publication table (Excel) saved:", pub_xlsx)
  }
  
  log_info("Publication table (CSV) saved:", pub_file)
}

# ─────────────────────────────────────────────────────────────────────────────
# Study-specific methodological notes
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== METHODOLOGICAL NOTES FOR ENDOMETRIOSIS → PREGNANCY OUTCOMES ===\n")
cat("Study Design: Two-sample Mendelian Randomization\n")
cat("Exposure: Endometriosis (41 genome-wide significant SNPs from Rahmioglu et al.)\n")
cat("Outcomes: 49 pregnancy and reproductive outcomes from multiple sources\n")
cat("\nData Sources:\n")
cat("- MR-PREG: Comprehensive pregnancy outcomes consortium\n")
cat("- FinnGen R12: Finnish population-based biobank\n") 
cat("- Westergaard (PPH): Postpartum hemorrhage subtypes\n")
cat("\nMR Methods:\n")
cat("- Primary: Inverse variance weighted (IVW)\n")
cat("- Sensitivity: Weighted median, MR-Egger\n")
cat("- Multiple testing correction: Bonferroni and FDR\n")
cat("==========================================\n")

log_info("=== TABLE GENERATION COMPLETE ===")
log_info("Tables saved to:", PATHS$results)
log_info("\nKey files generated:")
log_info("- Table1_pregnancy_sample_sizes.csv")
log_info("- ivw_results_pregnancy_table.xlsx/.csv")
log_info("- mr_sensitivity_pregnancy_summary.xlsx/.csv")
log_info("- comprehensive_mr_pregnancy_summary.csv")
log_info("- significant_associations_publication_ready.xlsx/.csv")