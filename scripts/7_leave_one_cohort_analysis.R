#!/usr/bin/env Rscript
###############################################################################
# LEAVE-ONE-COHORT-OUT MENDELIAN RANDOMIZATION (MR) ANALYSIS 
# Purpose: Test robustness of MR estimates by excluding each cohort in turn.
###############################################################################

message("=== LEAVE-ONE-COHORT-OUT MR ANALYSIS (refactored, labeled) ===")

### 1) SETUP -----------------------------------------------------------------
message("Loading packages...")
required_pkgs <- c("TwoSampleMR","dplyr","ggplot2","here","readr",
                   "data.table","gridExtra","scales","stringr","tidyr","purrr")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages(lapply(required_pkgs, function(p) library(p, character.only=TRUE)))
set.seed(1)

# Toggle: TRUE = ne plotter que les outcomes significatifs (p<0.05)
PLOT_ONLY_SIGNIFICANT <- TRUE
MAX_IN_COMBINED_PANEL <- 6   # nombre max de panneaux dans l’image combinée

results_dir <- here::here("results")
plots_dir   <- here::here("plots")
dir.create(results_dir, FALSE, TRUE); dir.create(plots_dir, FALSE, TRUE)

outcome_file <- here::here("data","OUTCOME_MR-PREG","stu_out_dat.txt")
harm_file    <- here::here("results","harmonised_rahmioglu_bpo.csv")

### 2) LOAD & VALIDATE -------------------------------------------------------
message("Loading outcome data with cohort information...")
if (!file.exists(outcome_file)) stop("Outcome file not found: ", outcome_file)
if (!file.exists(harm_file))    stop("Harmonized exposure file not found: ", harm_file)

mr_data <- readr::read_table(outcome_file, col_types = cols()) |> as.data.frame()

required_cols_out <- c("SNP","beta","se","pval","samplesize","study","Phenotype",
                       "effect_allele","other_allele","eaf","ncase","ncontrol")
miss_out <- setdiff(required_cols_out, names(mr_data))
if (length(miss_out)) stop("Missing required columns in outcome data: ", paste(miss_out, collapse=", "))

harm_data <- data.table::fread(harm_file) |> as.data.frame()
required_cols_exp <- c("SNP","beta.exposure","se.exposure","pval.exposure",
                       "effect_allele.exposure","other_allele.exposure","eaf.exposure",
                       "samplesize.exposure","id.exposure","exposure")
miss_exp <- setdiff(required_cols_exp, names(harm_data))
if (length(miss_exp)) stop("Exposure file lacks required columns: ", paste(miss_exp, collapse=", "))

message("Outcome rows: ", nrow(mr_data), " | Exposure rows (raw): ", nrow(harm_data))

### 3) COHORTS & LABELS ------------------------------------------------------
cohort_summary <- mr_data |>
  group_by(study) |>
  summarise(n_snps=n_distinct(SNP), n_phenotypes=n_distinct(Phenotype),
            mean_samplesize=mean(samplesize, na.rm=TRUE), .groups="drop") |>
  arrange(desc(n_snps))
message("Available cohorts: ", paste(cohort_summary$study, collapse=", "))
if (nrow(cohort_summary) < 2) stop("Need ≥2 cohorts for LOO; found: ", nrow(cohort_summary))

# Human-readable labels
outcome_labels <- c(
  pretb_all="Preterm birth (any)", el_cs="Elective caesarean section",
  rup_memb="Premature rupture of membranes", pretb_subsamp="Preterm birth (spontaneous)",
  apgar1="Apgar score at 1 min", vpretb_all="Very preterm birth (any)",
  em_cs="Emergency caesarean section", lowapgar1="Low Apgar score at 1 min",
  pe_subsamp="Preeclampsia", lbw_all="Low birthweight (<2500 g)",
  ga_all="Gestational age", gh_subsamp="Gestational hypertension",
  ga_subsamp="Gestational age (subset)", lga="Large for gestational age",
  zbw_all="Z-score birthweight", nvp_sev_subsamp="Severe nausea/vomiting (subset)",
  nvp_sev_all="Severe nausea/vomiting", hbw_all="High birthweight (>4000 g)",
  anaemia_preg_all="Pregnancy anaemia",
  finngen_R12_N14_FEMALEINFERT="Female infertility",
  finngen_R12_O15_PLAC_DISORD="Placental disorders",
  finngen_R12_O15_PLAC_PRAEVIA="Placenta praevia",
  finngen_R12_O15_PLAC_PREMAT_SEPAR="Premature placental separation",
  finngen_R12_O15_PREG_ECTOP="Ectopic pregnancy"
)

# Binary/continuous
outcome_types <- c(
  pretb_all="binary", el_cs="binary", rup_memb="binary", pretb_subsamp="binary",
  apgar1="continuous", vpretb_all="binary", em_cs="binary", lowapgar1="binary",
  pe_subsamp="binary", lbw_all="binary", ga_all="continuous", gh_subsamp="binary",
  ga_subsamp="continuous", lga="binary", zbw_all="continuous",
  nvp_sev_subsamp="binary", nvp_sev_all="binary", hbw_all="binary",
  anaemia_preg_all="binary", finngen_R12_N14_FEMALEINFERT="binary",
  finngen_R12_O15_PLAC_DISORD="binary", finngen_R12_O15_PLAC_PRAEVIA="binary",
  finngen_R12_O15_PLAC_PREMAT_SEPAR="binary", finngen_R12_O15_PREG_ECTOP="binary"
)

mr_data <- mr_data |>
  mutate(
    id.outcome    = Phenotype,
    outcome_label = dplyr::coalesce(outcome_labels[Phenotype], Phenotype),
    outcome_type  = dplyr::coalesce(outcome_types[Phenotype], "binary")
  )

### 4) EXPOSURE (UNIQUE SNP ROWS) — WITH phenotype_col -----------------------
message("Preparing exposure dataset...")
exposure_df_raw <- harm_data |>
  select(SNP, beta.exposure, se.exposure, pval.exposure,
         effect_allele.exposure, other_allele.exposure, eaf.exposure,
         samplesize.exposure, id.exposure, exposure) |>
  distinct(SNP, .keep_all=TRUE)

exposure_for_format <- exposure_df_raw |>
  transmute(
    SNP,
    beta = beta.exposure, se = se.exposure, pval = pval.exposure,
    effect_allele = effect_allele.exposure, other_allele = other_allele.exposure,
    eaf = eaf.exposure, samplesize = samplesize.exposure,
    ncase = NA_real_, ncontrol = NA_real_,
    id = id.exposure, exposure = exposure,
    phenotype = exposure                 # <= important: vrai nom pour l’exposition
  )

exposure_formatted <- TwoSampleMR::format_data(
  exposure_for_format, type="exposure", phenotype_col = "phenotype"
)
message("Exposure SNPs: ", dplyr::n_distinct(exposure_formatted$SNP))

### 5) HELPERS: META & HARMONISATION ----------------------------------------
meta_outcome_after_drop <- function(dat_one_outcome) {
  dat_one_outcome |>
    filter(!is.na(SNP), is.finite(beta.outcome), is.finite(se.outcome)) |>
    group_by(SNP) |>
    summarise(
      w = sum(1/(se.outcome^2), na.rm=TRUE),
      beta = ifelse(w > 0, sum(beta.outcome/(se.outcome^2), na.rm=TRUE) / w, NA_real_),
      se   = ifelse(w > 0, sqrt(1/w), NA_real_),
      pval = ifelse(w > 0, 2*pnorm(-abs(beta/se)), NA_real_),
      effect_allele = { tmp <- na.omit(effect_allele.outcome); if (length(tmp)) tmp[1] else NA_character_ },
      other_allele  = { tmp <- na.omit(other_allele.outcome);  if (length(tmp)) tmp[1] else NA_character_ },
      ss_sum = sum(samplesize.outcome, na.rm=TRUE),
      eaf = ifelse(ss_sum > 0 & any(!is.na(eaf.outcome)),
                   sum(eaf.outcome * samplesize.outcome, na.rm=TRUE) / ss_sum,
                   mean(eaf.outcome, na.rm=TRUE)),
      samplesize = ss_sum,
      ncase = sum(ncase.outcome, na.rm=TRUE),
      ncontrol = sum(ncontrol.outcome, na.rm=TRUE),
      .groups = "drop"
    ) |>
    select(-w, -ss_sum) |>
    filter(is.finite(beta), is.finite(se), !is.na(effect_allele), !is.na(other_allele))
}

build_outcome_formatted <- function(meta_df, outcome_id, outcome_label) {
  if (is.null(meta_df) || nrow(meta_df) == 0) return(NULL)
  tmp <- meta_df |>
    transmute(
      SNP, beta, se, pval,
      effect_allele, other_allele, eaf,
      samplesize, ncase, ncontrol,
      id = outcome_id,
      outcome = outcome_label,
      phenotype = outcome_label      # <= important: vrai nom lisible de l’outcome
    ) |>
    filter(!is.na(SNP), is.finite(beta), is.finite(se),
           !is.na(effect_allele), !is.na(other_allele))
  if (nrow(tmp) == 0) return(NULL)
  TwoSampleMR::format_data(tmp, type="outcome", phenotype_col = "phenotype")
}

### 6) MAIN + LEAVE-ONE-COHORT-OUT RUNNER -----------------------------------
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
  
  # MAIN
  main_meta <- meta_outcome_after_drop(one_outcome_df)
  main_outcome_fmt <- build_outcome_formatted(main_meta, outcome_id, outcome_name)
  if (is.null(main_outcome_fmt) || nrow(main_outcome_fmt) < 3) {
    message("   Skipping (main) — insufficient SNPs after meta/format.")
    return(NULL)
  }
  main_harm <- tryCatch(TwoSampleMR::harmonise_data(exposure_formatted, main_outcome_fmt),
                        error=function(e){ message("   Harmonisation error (main): ", e$message); NULL })
  if (is.null(main_harm) || dplyr::n_distinct(main_harm$SNP) < 3) {
    message("   Skipping (main) — <3 SNPs after harmonisation.")
    return(NULL)
  }
  main_ivw <- TwoSampleMR::mr(main_harm, method_list="mr_ivw")
  if (!nrow(main_ivw)) return(NULL)
  main_ivw$analysis_type <- "Main analysis"
  main_ivw$nsnp <- dplyr::n_distinct(main_harm$SNP)
  main_ivw$left_out_cohort <- "None"
  main_ivw$n_snps_remaining <- main_ivw$nsnp
  main_ivw$n_cohorts_remaining <- length(cohorts)
  # for safety, re-attach nice name
  main_ivw$outcome <- outcome_name
  main_ivw$id.outcome <- outcome_id
  
  # LOO
  loo_list <- lapply(cohorts, function(left_out) {
    sub_df <- dplyr::filter(one_outcome_df, study != left_out)
    meta   <- meta_outcome_after_drop(sub_df)
    out_fmt <- build_outcome_formatted(meta, outcome_id, outcome_name)
    if (is.null(out_fmt) || nrow(out_fmt) < 3) {
      message("   Skipping without ", left_out, " — insufficient SNPs after meta/format.")
      return(NULL)
    }
    harm <- tryCatch(TwoSampleMR::harmonise_data(exposure_formatted, out_fmt),
                     error=function(e){ message("   Harmonisation error without ", left_out, ": ", e$message); NULL })
    if (is.null(harm) || dplyr::n_distinct(harm$SNP) < 3) {
      message("   Skipping without ", left_out, " — <3 SNPs after harmonisation.")
      return(NULL)
    }
    res <- TwoSampleMR::mr(harm, method_list="mr_ivw")
    if (!nrow(res)) return(NULL)
    res$analysis_type <- paste0("Without ", left_out)
    res$left_out_cohort <- left_out
    res$n_snps_remaining <- dplyr::n_distinct(harm$SNP)
    res$n_cohorts_remaining <- dplyr::n_distinct(sub_df$study)
    # keep names as well
    res$outcome <- outcome_name
    res$id.outcome <- outcome_id
    res
  })
  
  dplyr::bind_rows(main_ivw, dplyr::bind_rows(loo_list))
}

message("Running main + LOO across outcomes...")
results_list   <- lapply(outcomes_list, run_loo_for_outcome)
all_results_raw <- dplyr::bind_rows(results_list)
if (!nrow(all_results_raw)) stop("No MR results produced. Check data and cohort composition.")

# Assure outcome label column is the readable one
all_results <- all_results_raw |>
  mutate(outcome = dplyr::coalesce(outcome, outcome_labels[id.outcome], id.outcome)) |>
  select(id.outcome, outcome, method, b, se, pval, nsnp,
         analysis_type, left_out_cohort, n_snps_remaining, n_cohorts_remaining) |>
  mutate(
    outcome_type = dplyr::coalesce(outcome_types[id.outcome], "binary"),
    effect_label = ifelse(outcome_type=="binary","OR","Beta"),
    est_display  = ifelse(outcome_type=="binary", exp(b), b),
    ci_low       = ifelse(outcome_type=="binary", exp(b - 1.96*se), b - 1.96*se),
    ci_high      = ifelse(outcome_type=="binary", exp(b + 1.96*se), b + 1.96*se),
    est_ci       = sprintf("%.3f (%.3f–%.3f)", est_display, ci_low, ci_high),
    p_display    = dplyr::case_when(is.na(pval) ~ NA_character_,
                                    pval < 0.001 ~ "<0.001",
                                    pval < 0.01  ~ sprintf("%.3f", pval),
                                    TRUE         ~ sprintf("%.3f", pval))
  )

message("MR results generated: ", nrow(all_results))

### 7) ROBUSTNESS SUMMARY ----------------------------------------------------
message("Computing robustness metrics...")
robustness_summary <- all_results |>
  group_by(id.outcome, outcome, outcome_type) |>
  summarise(
    main_b  = b[analysis_type=="Main analysis"][1],
    main_se = se[analysis_type=="Main analysis"][1],
    main_p  = pval[analysis_type=="Main analysis"][1],
    main_nsnp = nsnp[analysis_type=="Main analysis"][1],
    n_loo_analyses = sum(analysis_type!="Main analysis"),
    main_est_display = est_display[analysis_type=="Main analysis"][1],
    min_est_display  = suppressWarnings(min(est_display[analysis_type!="Main analysis"], na.rm=TRUE)),
    max_est_display  = suppressWarnings(max(est_display[analysis_type!="Main analysis"], na.rm=TRUE)),
    range_est        = max_est_display - min_est_display,
    relative_range   = (range_est / abs(main_est_display)) * 100,
    max_change_percent = suppressWarnings(max(abs(est_display[analysis_type!="Main analysis"] - main_est_display) /
                                                abs(main_est_display) * 100, na.rm=TRUE)),
    n_sig_loo = sum(pval[analysis_type!="Main analysis"] < 0.05, na.rm=TRUE),
    consistency_rate = ifelse(n_loo_analyses>0, n_sig_loo/n_loo_analyses*100, NA_real_),
    robust = dplyr::case_when(!is.na(max_change_percent) & !is.na(consistency_rate) &
                                max_change_percent < 20 & consistency_rate > 80 ~ "High",
                              !is.na(max_change_percent) & !is.na(consistency_rate) &
                                max_change_percent < 50 & consistency_rate > 60 ~ "Moderate",
                              TRUE ~ "Low"),
    .groups="drop"
  ) |>
  arrange(desc(ifelse(outcome_type=="binary", abs(log(main_est_display)), abs(main_est_display))))

### 8) EXPORTS ---------------------------------------------------------------
message("Exporting CSVs...")
readr::write_csv(all_results,        file.path(results_dir, "leave_one_cohort_results.csv"))
readr::write_csv(robustness_summary, file.path(results_dir, "leave_one_cohort_summary.csv"))
readr::write_csv(cohort_summary,     file.path(results_dir, "cohort_composition.csv"))
message("Saved CSVs to: ", results_dir)

### 9) PLOTTING --------------------------------------------------------------
message("Creating forest plots...")
create_forest_plot <- function(df_one_outcome) {
  if (!nrow(df_one_outcome)) return(NULL)
  out_name <- unique(df_one_outcome$outcome)[1]
  out_type <- unique(df_one_outcome$outcome_type)[1]
  
  plot_df <- df_one_outcome |>
    mutate(analysis_label = ifelse(analysis_type=="Main analysis","Main analysis",
                                   paste0("Without ", left_out_cohort)),
           is_main = analysis_type=="Main analysis") |>
    arrange(desc(is_main), est_display)
  
  p <- ggplot(plot_df, aes(x=est_display, y=reorder(analysis_label, est_display))) +
    { if (out_type=="binary") geom_vline(xintercept=1, linetype="dashed") else geom_vline(xintercept=0, linetype="dashed") } +
    geom_point(aes(shape=is_main), size=2.5) +
    geom_errorbarh(aes(xmin=ci_low, xmax=ci_high), height=0.2) +
    { if (out_type=="binary") scale_x_log10(labels=scales::number_format(accuracy=0.01)) else scale_x_continuous() } +
    scale_shape_manual(values=c(`TRUE`=17, `FALSE`=16), labels=c("Leave-one-out","Main analysis")) +
    labs(title=paste0("Leave-one-cohort-out: ", out_name),
         x=ifelse(out_type=="binary","Odds Ratio (95% CI)","Beta (95% CI)"),
         y="Analysis", shape="Point type") +
    theme_minimal() +
    theme(plot.title=element_text(size=12, face="bold"),
          axis.text.y=element_text(size=9),
          legend.position="bottom")
  p
}

if (PLOT_ONLY_SIGNIFICANT) {
  ids_to_plot <- robustness_summary |>
    filter(!is.na(main_p), main_p < 0.05) |>
    slice_head(n = MAX_IN_COMBINED_PANEL) |>
    pull(id.outcome)
} else {
  ids_to_plot <- unique(all_results$id.outcome)
}

plot_list <- list()
if (length(ids_to_plot) > 0) {
  for (oid in ids_to_plot) {
    df_sub <- all_results |> filter(id.outcome == oid)
    p <- create_forest_plot(df_sub)
    if (!is.null(p)) {
      plot_list[[oid]] <- p
      ggsave(file.path(plots_dir, paste0("loo_", oid, ".png")), plot=p, width=7, height=5, dpi=300)
    }
  }
  if (length(plot_list) == 1) {
    ggsave(file.path(plots_dir, "leave_one_cohort_forestplot.png"),
           plot=plot_list[[1]], width=7, height=5, dpi=300)
  } else if (length(plot_list) > 1) {
    combined <- do.call(gridExtra::grid.arrange, c(plot_list, ncol=2))
    ggsave(file.path(plots_dir, "leave_one_cohort_forestplot.png"),
           plot=combined, width=14, height=10, dpi=300)
  }
  message("Saved plots to: ", plots_dir)
} else {
  message("No outcomes selected for plotting — plots skipped.")
}

### 10) SUMMARY LOG ----------------------------------------------------------
message("\n=== SUMMARY ===")
message("Cohorts analysed: ", paste(cohort_summary$study, collapse=", "))
message("Total outcomes evaluated: ", nrow(robustness_summary))
message("Outcomes with main significant effect (p<0.05): ",
        sum(!is.na(robustness_summary$main_p) & robustness_summary$main_p < 0.05))
message("High robustness outcomes: ", sum(robustness_summary$robust=="High", na.rm=TRUE))
message("Moderate robustness outcomes: ", sum(robustness_summary$robust=="Moderate", na.rm=TRUE))
message("Low robustness outcomes: ", sum(robustness_summary$robust=="Low", na.rm=TRUE))

top_robust <- robustness_summary |>
  filter(!is.na(main_p), main_p < 0.05) |>
  arrange(desc(consistency_rate), max_change_percent) |>
  head(5)

if (nrow(top_robust) > 0) {
  message("\nTop robust significant associations:")
  print(top_robust |>
          transmute(
            outcome,
            main_est_display = ifelse(outcome_type=="binary",
                                      sprintf("OR=%.3f", main_est_display),
                                      sprintf("Beta=%.3f", main_est_display)),
            main_p,
            max_change_percent = round(max_change_percent, 1),
            consistency_rate = round(consistency_rate, 1),
            robust
          ))
}

message("\nAnalysis completed successfully.")
message("Results: ", results_dir)
message("Plots:   ", plots_dir)
