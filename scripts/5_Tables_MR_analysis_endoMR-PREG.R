# =========================
# Table 1: Sample sizes (restrict to 41 SNPs)
# =========================
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(readr)
})

results_dir <- if (exists("results_dir")) results_dir else "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Whitelist 41 SNPs ----
snp_whitelist <- unique(clumped$SNP)
stopifnot(length(snp_whitelist) == 41)

# ---- Labels for outcomes ----
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
  "bf_dur_4c"          = "Breastfeeding duration",
  "bf_est"             = "Breastfeeding ever",
  "bf_ini"             = "Breastfeeding initiation",
  "bf_sus"             = "Breastfeeding sustained",
  "cs"                 = "Caesarean delivery (all)",
  "el_cs"              = "Elective caesarean delivery",
  "em_cs"              = "Emergency caesarean delivery",
  "depr_subsamp"       = "Depression (peripartum)",
  "ga_all"             = "Gestational age (all)",
  "ga_subsamp"         = "Gestational age (subset)",
  "gdm_subsamp"        = "Gestational diabetes",
  "gh_subsamp"         = "Gestational hypertension",
  "hdp_subsamp"        = "Hypertensive disorders (all)",
  "hbw_all"            = "High birthweight (>4000 g)",
  "lbw_all"            = "Low birthweight (<2500 g)",
  "lga"                = "Large for gestational age",
  "sga"                = "Small for gestational age",
  "hyp"                = "Hypothyroidism in pregnancy",
  "induction"          = "Labour induction",
  "lowapgar1"          = "Low Apgar score at 1 min",
  "lowapgar5"          = "Low Apgar score at 5 min",
  "misc_subsamp"       = "Miscarriage (all)",
  "nicu"               = "NICU admission",
  "nvp_sev_all"        = "Nausea/vomiting (all)",
  "nvp_sev_subsamp"    = "Nausea/vomiting (subset)",
  "pe_subsamp"         = "Preeclampsia",
  "posttb_all"         = "Postterm birth",
  "pretb_all"          = "Preterm birth (all)",
  "pretb_subsamp"      = "Preterm birth (subset)",
  "vpretb_all"         = "Very preterm birth",
  "r_misc_subsamp"     = "Recurrent miscarriage",
  "s_misc_subsamp"     = "Secondary miscarriage",
  "sb_subsamp"         = "Stillbirth",
  "rup_memb"           = "Premature rupture of membranes",
  "zbw_all"            = "Z-score birthweight",
  
  # Westergaard
  "Antepartum_bleeding_filtered"                 = "Antepartum bleeding",
  "Early_bleeding_ending_in_live_birth_filtered" = "Early bleeding (ending in live birth)",
  "Early_bleeding_with_any_outcome_filtered"     = "Early bleeding (any outcome)",
  "Postpartum_hemorrhage_due_to_atony_filtered"  = "Postpartum hemorrhage (atony)",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered" = "Postpartum hemorrhage (retained placenta)",
  "Postpartum_hemorrhage_filtered"               = "Postpartum hemorrhage (all)"
)

# ---- MR-PREG summary ----
mrpreg_case_col <- names(mr_preg_dat)[str_detect(names(mr_preg_dat), "^n(case|_case)\\.?outcome$")] |> dplyr::first()
mrpreg_ctrl_col <- names(mr_preg_dat)[str_detect(names(mr_preg_dat), "^n(control|_control)\\.?outcome$")] |> dplyr::first()

mrpreg_summary <- mr_preg_dat %>%
  filter(SNP %in% snp_whitelist) %>%
  group_by(outcome) %>%
  summarise(
    `No. SNPs`  = n_distinct(SNP),
    `N total`   = suppressWarnings(as.integer(median(samplesize.outcome, na.rm = TRUE))),
    Cases       = if (!is.null(mrpreg_case_col)) suppressWarnings(as.integer(median(.data[[mrpreg_case_col]], na.rm = TRUE))) else NA_integer_,
    Controls    = if (!is.null(mrpreg_ctrl_col)) suppressWarnings(as.integer(median(.data[[mrpreg_ctrl_col]], na.rm = TRUE))) else NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(Source = "MR-PREG")

# ---- Westergaard PPH summary ----
pph_raw <- bind_rows(lapply(res_pph, `[[`, "raw"))
pph_summary <- pph_raw %>%
  filter(rsid %in% snp_whitelist) %>%
  group_by(outcome) %>%
  summarise(
    `No. SNPs` = n_distinct(rsid),
    `N total`  = NA_integer_,
    Cases      = NA_integer_,
    Controls   = NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(Source = "Westergaard (PPH)")

# ---- FinnGen summary ----
fg_raw <- bind_rows(lapply(res_fg, `[[`, "raw"))
finngen_summary <- fg_raw %>%
  filter(rsid %in% snp_whitelist) %>%
  group_by(outcome) %>%
  summarise(
    `No. SNPs` = n_distinct(rsid),
    `N total`  = NA_integer_,
    Cases      = NA_integer_,
    Controls   = NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(Source = "FinnGen R12")

# ---- Combine & tidy ----
table1 <- bind_rows(mrpreg_summary, pph_summary, finngen_summary) %>%
  rename(Outcome = outcome) %>%
  mutate(Outcome = dplyr::recode(Outcome, !!!outcome_labels)) %>%
  arrange(Source, Outcome)

# ---- Save & preview ----
write.csv(table1, file.path(results_dir, "Table1_sample_sizes_restricted41.csv"), row.names = FALSE)

if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(table1, align = "lrrrr", caption = "Table 1. Sample sizes and SNP counts (restricted to 41 instruments)"))
} else {
  print(table1, n = nrow(table1))
}


###############################################
# Table 3. Primary Mendelian Randomization (IVW) Estimates
# for the Association of Genetically Predicted Endometriosis
# With Maternal and Fetal Pregnancy Outcomes
###############################################

library(dplyr)
library(openxlsx)

# Clean outcome labels
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
  finngen_R12_N14_FEMALEINFERT_filtered = "Female infertility",
  finngen_R12_O15_PLAC_PRAEVIA_filtered = "Placenta praevia",
  finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered = "Placental abruption",
  finngen_R12_O15_PLAC_DISORD_filtered = "Other placental disorders",
  finngen_R12_O15_PREG_ECTOP_filtered   = "Ectopic pregnancy",
  Early_bleeding_with_any_outcome_filtered = "Early bleeding (any outcome)",
  Postpartum_hemorrhage_due_to_retained_placenta_filtered = "PPH – retained placenta",
  Early_bleeding_ending_in_live_birth_filtered = "Early bleeding – live birth",
  Postpartum_hemorrhage_filtered = "Postpartum haemorrhage",
  Postpartum_hemorrhage_due_to_atony_filtered = "PPH – atony",
  Antepartum_bleeding_filtered = "Antepartum bleeding"
)

# Check object
if (!exists("ivw_res")) stop("Object 'ivw_res' not found")

# Format and rename for publication
ivw_tab_pub <- ivw_res %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    `Data source` = ifelse(grepl("^finngen", outcome), "FinnGen", "MR-PREG"),
    OR_val = exp(b),
    CI_low = exp(b - 1.96 * se),
    CI_high = exp(b + 1.96 * se),
    `OR` = sprintf("%.2f", OR_val),
    `95%% CI` = sprintf("(%.2f–%.2f)", CI_low, CI_high),
    `P` = sprintf("%.2e", pval),
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
    `Data source`, Outcome, `No. SNPs` = nsnp, `OR`, `95%% CI`, `P`, `P (Bonferroni)`, `Q (FDR)`, `Significance`
  )

# Preview
print(ivw_tab_pub)

# Export
write.xlsx(ivw_tab_pub, file = "results/ivw_results_table3.xlsx", overwrite = TRUE)


###############################################
# Sensitivity Analyses Summary (IVW / WM / Egger)
# - uses existing ivw_res, wm_res, egger_res
# - formats OR (95% CI) and P
# - applies readable outcome labels
# - prints to console + exports to Excel
###############################################

library(dplyr)
library(stringr)
library(tibble)
library(openxlsx)

# ---- Helper formatters ----
fmt_or   <- function(b) sprintf("%.2f", exp(b))
fmt_ci   <- function(b,se) paste0(sprintf("%.2f", exp(b - 1.96*se)),
                                  "–",
                                  sprintf("%.2f", exp(b + 1.96*se)))
fmt_p    <- function(p)  formatC(p, format = "e", digits = 2)

# ---- Outcome labels (your dictionary) ----
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
  finngen_R12_N14_FEMALEINFERT_filtered = "Female infertility",
  finngen_R12_O15_PLAC_PRAEVIA_filtered = "Placenta praevia",
  finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered = "Placental abruption",
  finngen_R12_O15_PLAC_DISORD_filtered = "Other placental disorders",
  finngen_R12_O15_PREG_ECTOP_filtered   = "Ectopic pregnancy",
  Early_bleeding_with_any_outcome_filtered = "Early bleeding (any outcome)",
  Postpartum_hemorrhage_due_to_retained_placenta_filtered = "PPH – retained placenta",
  Early_bleeding_ending_in_live_birth_filtered = "Early bleeding – live birth",
  Postpartum_hemorrhage_filtered = "Postpartum haemorrhage",
  Postpartum_hemorrhage_due_to_atony_filtered = "PPH – atony",
  Antepartum_bleeding_filtered = "Antepartum bleeding"
)

# Optional: data source detector
detect_source <- function(x){
  case_when(
    str_detect(x, "^finngen_R12_") ~ "FinnGen R12",
    str_detect(x, "Postpartum|Antepartum|Early_bleeding") ~ "Westergaard (PPH)",
    TRUE ~ "MR-PREG"
  )
}

# ---- Prepare method-specific slim tables (no recompute) ----
ivw_slim <- ivw_res %>%
  transmute(
    outcome,
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


# ---- Join all and label outcomes ----
sens_tab <- ivw_slim %>%
  left_join(wm_slim,    by = "outcome") %>%
  left_join(egger_slim, by = "outcome") %>%
  mutate(
    `Data source` = detect_source(outcome),
    Outcome = recode(outcome, !!!outcome_labels)
  ) %>%
  select(`Data source`, Outcome, SNPs,
         `IVW OR` = IVW_OR, `IVW 95% CI` = IVW_CI, `IVW P` = IVW_P,
         `WM OR`  = WM_OR,  `WM 95% CI`  = WM_CI,  `WM P`  = WM_P,
         `MR-Egger OR` = Egger_OR, `MR-Egger 95% CI` = Egger_CI, `MR-Egger P` = Egger_P) %>%
  arrange(factor(`Data source`, levels = c("MR-PREG","FinnGen R12","Westergaard (PPH)")),
          Outcome)

# ---- Print to console ----
options(na.print = "NA")
sens_tab <- tibble::as_tibble(sens_tab)
cat("\n===== Sensitivity Analyses (IVW, Weighted Median, MR-Egger) =====\n")
print(sens_tab, n = Inf, width = Inf)

# ---- Export to Excel (simple) ----
write.xlsx(sens_tab, "mr_sensitivity_summary.xlsx", asTable = TRUE)
cat("\nSaved: mr_sensitivity_summary.xlsx\n")

# ─────────────────────────────────────────────────────────────────────────────
# 7) outcome summary table for MR datasets -------------------------------
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(stringr)

## ---- 1) MR-PREG ----
case_col <- names(mr_preg_dat)[str_detect(names(mr_preg_dat), "^n(case|_case)\\.?outcome$")] %>% dplyr::first()
ctrl_col <- names(mr_preg_dat)[str_detect(names(mr_preg_dat), "^n(control|_control)\\.?outcome$")] %>% dplyr::first()

mrpreg_summary <- mr_preg_dat %>%
  group_by(outcome) %>%
  summarise(
    nsnp         = n_distinct(SNP),
    N_total      = suppressWarnings(as.integer(median(samplesize.outcome, na.rm = TRUE))),
    N_total_min  = suppressWarnings(as.integer(min(samplesize.outcome, na.rm = TRUE))),
    N_total_max  = suppressWarnings(as.integer(max(samplesize.outcome, na.rm = TRUE))),
    cases        = if (!is.null(case_col)) suppressWarnings(as.integer(median(.data[[case_col]], na.rm = TRUE))) else NA_integer_,
    controls     = if (!is.null(ctrl_col)) suppressWarnings(as.integer(median(.data[[ctrl_col]], na.rm = TRUE))) else NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(source = "MR-PREG") %>%
  relocate(source)

## ---- 2) Westergaard PPH ----
pph_summary <- bind_rows(lapply(res_pph, `[[`, "raw")) %>%
  group_by(outcome) %>%
  summarise(
    nsnp    = n_distinct(rsid),
    cohorts = suppressWarnings(as.integer(median(n, na.rm = TRUE))),
    N_total  = NA_integer_,
    cases    = NA_integer_,
    controls = NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(source = "Westergaard (PPH)") %>%
  relocate(source)

## ---- 3) FinnGen ----
finngen_summary <- bind_rows(lapply(res_fg, `[[`, "raw")) %>%
  group_by(outcome) %>%
  summarise(
    nsnp = n_distinct(rsid),
    N_total  = NA_integer_,
    cases    = NA_integer_,
    controls = NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(source = "FinnGen R12") %>%
  relocate(source)

## ---- 4) Combine all ----
outcome_summary <- bind_rows(
  mrpreg_summary %>% select(source, outcome, nsnp, N_total, cases, controls),
  pph_summary %>% select(source, outcome, nsnp, N_total, cases, controls),
  finngen_summary %>% select(source, outcome, nsnp, N_total, cases, controls)
) %>%
  arrange(source, outcome)

## ---- 5) Print directly in console ----
print(outcome_summary, n = nrow(outcome_summary))

## ---- 6) Print ranges for manuscript text ----
mrpreg_range <- mrpreg_summary %>%
  summarise(
    minN = min(N_total_min, na.rm = TRUE),
    maxN = max(N_total_max, na.rm = TRUE)
  ) %>%
  mutate(
    range_text = paste0("N = ", scales::comma(minN), " to N = ", scales::comma(maxN))
  ) %>%
  pull(range_text)

cat("\n--- Section 2.1 MR-PREG ---\n")
cat("Total sample sizes ranged from ", mrpreg_range, ".\n", sep = "")

cat("\n--- Section 2.2 FinnGen ---\n")
cat("Sample sizes to be added from FinnGen metadata.\n")

cat("\n--- Section 2.3 Westergaard PPH ---\n")
cat("Provide subtype-specific totals or use overall study totals (overall PPH ≈ 175,000).\n")


