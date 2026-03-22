#!/usr/bin/env Rscript
###############################################################################
# Script: 04_main_analyses_endoMR-PREG.R
# Purpose: Two-sample MR analyses on harmonised pregnancy outcomes
# Input  : results/harmonised_rahmioglu_bpo.csv
# Outputs: IVW, Egger, WM, all MR methods, with FDR correction,
#          saved under results/
###############################################################################

### 1) SETUP ###################################################################

required_pkgs <- c(
  "TwoSampleMR", "MRPRESSO", "dplyr", "ggplot2",
  "here", "readr", "data.table", "devtools"
)

safe_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "MRPRESSO") {
      if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools", repos = "https://cloud.r-project.org")
      }
      devtools::install_github("rondolab/MR-PRESSO")
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(required_pkgs, safe_install))

# Paths
results_dir <- here::here("results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

harm_file <- here::here("results", "harmonised_rahmioglu_bpo.csv")

### 2) LOAD + RESTRICT OUTCOMES ################################################

dat <- data.table::fread(harm_file)

stopifnot(all(c("id.exposure", "beta.exposure", "beta.outcome") %in% colnames(dat)))
stopifnot("outcome" %in% colnames(dat))

# ---- 30 retained outcomes ---------------------------------------------------

vars_keep <- c(
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

# ---- Labels -----------------------------------------------------------------
  
  outcome_labels <- c(
    # Bleeding (4)
    Antepartum_bleeding                       = "Antepartum bleeding",
    Postpartum_hemorrhage                     = "Postpartum hemorrhage (any)",
    Postpartum_hemorrhage_due_to_atony        = "PPH due to atony",
    Postpartum_hemorrhage_due_to_retained_placenta = "PPH due to retained placenta",
    
    # Placenta (3)
    finngen_R12_O15_PLAC_PRAEVIA              = "Placenta praevia",
    finngen_R12_O15_PLAC_DISORD               = "Placental disorders",
    finngen_R12_O15_PLAC_PREMAT_SEPAR         = "Premature placental separation",
    
    # Membranes (1)
    rup_memb                                  = "Premature rupture of membranes",
    
    # Birth timing (4)
    pretb_all                                 = "Preterm birth <37 weeks (any)",
    vpretb_all                                = "Very preterm birth <34 weeks",
    posttb_all                                = "Post-term birth",
    ga_all                                    = "Gestational age",
    
    # Fetal Growth (5)
    sga                                       = "Small for gestational age",
    lbw_all                                   = "Low birthweight <2500g",
    hbw_all                                   = "High birthweight >4000g",
    lga                                       = "Large for gestational age",
    zbw_all                                   = "Z-score birthweight",
    
    # Neonatal adaptation / perinatal death (4)
    lowapgar1                                 = "Low Apgar score at 1 min",
    lowapgar5                                 = "Low Apgar score at 5 min",
    nicu                                      = "NICU admission",
    sb_subsamp                                = "Stillbirth",
    
    # Maternal complications (6)
    anaemia_preg_all                          = "Pregnancy anemia",
    gdm_subsamp                               = "Gestational diabetes",
    gh_subsamp                                = "Gestational hypertension",
    hdp_subsamp                               = "Hypertensive disorders of pregnancy",
    pe_subsamp                                = "Preeclampsia",
    depr_subsamp                              = "Postpartum depression",
    
    # Caesarean section (2)
    el_cs                                     = "Elective caesarean section",
    em_cs                                     = "Emergency caesarean section",
    
    # Other obstetric timing (1)
    induction                                 = "Labour induction"
  )
# ---- Restrict dataset -------------------------------------------------------

dat <- dat[dat$outcome %in% vars_keep, , drop = FALSE]

labels_df <- data.frame(
  outcome      = names(outcome_labels),
  outcome_full = unname(outcome_labels),
  stringsAsFactors = FALSE
)

dat <- merge(dat, labels_df, by = "outcome", all.x = TRUE, sort = FALSE)
dat$id.outcome <- dat$outcome

dat <- dat %>%
  dplyr::filter(
    !is.na(beta.outcome),
    !is.na(se.outcome),
    se.outcome > 0,
    !is.na(pval.outcome)
  )

message("Outcomes included (n = ", length(unique(dat$outcome)), ") [expected 30]:")
print(sort(unique(dat$outcome_full)))

### 3) HELPERS ##################################################################

run_mr_methods <- function(data, methods) {
  out <- mr(data, method_list = methods)
  as.data.frame(out)
}

export_csv <- function(df, stem) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    warning("Nothing to export for: ", stem)
    return(invisible(NULL))
  }
  
  if ("outcome" %in% colnames(df) && !"outcome_full" %in% colnames(df)) {
    df <- merge(df, labels_df, by = "outcome", all.x = TRUE, sort = FALSE)
    df$outcome_full[is.na(df$outcome_full)] <- df$outcome[is.na(df$outcome_full)]
  }
  
  num_cols <- vapply(df, is.numeric, logical(1L))
  if (any(num_cols)) {
    df[num_cols] <- lapply(df[num_cols], function(x) {
      sapply(x, function(val) {
        if (is.na(val)) return(NA)
        if (abs(val) < 1e-3) format(val, scientific = TRUE, digits = 3)
        else round(val, 3)
      })
    })
  }
  
  out_path <- file.path(results_dir, paste0(stem, ".csv"))
  write.csv(df, out_path, row.names = FALSE)
  message("Wrote: ", out_path)
}

### 4) MAIN MR ANALYSES #########################################################

ivw_res   <- run_mr_methods(dat, "mr_ivw")
egger_res <- run_mr_methods(dat, "mr_egger_regression")
wm_res    <- run_mr_methods(dat, "mr_weighted_median")
all_res   <- mr(dat)

### 4b) MULTIPLE TESTING FDR #####################################################

ivw_res$qval   <- p.adjust(ivw_res$pval,   method = "fdr")
egger_res$qval <- p.adjust(egger_res$pval, method = "fdr")
wm_res$qval    <- p.adjust(wm_res$pval,    method = "fdr")
# all_res contains all methods combined — FDR not applied (mixed denominators)


# Export
export_csv(ivw_res,   "ivw_results")
export_csv(egger_res, "egger_results")
export_csv(wm_res,    "weighted_median_results")
export_csv(all_res,   "all_mr_methods")
