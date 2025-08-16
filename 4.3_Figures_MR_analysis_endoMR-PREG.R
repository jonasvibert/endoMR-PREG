###############################################################################
# MR Scatter Plots for Significant Outcomes (JAMA-style figures)
# - Generates scatter plots for significant IVW results
# - Adds human-readable titles and effect size subtitles
# - Skips gracefully if insufficient SNPs or missing results
# - Saves high-resolution PNGs for manuscript submission
###############################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# ===== Title: Scatter plots with IVW + MR-Egger + Weighted Median + Modes =====
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  # # # # #

library(TwoSampleMR)
library(dplyr)
library(purrr)
library(ggplot2)

# 1) Make sure your harmonised data are loaded
# dat <- read.csv(harm_file)  # or data.table::fread(harm_file)

# 2) Choose the MR methods to run
method_list <- c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_weighted_mode",
  "mr_simple_mode"
)

# 3) Re-run MR so that res includes all these methods
res_all <- mr(dat, method_list = method_list)

# 4) Outcomes you want to plot
sig_outcomes <- c(
  "finngen_R12_N14_FEMALEINFERT_filtered",   # Female infertility
  "5sBPKR",                                # Premature rupture of membranes
  "finngen_R12_O15_PLAC_PRAEVIA_filtered",   # Placenta praevia
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered" # Placental abruption
)

# (Optional) pretty labels for plot titles
pretty_names <- c(
  finngen_R12_N14_FEMALEINFERT_filtered = "Female infertility (FinnGen)",
  "5sBPKR"     = "Premature rupture of membranes",
  finngen_R12_O15_PLAC_PRAEVIA_filtered  = "Placenta praevia (FinnGen)",
  finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered = "Placental abruption (FinnGen)"
)

# 5) Create output folder
out_dir <- "mr_scatter_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 6) Loop over outcomes: filter results to those outcomes (keep all methods),
#    then give the same subset of harmonised data to mr_scatter_plot
plots <- map(sig_outcomes, function(o){
  
  # rows in the harmonised dataset for this outcome
  dat_o <- dat %>% filter(id.outcome == o)
  
  # MR estimates for this outcome & all requested methods (drop empty ones)
  res_o <- res_all %>%
    filter(id.outcome == o, method %in% c(
      "Inverse variance weighted",
      "MR Egger",
      "Weighted median",
      "Simple mode",
      "Weighted mode"
    )) %>%
    arrange(method)
  
  # If nothing to plot (e.g., <3 SNPs so Egger or modes unavailable), skip gracefully
  if (nrow(res_o) == 0 || nrow(dat_o) == 0) return(NULL)
  
  # mr_scatter_plot returns a list of ggplots; for a single exposure–outcome,
  # it’s length 1. We’ll extract [[1]].
  gp_list <- mr_scatter_plot(res_o, dat_o)
  p <- gp_list[[1]] +
    ggtitle(paste0("MR estimate: ", pretty_names[[o]] %||% o)) +
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  # Save both PNG and PDF
  png_file <- file.path(out_dir, paste0(o, "_scatter.png"))
  pdf_file <- file.path(out_dir, paste0(o, "_scatter.pdf"))
  ggsave(png_file, p, width = 6, height = 5, dpi = 300)
  ggsave(pdf_file, p, width = 6, height = 5)
  
  p
})

# 7) (Optional) show which plots were created
names(plots) <- sig_outcomes
plots

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# ===== Forest Plots for Four Key Pregnancy Outcomes in MR Analysis
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# 1) Outcomes (use IDs exactly as in your data)
target_outcomes <- c(
  "finngen_R12_N14_FEMALEINFERT_filtered",   # Female infertility (FinnGen)
  "5sBPKR",                                  # Premature rupture of membranes (PROM)
  "finngen_R12_O15_PLAC_PRAEVIA_filtered",   # Placenta praevia (FinnGen)
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered" # Placental abruption (FinnGen)
)

# 2) Pretty titles for figure headers (use a named list -> names can be any string)
pretty_titles <- list(
  "finngen_R12_N14_FEMALEINFERT_filtered" = "Female infertility (FinnGen)",
  "5sBPKR"                                 = "Premature rupture of membranes",
  "finngen_R12_O15_PLAC_PRAEVIA_filtered"  = "Placenta praevia (FinnGen)",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered" = "Placental abruption (FinnGen)"
)

# 3) Output folder
out_dir <- "mr_forest_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 4) Helper: build & save a forest plot for one outcome
build_save_forest <- function(dat, outcome_id, title_text) {
  dat_o <- dat %>% filter(id.outcome == outcome_id)
  if (nrow(dat_o) == 0) {
    message("No harmonised SNPs for outcome: ", outcome_id, " — skipping.")
    return(NULL)
  }
  # Run singlesnp just for this outcome (keeps the IVW 'All' row)
  res_single_o <- mr_singlesnp(dat_o)
  
  # If still empty (e.g., after internal checks), skip gracefully
  if (nrow(res_single_o) == 0) {
    message("mr_singlesnp returned 0 rows for ", outcome_id, " — skipping.")
    return(NULL)
  }
  
  # Create forest plot (mr_forest_plot returns a list, length 1 here)
  p_list <- mr_forest_plot(res_single_o)
  p <- p_list[[1]] +
    ggtitle(title_text) +
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  # Save PNG + PDF
  png_file <- file.path(out_dir, paste0(outcome_id, "_forest.png"))
  pdf_file <- file.path(out_dir, paste0(outcome_id, "_forest.pdf"))
  ggsave(png_file, p, width = 7, height = 9, dpi = 300)
  ggsave(pdf_file, p, width = 7, height = 9)
  
  return(p)
}

# 5) Loop over the four outcomes and build/save plots
plots <- map(target_outcomes, ~build_save_forest(dat, .x, pretty_titles[[.x]]))
names(plots) <- target_outcomes

# 6) Optional: preview in the viewer (prints only those that were created)
invisible(lapply(plots[!sapply(plots, is.null)], print))