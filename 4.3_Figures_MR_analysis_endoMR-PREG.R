###############################################################################
# MR figures for “Endometriosis → Fertility & Pregnancy Outcomes”
###############################################################################

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(stringr)
})

if (!requireNamespace("svglite", quietly = TRUE)) {
  try(utils::install.packages("svglite"), silent = TRUE)
}

#========================
# 0) Inputs & dictionaries
#========================

mr_methods <- c("mr_ivw","mr_egger_regression","mr_weighted_median",
                "mr_weighted_mode","mr_simple_mode")

core_ids <- c(
  "finngen_R12_N14_FEMALEINFERT_filtered",
  "5sBPKR",
  "finngen_R12_O15_PLAC_PRAEVIA_filtered",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered"
)

extra_ids <- c(
  "oAOiBt","Rwbyma","zMUPe9","vywWPW",
  "O2sj1s","9QjNDJ","i3cwH7",
  "TTnKdM","mxO2D8","OYxJPL",
  "Postpartum_hemorrhage_filtered",
  "Postpartum_hemorrhage_due_to_atony_filtered",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered",
  "Antepartum_bleeding_filtered",
  "Early_bleeding_with_any_outcome_filtered",
  "Early_bleeding_ending_in_live_birth_filtered",
  "zKkCNG","ycgY62","yoKaH0","iS2FMp","dhBaTu",
  "z2noFK","dBcqKG","wvATQb","RJauGo","6g03jM"
)

label_map <- c(
  "finngen_R12_N14_FEMALEINFERT_filtered"      = "Female infertility (FinnGen)",
  "5sBPKR"                                     = "Premature rupture of membranes",
  "finngen_R12_O15_PLAC_PRAEVIA_filtered"      = "Placenta praevia (FinnGen)",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered" = "Placental abruption (FinnGen)",
  "oAOiBt" = "Hypertensive disorders of pregnancy (subsample)",
  "Rwbyma" = "Preeclampsia (subsample)",
  "zMUPe9" = "Gestational hypertension (subsample)",
  "vywWPW" = "Gestational diabetes (subsample)",
  "O2sj1s" = "Preterm birth (all)",
  "9QjNDJ" = "Preterm birth (subsample)",
  "i3cwH7" = "Very preterm birth (all)",
  "TTnKdM" = "Caesarean section (all)",
  "mxO2D8" = "Elective caesarean section",
  "OYxJPL" = "Emergency caesarean section",
  "Postpartum_hemorrhage_filtered"                          = "Postpartum haemorrhage (overall)",
  "Postpartum_hemorrhage_due_to_atony_filtered"             = "PPH – uterine atony",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered" = "PPH – retained placenta",
  "Antepartum_bleeding_filtered"                            = "Antepartum bleeding",
  "Early_bleeding_with_any_outcome_filtered"                = "Early bleeding (any outcome)",
  "Early_bleeding_ending_in_live_birth_filtered"            = "Early bleeding (ending in live birth)",
  "zKkCNG" = "Low birthweight (<2500 g)",
  "ycgY62" = "Large for gestational age",
  "yoKaH0" = "Small for gestational age",
  "iS2FMp" = "Birthweight Z-score",
  "dhBaTu" = "High birthweight (>4000 g)",
  "z2noFK" = "Low Apgar at 1 min",
  "dBcqKG" = "Low Apgar at 5 min",
  "wvATQb" = "Apgar score at 1 min",
  "RJauGo" = "Apgar score at 5 min",
  "6g03jM" = "NICU admission"
)

group_map <- c(
  "finngen_R12_N14_FEMALEINFERT_filtered"      = "Fertility",
  "5sBPKR"                                     = "Preterm & PROM",
  "finngen_R12_O15_PLAC_PRAEVIA_filtered"      = "Placental disorders",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered" = "Placental disorders",
  "oAOiBt"="Hypertensive disorders","Rwbyma"="Hypertensive disorders","zMUPe9"="Hypertensive disorders",
  "vywWPW"="Metabolic",
  "O2sj1s"="Preterm & PROM","9QjNDJ"="Preterm & PROM","i3cwH7"="Preterm & PROM",
  "TTnKdM"="Delivery mode","mxO2D8"="Delivery mode","OYxJPL"="Delivery mode",
  "Postpartum_hemorrhage_filtered"="Bleeding/PPH",
  "Postpartum_hemorrhage_due_to_atony_filtered"="Bleeding/PPH",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered"="Bleeding/PPH",
  "Antepartum_bleeding_filtered"="Bleeding/PPH",
  "Early_bleeding_with_any_outcome_filtered"="Bleeding/PPH",
  "Early_bleeding_ending_in_live_birth_filtered"="Bleeding/PPH",
  "zKkCNG"="Growth & size","ycgY62"="Growth & size","yoKaH0"="Growth & size",
  "iS2FMp"="Growth & size","dhBaTu"="Growth & size",
  "z2noFK"="Apgar & NICU","dBcqKG"="Apgar & NICU","wvATQb"="Apgar & NICU",
  "RJauGo"="Apgar & NICU","6g03jM"="Apgar & NICU"
)
group_levels <- c("Fertility","Placental disorders","Preterm & PROM","Delivery mode",
                  "Bleeding/PPH","Hypertensive disorders","Metabolic","Growth & size","Apgar & NICU")

#========================
# 1) MR once + per-SNP once
#========================

message("MR across methods…")
res_all <- mr(dat, method_list = mr_methods)

message("Single-SNP effects…")
res_single_all <- mr_singlesnp(dat)

safe_filename <- function(x) {
  x |> str_replace_all("[^A-Za-z0-9_-]+", "_") |> str_squish()
}

#========================
# Output folders (under plot/)
#========================
base_dir   <- "plot"
scatter_dir<- file.path(base_dir, "fig_scatter")
forest_dir <- file.path(base_dir, "fig_forest")
loo_dir    <- file.path(base_dir, "fig_leaveoneout_egger")
funnel_dir <- file.path(base_dir, "fig_funnel")
overview_dir <- file.path(base_dir, "fig_1toMany_forest")

dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(forest_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(loo_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(loo_dir,"tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(funnel_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(overview_dir,recursive = TRUE, showWarnings = FALSE)

#========================
# 2) Scatter plots (core 4)
#========================

make_scatter <- function(outcome_id) {
  dat_o <- dplyr::filter(dat, id.outcome == outcome_id)
  res_o <- res_all |>
    filter(id.outcome == outcome_id,
           method %in% c("Inverse variance weighted","MR Egger",
                         "Weighted median","Weighted mode","Simple mode")) |>
    arrange(method)
  if (!nrow(dat_o) || !nrow(res_o)) return(invisible(NULL))
  p <- mr_scatter_plot(res_o, dat_o)[[1]] +
    ggtitle(paste0("MR: ", label_map[[outcome_id]])) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 13, face = "bold"))
  base <- safe_filename(label_map[[outcome_id]])
  ggsave(file.path(scatter_dir, paste0(base, "_scatter.png")), p, w = 7.5, h = 6, dpi = 300)
  ggsave(file.path(scatter_dir, paste0(base, "_scatter.pdf")),  p, w = 7.5, h = 6)
  p
}
scatter_plots <- map(core_ids, make_scatter)

#========================
# 3) Forest plots (core 4)
#========================

make_forest <- function(outcome_id) {
  res_single_o <- res_single_all |> filter(id.outcome == outcome_id)
  if (!nrow(res_single_o)) return(invisible(NULL))
  p <- mr_forest_plot(res_single_o)[[1]] +
    ggtitle(label_map[[outcome_id]]) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 13, face = "bold"))
  base <- safe_filename(label_map[[outcome_id]])
  ggsave(file.path(forest_dir, paste0(base, "_forest.png")), p, w = 7.5, h = 10, dpi = 300)
  ggsave(file.path(forest_dir, paste0(base, "_forest.pdf")),  p, w = 7.5, h = 10)
  p
}
forest_plots <- map(core_ids, make_forest)

#========================
# 4) Leave-one-out (MR-Egger) + tables (core 4)
#========================

make_loo_egger <- function(outcome_id) {
  dat_o <- dat |> filter(id.outcome == outcome_id)
  if (nrow(dat_o) < 3) return(invisible(NULL))
  loo <- try(mr_leaveoneout(dat_o, method = mr_egger_regression), silent = TRUE)
  if (inherits(loo, "try-error") || !nrow(loo)) return(invisible(NULL))
  write.csv(loo, file.path(loo_dir, "tables",
                           paste0(safe_filename(label_map[[outcome_id]]), "_LOO_Egger.csv")),
            row.names = FALSE)
  p <- mr_leaveoneout_plot(loo)[[1]] +
    ggtitle(paste0("Leave-one-out (MR-Egger): ", label_map[[outcome_id]])) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 12, face = "bold"))
  base <- safe_filename(label_map[[outcome_id]])
  ggsave(file.path(loo_dir, paste0(base, "_leaveoneout_egger.png")), p, w = 7.5, h = 6.5, dpi = 300)
  ggsave(file.path(loo_dir, paste0(base, "_leaveoneout_egger.pdf")),  p, w = 7.5, h = 6.5)
  p
}
loo_plots <- map(core_ids, make_loo_egger)

#========================
# 5) Funnel plots (core 4)
#========================

make_funnel <- function(outcome_id) {
  res_single_o <- res_single_all |> filter(id.outcome == outcome_id)
  if (!nrow(res_single_o)) return(invisible(NULL))
  p <- mr_funnel_plot(res_single_o)[[1]] +
    ggtitle(paste0("Funnel: ", label_map[[outcome_id]])) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 12, face = "bold"))
  base <- safe_filename(label_map[[outcome_id]])
  ggsave(file.path(funnel_dir, paste0(base, "_funnel.png")), p, w = 7, h = 6, dpi = 300)
  ggsave(file.path(funnel_dir, paste0(base, "_funnel.pdf")),  p, w = 7, h = 6)
  p
}
funnel_plots <- map(core_ids, make_funnel)

#========================
# 6) 1-to-many forest (overview)
#========================

keep_ids <- c(core_ids, extra_ids)

res_overview <- res_all |>
  filter(id.outcome %in% keep_ids) |>
  subset_on_method() |>
  mutate(
    outcome_label = dplyr::recode(id.outcome, !!!label_map),  # force clean names
    group         = dplyr::recode(id.outcome, !!!group_map),
    OR  = exp(b), LCL = exp(b - 1.96*se), UCL = exp(b + 1.96*se),
    Effect_CI = sprintf("%.2f (%.2f–%.2f)", OR, LCL, UCL),
    P_value   = formatC(pval, format = "e", digits = 2),
    Method    = dplyr::recode(method,
                              "Inverse variance weighted"="IVW",
                              "MR Egger"="MR-Egger",
                              "Wald ratio"="Wald",
                              .default = method),
    weight_plot = 0.6 * (1/se) / max(1/se, na.rm = TRUE)
  )

res_overview$outcome_label <- str_wrap(res_overview$outcome_label, width = 40)
res_overview$group <- factor(res_overview$group, levels = group_levels)
res_overview <- res_overview |> arrange(group, desc(abs(b)))

lo <- max(0.30, min(res_overview$LCL, na.rm = TRUE) * 0.90)
up <- min(3.00,  max(res_overview$UCL, na.rm = TRUE) * 1.10)

p_overview <- forest_plot_1_to_many(
  res_overview,
  b = "b", se = "se",
  exponentiate = TRUE,
  ao_slc = FALSE,
  lo = lo, up = up,
  TraitM = "outcome_label",
  by = "group",
  trans = "log2",
  xlab = "Odds ratio per unit increase in genetic liability to endometriosis (95% CI)",
  weight = "weight_plot",
  subheading_size = 10,
  col1_title = "Outcome",
  col1_width = 4.8,
  col_text_size = 3.2,
  addcols = c("nsnp","Method","Effect_CI","P_value"),
  addcol_widths = c(0.9, 1.0, 2.1, 1.2),
  addcol_titles = c("No. SNPs","Method","Effect","P-value")
) +
  ggtitle("Mendelian randomization: Endometriosis → Pregnancy outcomes") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    strip.text.y = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 11),
    axis.text.y  = element_text(size = 9),
    plot.margin  = margin(8, 14, 8, 10)
  )

n_rows <- nrow(res_overview)
fig_h  <- max(10, min(0.35 * n_rows, 18))
fig_w  <- 11

ggsave(file.path(overview_dir, "MR_overview_1toMany.png"), p_overview, width = fig_w, height = fig_h, dpi = 300)
ggsave(file.path(overview_dir, "MR_overview_1toMany.pdf"),  p_overview, width = fig_w, height = fig_h)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(overview_dir, "MR_overview_1toMany.svg"),  p_overview, width = fig_w, height = fig_h)
}

#========================
# 7) Quick preview
#========================
print(p_overview)

message("\nSaved to:\n- ", scatter_dir,
        "\n- ", forest_dir,
        "\n- ", loo_dir,
        "\n- ", funnel_dir,
        "\n- ", overview_dir)
