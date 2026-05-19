#!/usr/bin/env Rscript
###############################################################################
# Script: 07_plots_endoMR-PREG.R
# Project: endoMR-PREG
#
# Purpose:
#   Generate all main and supplementary figures for the manuscript.
#
#   Main figures:
#     Figure 2     — Two-panel forest plot (binary OR + continuous β,
#                    all 30 outcomes, IVW primary results)
#     Figure 3     — Trios overview: maternal vs fetal vs paternal genetic
#                    effects (DONUTS estimates, conditioned)
#
#   Supplementary figures:
#     S1–S9        — Domain-specific multi-method forest plots (IVW + Egger +
#                    Weighted Median + Weighted Mode), one per domain
#                    S1 Placental disorders, S2 Bleeding & haemorrhage,
#                    S3 Pregnancy timing, S4 Labour & delivery,
#                    S5 Hypertensive disorders, S6 Fetal growth & birthweight,
#                    S7 Maternal metabolic/haematologic, S8 Maternal mental health,
#                    S9 Neonatal condition
#     S10–S15      — SNP-level leave-one-out for nominally significant outcomes:
#                    S10 Placenta praevia, S11 Premature placental separation,
#                    S12 PROM, S13 Preterm birth, S14 Very preterm birth,
#                    S15 Elective caesarean section
#
# Inputs:
#   results/ivw_results.csv
#   results/all_mr_methods.csv
#   results/trios_adj_mr_results_by_outcome_long.csv
#   results/mr_leaveoneout_snp.csv
#
# Outputs (results/plots/):
#   Figure_2_forest_plot_two_panels.png / .pdf
#   Figure_3_trios_overview.png         / .pdf
#   Supplementary_Figure_S{1-9}_*.png  / .pdf
#   Supplementary_Figure_S{10-15}_LOO_SNP_*.png / .pdf
###############################################################################

### 0) SETUP ##################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(here)
  library(patchwork)
  library(tibble)
})

results_dir <- here::here("results")
figures_dir <- file.path(results_dir, "plots")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1) OUTCOME METADATA  (30 outcomes, ordered by domain)
###############################################################################

domain_levels <- c(
  "Placental disorders",
  "Bleeding & haemorrhage",
  "Pregnancy timing",
  "Labour & delivery",
  "Hypertensive disorders",
  "Fetal growth & birthweight",
  "Maternal metabolic/haematologic",
  "Maternal mental health",
  "Neonatal condition"
)

outcome_meta <- tibble::tribble(
  ~outcome_id,                                      ~label,                               ~domain,                           ~scale,       ~order_within,

  # Placental disorders
  "finngen_R12_O15_PLAC_PRAEVIA",                   "Placenta praevia",                   "Placental disorders",             "binary",     1,
  "finngen_R12_O15_PLAC_DISORD",                    "Placental disorders (other)",         "Placental disorders",             "binary",     2,
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",              "Premature placental separation",      "Placental disorders",             "binary",     3,

  # Bleeding & haemorrhage
  "Antepartum_bleeding",                            "Antepartum bleeding",                 "Bleeding & haemorrhage",          "binary",     1,
  "Postpartum_hemorrhage",                          "Postpartum haemorrhage (any)",        "Bleeding & haemorrhage",          "binary",     2,
  "Postpartum_hemorrhage_due_to_atony",             "PPH due to uterine atony",            "Bleeding & haemorrhage",          "binary",     3,
  "Postpartum_hemorrhage_due_to_retained_placenta", "PPH due to retained placenta",        "Bleeding & haemorrhage",          "binary",     4,

  # Pregnancy timing
  "pretb_all",                                      "Preterm birth <37 weeks",             "Pregnancy timing",                "binary",     1,
  "vpretb_all",                                     "Very preterm birth <34 weeks",        "Pregnancy timing",                "binary",     2,
  "posttb_all",                                     "Post-term birth",                     "Pregnancy timing",                "binary",     3,
  "ga_all",                                         "Gestational age",                     "Pregnancy timing",                "continuous", 4,
  "rup_memb",                                       "Premature rupture of membranes",      "Pregnancy timing",                "binary",     5,

  # Labour & delivery
  "induction",                                      "Labour induction",                    "Labour & delivery",               "binary",     1,
  "el_cs",                                          "Elective caesarean section",          "Labour & delivery",               "binary",     2,
  "em_cs",                                          "Emergency caesarean section",         "Labour & delivery",               "binary",     3,

  # Hypertensive disorders
  "gh_subsamp",                                     "Gestational hypertension",            "Hypertensive disorders",          "binary",     1,
  "hdp_subsamp",                                    "Hypertensive disorders of pregnancy", "Hypertensive disorders",          "binary",     2,
  "pe_subsamp",                                     "Preeclampsia",                        "Hypertensive disorders",          "binary",     3,

  # Fetal growth & birthweight
  "sga",                                            "Small for gestational age",           "Fetal growth & birthweight",      "binary",     1,
  "lga",                                            "Large for gestational age",           "Fetal growth & birthweight",      "binary",     2,
  "lbw_all",                                        "Low birthweight <2500 g",             "Fetal growth & birthweight",      "binary",     3,
  "hbw_all",                                        "High birthweight >4000 g",            "Fetal growth & birthweight",      "binary",     4,
  "zbw_all",                                        "Birthweight z-score",                 "Fetal growth & birthweight",      "continuous", 5,

  # Maternal metabolic/haematologic
  "gdm_subsamp",                                    "Gestational diabetes",                "Maternal metabolic/haematologic", "binary",     1,
  "anaemia_preg_all",                               "Pregnancy anaemia",                   "Maternal metabolic/haematologic", "binary",     2,

  # Maternal mental health
  "depr_subsamp",                                   "Postpartum depression",               "Maternal mental health",          "binary",     1,

  # Neonatal condition
  "lowapgar1",                                      "Low Apgar score at 1 min",            "Neonatal condition",              "binary",     1,
  "lowapgar5",                                      "Low Apgar score at 5 min",            "Neonatal condition",              "binary",     2,
  "nicu",                                           "NICU admission",                      "Neonatal condition",              "binary",     3,
  "sb_subsamp",                                     "Stillbirth",                          "Neonatal condition",              "binary",     4
)

###############################################################################
# 2) DOMAIN COLOUR PALETTE
###############################################################################

domain_colours <- c(
  "Placental disorders"             = "#E69F00",
  "Bleeding & haemorrhage"          = "#D55E00",
  "Pregnancy timing"                = "#56B4E9",
  "Labour & delivery"               = "#0072B2",
  "Hypertensive disorders"          = "#CC79A7",
  "Fetal growth & birthweight"      = "#009E73",
  "Maternal metabolic/haematologic" = "#F0E442",
  "Maternal mental health"          = "#999999",
  "Neonatal condition"              = "#332288"
)

###############################################################################
# 3) LOAD IVW RESULTS
###############################################################################

ivw_path <- file.path(results_dir, "ivw_results.csv")
stopifnot(file.exists(ivw_path))
ivw_raw <- readr::read_csv(ivw_path, show_col_types = FALSE)

###############################################################################
# 4) BASE DATA PREPARATION
###############################################################################

df_base <- ivw_raw %>%
  inner_join(outcome_meta, by = c("outcome" = "outcome_id")) %>%
  mutate(
    domain = factor(domain, levels = domain_levels),

    # x-axis positions (log-scale for binary, linear for continuous)
    x_pos = b,
    x_lci = b - 1.96 * se,
    x_uci = b + 1.96 * se,

    # Effect display: OR for binary, β for continuous
    effect_est = if_else(scale == "binary", exp(b),             b),
    effect_lci = if_else(scale == "binary", exp(b - 1.96 * se), b - 1.96 * se),
    effect_uci = if_else(scale == "binary", exp(b + 1.96 * se), b + 1.96 * se),

    # Significance flag (nominal)
    significant = as.numeric(pval) < 0.05
  ) %>%
  arrange(domain, order_within) %>%
  mutate(
    est_text = if_else(scale == "binary",
                       sprintf("%.2f",           effect_est),
                       sprintf("%.3f",           effect_est)),
    ci_text  = if_else(scale == "binary",
                       sprintf("(%.2f, %.2f)",   effect_lci, effect_uci),
                       sprintf("(%.3f, %.3f)",   effect_lci, effect_uci)),
    p_text   = dplyr::case_when(
      as.numeric(pval) < 0.001 ~ formatC(as.numeric(pval), format = "e", digits = 1),
      TRUE                     ~ sprintf("%.3f", as.numeric(pval))
    ),
    q_text   = dplyr::case_when(
      as.numeric(qval) < 0.001 ~ formatC(as.numeric(qval), format = "e", digits = 1),
      TRUE                     ~ sprintf("%.3f", as.numeric(qval))
    ),
    n_text   = as.character(nsnp)
  )

stopifnot(nrow(df_base) == 30)

###############################################################################
# 5) SHARED LAYOUT CONSTANTS
###############################################################################

# Binary / OR x-axis (log scale, displayed as OR)
x_axis_breaks_bin <- log(c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0))
x_axis_labels_bin <- c("0.50", "0.75", "1.00", "1.25", "1.50", "2.00")
x_plot_min_bin    <- log(0.40)
x_plot_max_bin    <- log(2.50)

# Continuous / β x-axis (linear scale)
x_axis_breaks_cont <- c(-0.50, -0.25, 0.00, 0.25, 0.50)
x_axis_labels_cont <- c("-0.50", "-0.25", "0.00", "0.25", "0.50")
x_plot_min_cont    <- -0.70
x_plot_max_cont    <-  0.70

# Text column x positions
gap           <- 0.10
col_width     <- 0.55
x_or          <- x_plot_max_bin + gap
x_ci          <- x_or  + col_width * 1.3
x_p           <- x_ci  + col_width * 1.2
x_q           <- x_p   + col_width * 0.9
x_n           <- x_q   + col_width * 0.8
x_right_limit <- x_n   + col_width * 0.6

# Outcome label and domain label x positions
x_label        <- x_plot_min_bin - 0.05
x_domain_label <- x_plot_min_bin - 2.20

###############################################################################
# 6) HELPER: compute_domain_bounds()
###############################################################################

compute_domain_bounds <- function(df_panel) {
  df_panel %>%
    group_by(domain) %>%
    summarise(
      y_min    = min(y_pos) - 0.5,
      y_max    = max(y_pos) + 0.5,
      y_center = mean(y_pos),
      .groups = "drop"
    ) %>%
    mutate(band_fill = ifelse(row_number() %% 2 == 0, "grey93", "white"))
}

###############################################################################
# 7) HELPER: build_forest_panel()
###############################################################################

build_forest_panel <- function(df_panel, db_panel,
                               x_forest_min, x_forest_max,
                               x_breaks, x_labels,
                               x_axis_title,
                               col_header,
                               parse_col_header = FALSE,
                               parse_x_title    = FALSE) {

  y_top    <- max(df_panel$y_pos) + 2.4
  y_bottom <- min(df_panel$y_pos) - 0.5

  ggplot(df_panel, aes(y = y_pos)) +

    # Alternating domain bands
    geom_rect(
      data        = db_panel,
      inherit.aes = FALSE,
      aes(xmin = x_forest_min, xmax = x_forest_max,
          ymin = y_min, ymax = y_max, fill = band_fill),
      alpha = 0.6
    ) +
    scale_fill_identity() +

    # Null line
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = 0.5) +

    # Domain labels (left margin)
    geom_text(
      data        = db_panel,
      inherit.aes = FALSE,
      aes(x      = x_domain_label,
          y      = y_center,
          label  = as.character(domain),
          colour = as.character(domain)),
      hjust = 0, vjust = 0.5, size = 4.7, fontface = "bold"
    ) +

    # CI lines
    geom_segment(
      aes(x = x_lci, xend = x_uci, y = y_pos, yend = y_pos,
          colour = domain),
      linewidth = 0.55
    ) +

    # Points: binary = filled circle, continuous = diamond
    geom_point(
      data = filter(df_panel, scale == "binary"),
      aes(x = x_pos, y = y_pos, colour = domain),
      shape = 16, size = 2.5
    ) +
    geom_point(
      data = filter(df_panel, scale == "continuous"),
      aes(x = x_pos, y = y_pos, colour = domain),
      shape = 18, size = 3.2
    ) +

    # Colour scale
    scale_colour_manual(values = domain_colours, guide = "none") +

    # Outcome labels
    geom_text(
      aes(x = x_label, y = y_pos, label = label),
      hjust = 1, vjust = 0.5, size = 4.4
    ) +

    # Text columns
    geom_text(aes(x = x_or, y = y_pos, label = est_text),
              hjust = 0.5, vjust = 0.5, size = 4.2) +
    geom_text(aes(x = x_ci, y = y_pos, label = ci_text),
              hjust = 0.5, vjust = 0.5, size = 4.2) +
    geom_text(aes(x = x_p,  y = y_pos, label = p_text),
              hjust = 0.5, vjust = 0.5, size = 4.2) +
    geom_text(aes(x = x_q,  y = y_pos, label = q_text),
              hjust = 0.5, vjust = 0.5, size = 4.2) +
    geom_text(aes(x = x_n,  y = y_pos, label = n_text),
              hjust = 0.5, vjust = 0.5, size = 4.2) +

    # Column headers
    annotate("text", x = x_label,        y = y_top, label = "Outcome",
             hjust = 1,   vjust = 0.5, size = 4.7, fontface = "bold") +
    annotate("text", x = x_domain_label, y = y_top, label = "Domain",
             hjust = 0,   vjust = 0.5, size = 4.7, fontface = "bold") +
    annotate("text", x = x_or,           y = y_top, label = col_header,
             hjust = 0.5, vjust = 0.5, size = 4.7, fontface = "bold",
             parse = parse_col_header) +
    annotate("text", x = x_ci,           y = y_top, label = "95% CI",
             hjust = 0.5, vjust = 0.5, size = 4.7, fontface = "bold") +
    annotate("text", x = x_p,            y = y_top, label = "P value",
             hjust = 0.5, vjust = 0.5, size = 4.7, fontface = "bold") +
    annotate("text", x = x_q,            y = y_top, label = "Q value",
             hjust = 0.5, vjust = 0.5, size = 4.7, fontface = "bold") +
    annotate("text", x = x_n,            y = y_top, label = "SNPs (no.)",
             hjust = 0.5, vjust = 0.5, size = 4.7, fontface = "bold") +

    # Manual x-axis baseline
    annotate("segment",
             x = x_forest_min, xend = x_forest_max,
             y = y_bottom - 0.35, yend = y_bottom - 0.35,
             colour = "grey50", linewidth = 0.4) +

    # Manual x-axis ticks
    annotate("segment",
             x    = x_breaks, xend = x_breaks,
             y    = y_bottom - 0.20, yend = y_bottom - 0.50,
             colour = "grey50", linewidth = 0.35) +

    # Manual x-axis tick labels
    annotate("text",
             x     = x_breaks,
             y     = y_bottom - 0.75,
             label = x_labels,
             size  = 3.7, colour = "grey30", hjust = 0.5, vjust = 1) +

    # Manual x-axis title
    annotate("text",
             x     = (x_forest_min + x_forest_max) / 2,
             y     = y_bottom - 1.7,
             label = x_axis_title,
             size  = 4.7, fontface = "bold",
             parse = parse_x_title) +

    # Scales
    scale_x_continuous(
      limits = c(x_domain_label - 0.05, x_right_limit),
      breaks = NULL,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(y_bottom - 2.2, y_top + 0.5),
      expand = c(0, 0)
    ) +

    # Theme
    theme_void() +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.line.x      = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(t = 12, r = 10, b = 10, l = 10, unit = "mm")
    ) +
    labs(title = NULL, subtitle = NULL)
}

###############################################################################
# 8) FIGURE 2 — TWO-PANEL FOREST PLOT
#
#   Top panel   : 28 binary outcomes   → OR scale
#   Bottom panel:  2 continuous outcomes → β scale (separate x-axis)
###############################################################################

# Top panel: binary outcomes
df_binary <- df_base %>%
  filter(scale == "binary") %>%
  mutate(
    y_pos = rev(row_number()),
    label = factor(label, levels = rev(label))
  )

db_binary <- compute_domain_bounds(df_binary)

p_binary <- build_forest_panel(
  df_panel     = df_binary,
  db_panel     = db_binary,
  x_forest_min = x_plot_min_bin,
  x_forest_max = x_plot_max_bin,
  x_breaks     = x_axis_breaks_bin,
  x_labels     = x_axis_labels_bin,
  x_axis_title = "Odds ratio (95% CI)",
  col_header   = "OR"
)

# Bottom panel: continuous outcomes
df_cont <- df_base %>%
  filter(scale == "continuous") %>%
  mutate(
    y_pos = rev(row_number()),
    label = factor(label, levels = rev(label))
  )

db_cont <- compute_domain_bounds(df_cont)

p_cont <- build_forest_panel(
  df_panel         = df_cont,
  db_panel         = db_cont,
  x_forest_min     = x_plot_min_cont,
  x_forest_max     = x_plot_max_cont,
  x_breaks         = x_axis_breaks_cont,
  x_labels         = x_axis_labels_cont,
  x_axis_title     = "bold(beta)~bold('(95% CI)')",
  col_header       = "bold(beta)",
  parse_col_header = TRUE,
  parse_x_title    = TRUE
)

# Combine with patchwork (heights proportional to number of outcomes)
n_bin  <- nrow(df_binary)   # 28
n_cont <- nrow(df_cont)     #  2

h_bin  <- n_bin  + 4
h_cont <- n_cont + 7

p_fig2 <- p_binary / p_cont +
  patchwork::plot_layout(heights = c(h_bin, h_cont))

ggsave(file.path(figures_dir, "Figure_2_forest_plot_two_panels.png"),
       plot = p_fig2, width = 18, height = 18, dpi = 300, bg = "white")
ggsave(file.path(figures_dir, "Figure_2_forest_plot_two_panels.pdf"),
       plot = p_fig2, width = 18, height = 18, bg = "white")

message("Figure 2 saved.")

###############################################################################
# 9) FIGURE 3 — TRIOS OVERVIEW
#
#   Three components per outcome (per Caro's suggestion):
#     1. Maternal (unadj.)    — beta_mat, conventional maternal GWAS estimate
#     2. Maternal (adj. fetal)— beta_mat_donuts, adjusted for fetal genotype
#     3. Paternal (adj.)      — beta_pat_donuts, negative control
#   (Fetal (adj.) kept in Supplementary Table S4 but not in main figure)
#
#   Two panels like Figure 2:
#     Top    : 28 binary outcomes  → OR scale (log)
#     Bottom :  2 continuous outcomes → β scale (linear)
###############################################################################

trios_path <- file.path(results_dir, "trios_adj_mr_results_by_outcome_long.csv")
stopifnot(file.exists(trios_path))

trios_raw <- readr::read_csv(trios_path, show_col_types = FALSE)

# Colour palette — 3 components for the main figure
trio_colours <- c(
  "Maternal (unadj.)"     = "#636363",
  "Maternal (adj. fetal)" = "#1f78b4",
  "Paternal (adj.)"       = "#e31a1c"
)
trio_shapes <- c(
  "Maternal (unadj.)"     = 16,
  "Maternal (adj. fetal)" = 17,
  "Paternal (adj.)"       = 15
)

# Merge domain info; keep only the 3 main-figure components
trios_df <- trios_raw %>%
  left_join(
    outcome_meta %>% select(outcome_id, label, domain, order_within),
    by = "outcome_id"
  ) %>%
  filter(
    origin %in% c("Maternal (unadj.)", "Maternal (adj. fetal)", "Paternal (adj.)"),
    !is.na(domain)
  ) %>%
  mutate(
    domain = factor(domain, levels = domain_levels),
    origin = factor(origin,
                    levels = c("Maternal (unadj.)", "Maternal (adj. fetal)", "Paternal (adj.)"))
  ) %>%
  arrange(domain, order_within)

# Helper: build one trios panel (binary OR continuous)
#
# Architecture mirrors build_forest_panel: theme_void() + manual x-axis drawn
# with annotate(). The x scale spans [x_domain_label-0.05, x_forest_lim[2]+0.02],
# giving a ~35% label zone and ~65% forest zone in pixel width.
#
#   x_domain_label : x position where domain labels start (left-aligned)
#   x_label        : x position where outcome labels end  (right-aligned)
#   x_forest_lim   : c(min, max) of the visible forest area
#   x_breaks / x_labels : axis tick positions and text (in original units)
#   parse_x_title  : if TRUE, x_title is parsed as plotmath expression
build_trios_panel <- function(df_panel,
                              x_null, x_forest_lim,
                              x_breaks, x_labels,
                              x_title,
                              x_label,
                              x_domain_label,
                              is_log = TRUE,
                              parse_x_title = FALSE,
                              x_margin = 0.05) {

  # Transform data to plot coordinates
  if (is_log) {
    df_panel <- df_panel %>%
      mutate(x_est = log(OR), x_lci = log(LCL), x_uci = log(UCL))
    x_null_plot   <- log(x_null)
    x_breaks_plot <- log(x_breaks)
  } else {
    df_panel <- df_panel %>%
      mutate(x_est = OR, x_lci = LCL, x_uci = UCL)
    x_null_plot   <- x_null
    x_breaks_plot <- x_breaks
  }

  # y positions (one per outcome, top-to-bottom)
  outcome_order <- df_panel %>%
    filter(origin == "Maternal (unadj.)") %>%
    arrange(domain, order_within) %>%
    mutate(y_pos = rev(seq_len(n()))) %>%
    select(outcome_id, y_pos)

  df_panel <- df_panel %>% left_join(outcome_order, by = "outcome_id")

  # Domain bands
  band <- df_panel %>%
    filter(origin == "Maternal (unadj.)") %>%
    arrange(domain, order_within) %>%
    mutate(y_pos = rev(seq_len(n()))) %>%
    group_by(domain) %>%
    summarise(
      y_min    = min(y_pos) - 0.6,
      y_max    = max(y_pos) + 0.6,
      y_center = median(y_pos),
      .groups  = "drop"
    ) %>%
    mutate(band_fill = ifelse(row_number() %% 2 == 0, "grey93", "white"))

  y_top    <- max(df_panel$y_pos, na.rm = TRUE) + 2.4
  y_bottom <- min(df_panel$y_pos, na.rm = TRUE) - 0.5

  x_left  <- x_domain_label - x_margin
  # Proportional right pad (2% of forest width) so null fraction is the same
  # in any panel that uses the same f_label and a symmetric forest.
  x_right <- x_forest_lim[2] + 0.02 * (x_forest_lim[2] - x_forest_lim[1])

  # dodge offset per origin level
  n_orig  <- nlevels(df_panel$origin)
  dodge_w <- 0.25

  ggplot(df_panel, aes(y = y_pos)) +

    # Alternating domain bands (forest area only)
    geom_rect(
      data = band, inherit.aes = FALSE,
      aes(xmin = x_forest_lim[1], xmax = x_forest_lim[2],
          ymin = y_min, ymax = y_max, fill = band_fill),
      alpha = 0.55
    ) +
    scale_fill_identity() +

    # Vertical gridlines at breaks
    annotate("segment",
             x = x_breaks_plot, xend = x_breaks_plot,
             y = y_bottom - 0.3, yend = y_top - 0.5,
             colour = "grey90", linewidth = 0.3) +

    # Null reference line
    geom_vline(xintercept = x_null_plot, linetype = "dashed",
               colour = "grey40", linewidth = 0.45) +

    # CI bars
    geom_errorbarh(
      aes(xmin   = x_lci,
          xmax   = x_uci,
          y      = y_pos + (as.numeric(origin) - (n_orig + 1) / 2) * dodge_w,
          colour = origin),
      height = 0.15, linewidth = 0.55, lineend = "round"
    ) +

    # Estimate points
    geom_point(
      aes(x      = x_est,
          y      = y_pos + (as.numeric(origin) - (n_orig + 1) / 2) * dodge_w,
          colour = origin,
          shape  = origin),
      size = 2.8
    ) +

    # Column headers
    annotate("text", x = x_domain_label, y = y_top,
             label = "Domain", hjust = 0, vjust = 0.5,
             size = 4.6, fontface = "bold") +
    annotate("text", x = x_label, y = y_top,
             label = "Outcome", hjust = 1, vjust = 0.5,
             size = 4.6, fontface = "bold") +

    # Domain labels (left-aligned, far left)
    annotate("text",
             x      = x_domain_label,
             y      = band$y_center,
             label  = as.character(band$domain),
             colour = domain_colours[as.character(band$domain)],
             hjust  = 0, vjust = 0.5, size = 4.2, fontface = "bold") +

    # Outcome labels (right-aligned, just left of forest)
    geom_text(
      data = df_panel %>% filter(origin == "Maternal (unadj.)"),
      aes(x = x_label, y = y_pos, label = label),
      hjust = 1, vjust = 0.5, size = 4.2, colour = "grey20"
    ) +

    # Manual x-axis: baseline, ticks, tick labels, title
    annotate("segment",
             x = x_forest_lim[1], xend = x_forest_lim[2],
             y = y_bottom - 0.30, yend = y_bottom - 0.30,
             colour = "grey50", linewidth = 0.4) +
    annotate("segment",
             x = x_breaks_plot, xend = x_breaks_plot,
             y = y_bottom - 0.20, yend = y_bottom - 0.45,
             colour = "grey50", linewidth = 0.35) +
    annotate("text",
             x = x_breaks_plot, y = y_bottom - 0.75,
             label = x_labels,
             size = 3.7, colour = "grey30", hjust = 0.5, vjust = 1) +
    annotate("text",
             x     = (x_forest_lim[1] + x_forest_lim[2]) / 2,
             y     = y_bottom - 1.85,
             label = x_title,
             size  = 4.7, fontface = "bold", hjust = 0.5,
             parse = parse_x_title) +

    scale_colour_manual(values = trio_colours, name = "Genetic effect") +
    scale_shape_manual(values  = trio_shapes,  name = "Genetic effect") +

    scale_x_continuous(limits = c(x_left, x_right), breaks = NULL, expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_bottom - 2.2, y_top + 0.5), expand = c(0, 0)) +

    labs(x = NULL) +

    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title    = element_text(size = 13, face = "bold"),
      legend.text     = element_text(size = 13),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin     = margin(t = 6, r = 15, b = 4, l = 5, unit = "mm")
    )
}

# Split binary / continuous
trios_bin  <- trios_df %>% filter(scale == "binary")
trios_cont <- trios_df %>% filter(scale == "continuous")

# ─── Layout helper ────────────────────────────────────────────────────────────
# For a symmetric forest c(-A, A) and label-zone fraction f_label, the formula
# x_domain_aligned(A) gives x_domain_label so that:
#   (a) label zone = f_label * total panel width
#   (b) null (at 0) falls at the SAME pixel fraction in every panel with the
#       same f_label — regardless of A — because the null fraction simplifies to
#       (1 + 1.04*f_label) / (1 + 1.04*f_label + 1.04*(1-f_label)) = const.
#
# Derivation: x_right = A + 0.02*2A = 1.04A (2% proportional right pad)
#   label_zone / total = f_label
#   (-A - x_domain + 0.05) / (1.04A - x_domain + 0.05) = f_label
#   → x_domain = -A * (1 + 1.04*f_label)/(1-f_label) + 0.05
f_label_trios <- 0.44
x_domain_aligned <- function(A, f = f_label_trios) {
  -A * (1 + 1.04 * f) / (1 - f) + 0.05
}

# Dynamic forest half-widths from data
# Binary : A = max |log(CI)| rounded up to 0.1, but at least log(4)≈1.386 so
#          standard breaks 0.25/0.50/1.00/2.00/4.00 always fit inside.
# Continuous: A = max |CI| rounded up to 0.05 with 15% padding.
A_bin_fig3  <- max(
  ceiling(max(abs(log(c(trios_bin$LCL,  trios_bin$UCL))),  na.rm = TRUE) / 0.1) * 0.1,
  log(4)   # ensures breaks at 0.25 / 4.00 are within the forest
)
A_cont_fig3 <- ceiling(
  max(abs(c(trios_cont$LCL, trios_cont$UCL)), na.rm = TRUE) * 1.15 / 0.05
) * 0.05

trios_x_domain_bin  <- x_domain_aligned(A_bin_fig3)
trios_x_label_bin   <- -A_bin_fig3  - 0.08   # right-aligned, just left of forest
trios_x_domain_cont <- x_domain_aligned(A_cont_fig3)
trios_x_label_cont  <- -A_cont_fig3 - 0.04   # right-aligned, just left of forest

# x_breaks for binary chosen so all fall within c(-A_bin_fig3, A_bin_fig3):
#   log(0.25)= -1.39, log(4.00)= 1.39 — guaranteed by max(., log(4)) above.
# x_breaks for continuous chosen to fit within c(-A_cont_fig3, A_cont_fig3).
x_breaks_cont_fig3 <- sort(unique(c(
  seq(0, -A_cont_fig3, by = -0.25), seq(0, A_cont_fig3, by = 0.25)
)))
x_breaks_cont_fig3 <- x_breaks_cont_fig3[
  x_breaks_cont_fig3 >= -A_cont_fig3 & x_breaks_cont_fig3 <= A_cont_fig3
]
x_labels_cont_fig3 <- sprintf("%.2f", x_breaks_cont_fig3)

p_t_bin <- build_trios_panel(
  trios_bin,
  x_null         = 1,
  x_forest_lim   = c(-A_bin_fig3,  A_bin_fig3),
  x_breaks       = c(0.25, 0.50, 1.00, 2.00, 4.00),
  x_labels       = c("0.25", "0.50", "1.00", "2.00", "4.00"),
  x_title        = "Odds ratio (95% CI, log scale)",
  x_label        = trios_x_label_bin,
  x_domain_label = trios_x_domain_bin,
  is_log         = TRUE
)

p_t_cont <- build_trios_panel(
  trios_cont,
  x_null         = 0,
  x_forest_lim   = c(-A_cont_fig3, A_cont_fig3),
  x_breaks       = x_breaks_cont_fig3,
  x_labels       = x_labels_cont_fig3,
  x_title        = "beta~'(95% CI)'",
  x_label        = trios_x_label_cont,
  x_domain_label = trios_x_domain_cont,
  is_log         = FALSE,
  parse_x_title  = TRUE,
  x_margin       = 0.05 * A_cont_fig3 / A_bin_fig3
)

n_bin_t  <- length(unique(trios_bin$outcome_id))
n_cont_t <- length(unique(trios_cont$outcome_id))
h_bin_t  <- n_bin_t  + 4
h_cont_t <- n_cont_t + 7

# guides = "collect" pools the legend from both panels into one shared legend,
# ensuring both panels have the same pixel width (no legend distortion in cont panel).
p_fig3 <- (p_t_bin / p_t_cont) +
  patchwork::plot_layout(heights = c(h_bin_t, h_cont_t), guides = "collect") &
  theme(legend.position = "bottom")

fig3_h <- max(10, n_bin_t * 0.42 + n_cont_t * 0.9 + 4)

ggsave(file.path(figures_dir, "Figure_3_trios_overview.png"),
       plot = p_fig3, width = 18, height = fig3_h, dpi = 300, bg = "white")
ggsave(file.path(figures_dir, "Figure_3_trios_overview.pdf"),
       plot = p_fig3, width = 18, height = fig3_h, bg = "white")

message("Figure 3 saved.")
###############################################################################
# 10) SUPPLEMENTARY FIGURES S1–S9 — DOMAIN-SPECIFIC MULTI-METHOD FOREST PLOTS
#
#   One figure per domain (9 domains = S1–S9). Shows all outcomes within the
#   domain with four MR methods: IVW (primary), Egger, Weighted Median, Weighted
#   Mode — colour- and shape-coded.  Same domain colour palette as Figure 2.
#   Domains with mixed binary/continuous outcomes use the β scale throughout.
###############################################################################

# Load all-methods results
all_methods_path <- file.path(results_dir, "all_mr_methods.csv")
stopifnot(file.exists(all_methods_path))
all_raw <- readr::read_csv(all_methods_path, show_col_types = FALSE)

# Standardise method names
method_map <- c(
  "Inverse variance weighted" = "IVW",
  "MR Egger"                 = "Egger",
  "Weighted median"          = "Weighted median",
  "Weighted mode"            = "Weighted mode"
)

method_colours <- c(
  "IVW"             = "#d62728",
  "Egger"           = "#1f77b4",
  "Weighted median" = "#2ca02c",
  "Weighted mode"   = "#9467bd"
)
method_shapes <- c(
  "IVW"             = 16,
  "Egger"           = 17,
  "Weighted median" = 15,
  "Weighted mode"   = 18
)

all_df <- all_raw %>%
  inner_join(outcome_meta, by = c("outcome" = "outcome_id")) %>%
  filter(method %in% names(method_map)) %>%
  mutate(
    method_short = dplyr::recode(method, !!!method_map),
    method_short = factor(method_short,
                          levels = c("IVW", "Egger", "Weighted median", "Weighted mode")),
    domain       = factor(domain, levels = domain_levels),
    # For binary: x = OR (exp(b)); for continuous: x = b
    x_est = if_else(scale == "binary", exp(b),             b),
    x_lci = if_else(scale == "binary", exp(b - 1.96 * se), b - 1.96 * se),
    x_uci = if_else(scale == "binary", exp(b + 1.96 * se), b + 1.96 * se),
    # Domains with mixed scales → flag to use b for all
    has_continuous = scale == "continuous"
  )

# Helper: one domain supplementary figure
build_domain_figure <- function(df_domain, dom_name, dom_colour) {

  df <- df_domain %>% arrange(order_within)
  is_mixed   <- any(df$scale == "continuous") && any(df$scale == "binary")
  is_cont    <- all(df$scale == "continuous")

  # Mixed domain: two panels (binary → OR, continuous → β)
  if (is_mixed) {
    make_panel <- function(dat, x_null, x_lab) {
      o_order <- dat %>% filter(method_short == "IVW") %>%
        arrange(order_within) %>% pull(label)
      dat <- dat %>% mutate(label = factor(label, levels = rev(o_order)))
      band <- dat %>% filter(method_short == "IVW") %>% arrange(order_within) %>%
        mutate(band_fill = ifelse(seq_len(n()) %% 2 == 0, "grey93", "white")) %>%
        select(label, band_fill)
      ggplot(dat, aes(x = x_est, y = label,
                      colour = method_short, shape = method_short)) +
        geom_rect(data = band, inherit.aes = FALSE,
                  aes(ymin = as.numeric(label) - 0.5,
                      ymax = as.numeric(label) + 0.5,
                      xmin = -Inf, xmax = Inf, fill = band_fill),
                  alpha = 0.5) +
        scale_fill_identity() +
        geom_vline(xintercept = x_null, linetype = "dashed",
                   colour = "grey40", linewidth = 0.5) +
        geom_errorbarh(aes(xmin = x_lci, xmax = x_uci),
                       height = 0.18, linewidth = 0.55,
                       position = position_dodge(width = 0.6)) +
        geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
        scale_colour_manual(values = method_colours, name = "MR method") +
        scale_shape_manual(values  = method_shapes,  name = "MR method") +
        scale_x_continuous(expand = expansion(mult = 0.08)) +
        labs(x = x_lab, y = NULL) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text.y        = element_text(size = 10, colour = "grey20"),
          axis.text.x        = element_text(size = 10),
          axis.title.x       = element_text(size = 11, face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.35),
          legend.position    = "bottom",
          legend.title       = element_text(size = 10, face = "bold"),
          legend.text        = element_text(size = 10),
          plot.background    = element_rect(fill = "white", colour = NA),
          panel.background   = element_rect(fill = "white", colour = NA),
          plot.margin        = margin(t = 6, r = 10, b = 6, l = 10, unit = "mm")
        )
    }
    df_bin  <- df %>% filter(scale == "binary")
    df_cont <- df %>% filter(scale == "continuous") %>%
      mutate(x_est = b, x_lci = b - 1.96 * se, x_uci = b + 1.96 * se)
    n_bin  <- length(unique(df_bin$label))
    n_cont <- length(unique(df_cont$label))
    p_bin  <- make_panel(df_bin,  x_null = 1, x_lab = "Odds ratio (95% CI)")
    p_cont <- make_panel(df_cont, x_null = 0, x_lab = expression(beta ~ "(95% CI)"))
    return(
      (p_bin + theme(legend.position = "none")) /
      (p_cont + theme(legend.position = "bottom")) +
      patchwork::plot_layout(heights = c(n_bin, n_cont)) +
      patchwork::plot_annotation(
        caption = "IVW = primary analysis. Egger, weighted median, weighted mode = sensitivity analyses."
      )
    )
  }

  # All-continuous domain: β scale throughout
  if (is_cont) {
    df <- df %>% mutate(x_est = b, x_lci = b - 1.96 * se, x_uci = b + 1.96 * se)
    x_null  <- 0
    x_lab   <- expression(beta ~ "(95% CI)")
    x_note  <- ""
  } else {
    x_null <- 1
    x_lab  <- "Odds ratio (95% CI)"
    x_note <- ""
  }

  # Outcome label ordered top-to-bottom
  outcome_order <- df %>%
    filter(method_short == "IVW") %>%
    arrange(order_within) %>%
    pull(label)
  df <- df %>%
    mutate(label = factor(label, levels = rev(outcome_order)))

  # Alternating row bands by outcome
  band_df <- df %>%
    filter(method_short == "IVW") %>%
    arrange(order_within) %>%
    mutate(
      y_num    = rev(seq_len(n())),
      band_fill = ifelse(seq_len(n()) %% 2 == 0, "grey93", "white")
    ) %>%
    select(label, y_num, band_fill)

  ggplot(df, aes(x = x_est, y = label,
                 colour = method_short, shape = method_short)) +

    # Band stripes
    geom_rect(
      data        = band_df,
      inherit.aes = FALSE,
      aes(ymin = as.numeric(label) - 0.5,
          ymax = as.numeric(label) + 0.5,
          xmin = -Inf, xmax = Inf,
          fill = band_fill),
      alpha = 0.5
    ) +
    scale_fill_identity() +

    # Null line
    geom_vline(xintercept = x_null, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +

    # CI + points (dodged)
    geom_errorbarh(
      aes(xmin = x_lci, xmax = x_uci),
      height   = 0.18,
      linewidth = 0.55,
      position = position_dodge(width = 0.6)
    ) +
    geom_point(
      size     = 2.4,
      position = position_dodge(width = 0.6)
    ) +

    scale_colour_manual(values = method_colours, name = "MR method") +
    scale_shape_manual(values  = method_shapes,  name = "MR method") +

    scale_x_continuous(expand = expansion(mult = 0.08)) +

    labs(
      x       = x_lab,
      y       = NULL,
      caption = paste0(
        "IVW = primary analysis. Egger, weighted median, weighted mode = sensitivity analyses.",
        x_note
      )
    ) +

    theme_minimal(base_size = 11) +
    theme(
      axis.text.y        = element_text(size = 10, colour = "grey20"),
      axis.text.x        = element_text(size = 10),
      axis.title.x       = element_text(size = 11, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.35),
      legend.position    = "bottom",
      legend.title       = element_text(size = 10, face = "bold"),
      legend.text        = element_text(size = 10),
      plot.caption       = element_text(size = 8, colour = "grey45", hjust = 0),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA),
      plot.margin        = margin(t = 6, r = 10, b = 6, l = 10, unit = "mm")
    )
}

# Supplementary figure order (S1–S9) — independent of domain_levels (Figure 2 order)
supp_domain_order <- c(
  "Placental disorders",             # S1
  "Bleeding & haemorrhage",          # S2
  "Pregnancy timing",                # S3
  "Labour & delivery",               # S4
  "Hypertensive disorders",          # S5
  "Fetal growth & birthweight",      # S6
  "Maternal metabolic/haematologic", # S7
  "Maternal mental health",          # S8
  "Neonatal condition"               # S9
)

domain_fig_map <- tibble::tibble(
  domain     = supp_domain_order,
  supp_num   = paste0("S", seq_along(supp_domain_order)),
  file_stem  = paste0("Supplementary_Figure_", supp_num, "_",
                      gsub("[/ &]", "_", supp_domain_order))
)

for (i in seq_len(nrow(domain_fig_map))) {
  dom      <- domain_fig_map$domain[i]
  stem     <- domain_fig_map$file_stem[i]
  dom_col  <- domain_colours[dom]
  df_dom   <- all_df %>% filter(domain == dom)

  if (nrow(df_dom) == 0) {
    message("No data for domain: ", dom, " — skipping.")
    next
  }

  n_out <- length(unique(df_dom$label))
  fig_h <- max(3.5, n_out * 1.1 + 2.5)

  p_dom <- build_domain_figure(df_dom, dom, dom_col)

  ggsave(file.path(figures_dir, paste0(stem, ".png")),
         plot = p_dom, width = 12, height = fig_h, dpi = 300, bg = "white")
  ggsave(file.path(figures_dir, paste0(stem, ".pdf")),
         plot = p_dom, width = 12, height = fig_h, bg = "white")

  message("Supplementary Figure ", domain_fig_map$supp_num[i], " (", dom, ") saved.")
}

###############################################################################
# 11) SUPPLEMENTARY FIGURES S10–S15 — SNP-LEVEL LEAVE-ONE-OUT
#     FOR NOMINALLY SIGNIFICANT OUTCOMES
#
#   One figure per outcome. Outcomes:
#     S10 — Placenta praevia              (finngen_R12_O15_PLAC_PRAEVIA)
#     S11 — Premature placental separation (finngen_R12_O15_PLAC_PREMAT_SEPAR)
#     S12 — Premature rupture of membranes (rup_memb)
#     S13 — Preterm birth <37 weeks        (pretb_all)
#     S14 — Very preterm birth <34 weeks   (vpretb_all)
#     S15 — Elective caesarean section     (el_cs)
###############################################################################

# Script 05 appends a timestamp: mr_leaveoneout_snp_YYYYMMDD-HHMMSS.csv
# Pick the most recently modified matching file
loo_candidates <- list.files(results_dir,
                             pattern = "^mr_leaveoneout_snp.*\\.csv$",
                             full.names = TRUE)
if (length(loo_candidates) == 0)
  stop("No mr_leaveoneout_snp*.csv found in results/. Run script 05 first.")
loo_snp_path <- loo_candidates[which.max(file.mtime(loo_candidates))]
message("Using LOO SNP file: ", basename(loo_snp_path))
loo_raw <- readr::read_csv(loo_snp_path, show_col_types = FALSE)

loo_outcomes <- tibble::tibble(
  outcome_id = c(
    "finngen_R12_O15_PLAC_PRAEVIA",
    "finngen_R12_O15_PLAC_PREMAT_SEPAR",
    "rup_memb",
    "pretb_all",
    "vpretb_all",
    "el_cs"
  ),
  supp_num = paste0("S", 10:15)
)

for (i in seq_len(nrow(loo_outcomes))) {
  oid      <- loo_outcomes$outcome_id[i]
  snum     <- loo_outcomes$supp_num[i]

  df_i <- loo_raw %>% filter(outcome == oid)
  if (nrow(df_i) == 0) { message(snum, ": no LOO data for ", oid, " — skipping."); next }

  lbl      <- outcome_meta$label[outcome_meta$outcome_id == oid]
  n_snps_i <- sum(df_i$SNP != "All")
  fig_h    <- max(5, n_snps_i * 0.22 + 2.5)
  stem_i   <- paste0("Supplementary_Figure_", snum, "_LOO_SNP_",
                     gsub("[^a-zA-Z0-9]", "_", lbl))

  p_loo_list <- TwoSampleMR::mr_leaveoneout_plot(df_i)
  p_loo      <- p_loo_list[[1]]

  ggsave(file.path(figures_dir, paste0(stem_i, ".png")),
         plot = p_loo, width = 10, height = fig_h, dpi = 300, bg = "white")
  ggsave(file.path(figures_dir, paste0(stem_i, ".pdf")),
         plot = p_loo, width = 10, height = fig_h, bg = "white")

  message("Supplementary Figure ", snum, " (LOO — ", lbl, ") saved.")
}

message("All figures written to: ", figures_dir)
