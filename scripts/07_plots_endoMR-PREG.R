#!/usr/bin/env Rscript
###############################################################################
# Script: 07_plots_endoMR-PREG.R
# Project: endoMR-PREG
#
# Purpose:
#   Generate figures for the MR manuscript:
#     - Figure 1: BP-style flowchart of the MR pipeline
#     - Figure 2: Multi-method forest plot (top FDR-significant outcomes)
#     - Figure 3: Multi-method forest plot (all outcomes)
#     - Figure 4: Maternal vs fetal vs paternal genetic effects
#
# Inputs:
#   - results/ivw_results.csv
#   - results/all_mr_methods.csv
#   - results/harmonised_pregnancy_outcomes.rds    (optional, for other plots)
#   - results/trios_mr_results_comparison.csv      (for trios forest plot)
#
# Outputs (all under results/plots/):
#   - Fig1_flowchart_BPstyle_endoMR-PREG_improved.png
#   - forest_endoMR-PREG_multiMethod_topFDR.png
#   - forest_endoMR-PREG_multiMethod_allOutcomes.png
#   - Figure_trios_maternal_fetal_paternal.png / .pdf
#
# Author: Jonas Vibert
###############################################################################

### 1) SETUP ###################################################################

required_pkgs <- c(
  "TwoSampleMR",
  "dplyr",
  "ggplot2",
  "here",
  "readr",
  "stringr",
  "purrr"
)

safe_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(required_pkgs, safe_install))

# Paths
results_dir <- here::here("results")
plots_dir   <- file.path(results_dir, "plots")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,   showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

read_csv_safe <- function(path) {
  if (!file.exists(path)) {
    message("File not found: ", path)
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

message("=== 07_plots_endoMR-PREG.R ===")

### 2) LOAD MR RESULTS + HARMONISED DATA #######################################

ivw_results_path <- file.path(results_dir, "ivw_results.csv")
all_results_path <- file.path(results_dir, "all_mr_methods.csv")
harm_path        <- file.path(results_dir, "harmonised_pregnancy_outcomes.rds")
trios_path       <- file.path(results_dir, "trios_mr_results_comparison.csv")

ivw_results <- read_csv_safe(ivw_results_path)
all_results <- read_csv_safe(all_results_path)

if (is.null(ivw_results) || is.null(all_results)) {
  stop("Core MR results not found. Run 04_main_analyses_endoMR-PREG.R first.")
}

# Ensure outcome_full present (export_csv should have added it)
if (!"outcome_full" %in% names(ivw_results)) {
  ivw_results$outcome_full <- ivw_results$outcome
}
if (!"outcome_full" %in% names(all_results)) {
  all_results$outcome_full <- all_results$outcome
}

### Harmonised data (optional; used only for scatter / funnel / LOO plots) #####

if (file.exists(harm_path)) {
  harm_list <- readRDS(harm_path)
  message("Loaded harmonised datasets (n = ", length(harm_list), ")")
} else {
  message("Harmonised RDS not found at: ", harm_path)
  message("Scatter / leave-one-SNP / funnel plots will be skipped.")
  harm_list <- list()   # Empty list instead of NULL → avoids errors later
}

message("Loaded IVW results: ", nrow(ivw_results))
message("Loaded all-methods results: ", nrow(all_results))

### Outcomes labels
outcome_labels <- c( Antepartum_bleeding = "Antepartum bleeding", Postpartum_hemorrhage = "Postpartum hemorrhage", Postpartum_hemorrhage_due_to_atony = "PPH due to atony", Postpartum_hemorrhage_due_to_retained_placenta = "PPH due to retained placenta", finngen_R12_O15_PLAC_PRAEVIA = "Placenta praevia", finngen_R12_O15_PLAC_DISORD = "Placental disorders", finngen_R12_O15_PLAC_PREMAT_SEPAR = "Premature placental separation", rup_memb = "Premature rupture of membranes", pretb_all = "Preterm birth (all)", vpretb_all = "Very preterm birth", ga_all = "Gestational age (all)", ga_subsamp = "Gestational age (subsample)", sga = "Small for gestational age", lbw_all = "Low birthweight", hbw_all = "High birthweight", lga = "Large for gestational age", zbw_all = "Z-score birthweight", lowapgar1 = "Low Apgar score at 1 min", lowapgar5 = "Low Apgar score at 5 min", nicu = "NICU admission", anaemia_preg_all = "Pregnancy anemia", gdm_subsamp = "Gestational diabetes", gh_subsamp = "Gestational hypertension", hdp_subsamp = "Hypertensive disorders of pregnancy", pe_subsamp = "Preeclampsia", induction = "Labour induction", posttb_all = "Post-term birth" )

###########################
### FIGURE 1.FLOWCHART  ###
###########################

message("Figure 1: Study flowchart created on draw.io")


#################################################################################################################################################################
### FIGURE 2. Inverse-variance weighted Mendelian randomization estimates for genetically predicted effects of endometriosis liability on pregnancy outcomes  ###
#################################################################################################################################################################
message("Creating Figure : forest plot IVW")

# Prepare data
df_ivw_layout <- ivw_filtered %>%
  filter(method == "Inverse variance weighted") %>%
  # garder uniquement les outcomes d'intérêt
  filter(outcome %in% names(outcome_labels)) %>%
  mutate(
    OR       = exp(b),
    OR_lower = exp(b - 1.96 * se),
    OR_upper = exp(b + 1.96 * se),
    
    # Transform to log2 space for symmetric visual scaling
    OR_t       = log2(OR),
    OR_lower_t = log2(OR_lower),
    OR_upper_t = log2(OR_upper),
    
    # FDR q-value
    qval   = p.adjust(pval, method = "fdr"),
    
    OR_text  = sprintf("%.2f", OR),
    CI_text  = sprintf("(%.2f-%.2f)", OR_lower, OR_upper),
    P_text   = ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval)),
    Q_text   = ifelse(qval < 0.001, "<0.001", sprintf("%.3f", qval)),
    SNP_text = paste0("n=", nsnp),
    
    color_group = case_when(
      pval < 0.001 ~ "p < 0.001",
      pval < 0.01  ~ "p < 0.01",
      pval < 0.05  ~ "p < 0.05",
      TRUE         ~ "p ≥ 0.05"
    )
  ) %>%
  arrange(desc(OR)) %>%
  mutate(y_pos = n() - row_number() + 1) %>%
  # appliquer les labels jolis
  mutate(outcome_clean = dplyr::recode(outcome, !!!outcome_labels))

if (nrow(df_ivw_layout) == 0) {
  message("No data available for IVW forest plot")
} else {
  
  forest_data <- df_ivw_layout
  
  # Dynamic height
  n_outcomes <- nrow(forest_data)
  base_height <- 8
  height_adjustment <- max(0, (n_outcomes - 20) * 0.3)
  final_height <- base_height + height_adjustment
  
  # X limits
  x_min <- 0.5
  x_max <- 2.5
  
  # Transform to log2
  x_min_t <- log2(x_min)
  x_max_t <- log2(x_max)
  
  text_start_t   <- x_max_t + 0.05
  text_spacing_t <- 0.25
  
  or_pos_t  <- text_start_t
  ci_pos_t  <- text_start_t + text_spacing_t
  p_pos_t   <- text_start_t + 2 * text_spacing_t
  q_pos_t   <- text_start_t + 3 * text_spacing_t
  snp_pos_t <- text_start_t + 4 * text_spacing_t
  
  y_bottom <- 0.5
  y_top    <- max(forest_data$y_pos) + 0.5
  y_header <- max(forest_data$y_pos) + 1
  
  panel_df <- data.frame(
    xmin = x_min_t, xmax = x_max_t,
    ymin = y_bottom, ymax = y_top
  )
  
  ref1_df <- data.frame(
    x = 0, xend = 0,
    y = y_bottom, yend = y_top
  )
  
  axis_base_df <- data.frame(
    x = x_min_t, xend = x_max_t,
    y = 0.3, yend = 0.3
  )
  
  tick_vals <- c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
  tick_df <- data.frame(x = log2(tick_vals))
  tick_lab_df <- transform(tick_df, lab = sprintf("%g", 2^x))
  
  header_df <- data.frame(
    x   = c(
      x_min_t - 0.15,
      or_pos_t,
      ci_pos_t,
      p_pos_t,
      q_pos_t,
      snp_pos_t
    ),
    lab = c(
      "Outcome",
      "OR",
      "95% CI",
      "P value",
      "Q value",
      "SNPs"
    )
  )
  
  x_limit_right <- snp_pos_t + 0.30
  
  # ---- Plot -------------------------------------------------------------------
  
  p <- ggplot(forest_data, aes(y = y_pos)) +
    geom_rect(
      data = panel_df, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "gray98", color = "gray85", alpha = 0.3
    ) +
    
    geom_segment(
      data = ref1_df, inherit.aes = FALSE,
      aes(x = x, xend = xend, y = y, yend = yend),
      linetype = "dashed", color = "gray60", linewidth = 0.5
    ) +
    
    geom_pointrange(
      aes(x = OR_t, xmin = OR_lower_t, xmax = OR_upper_t, color = color_group),
      size = 0.8, linewidth = 0.6
    ) +
    
    geom_segment(
      data = axis_base_df, inherit.aes = FALSE,
      aes(x = x, xend = xend, y = y, yend = yend),
      color = "black", linewidth = 0.5
    ) +
    
    geom_segment(
      data = transform(tick_df, xend = x, y = 0.25, yend = 0.35),
      inherit.aes = FALSE,
      aes(x = x, xend = xend, y = y, yend = yend),
      color = "black", linewidth = 0.4
    ) +
    
    geom_text(
      data = transform(tick_lab_df, y = 0.1),
      inherit.aes = FALSE,
      aes(x = x, y = y, label = lab),
      size = 3.2, hjust = 0.5, vjust = 0.5
    ) +
    
    geom_text(
      aes(x = x_min_t - 0.15, label = outcome_clean),
      hjust = 1, vjust = 0.5, size = 3.5
    ) +
    
    geom_text(aes(x = or_pos_t,  label = OR_text),  hjust = 0, vjust = 0.5, size = 3.2) +
    geom_text(aes(x = ci_pos_t,  label = CI_text),  hjust = 0, vjust = 0.5, size = 3.2) +
    geom_text(aes(x = p_pos_t,   label = P_text),   hjust = 0, vjust = 0.5, size = 3.2) +
    geom_text(aes(x = q_pos_t,   label = Q_text),   hjust = 0, vjust = 0.5, size = 3.2) +
    geom_text(aes(x = snp_pos_t, label = SNP_text), hjust = 0, vjust = 0.5, size = 3.2) +
    
    geom_text(
      data = transform(header_df, y = y_header),
      inherit.aes = FALSE,
      aes(x = x, y = y, label = lab),
      hjust = c(1, 0, 0, 0, 0, 0),
      vjust = 0.5, size = 3.8, fontface = "bold"
    ) +
    
    annotate(
      "text",
      x = (x_min_t + x_max_t) / 2, y = -0.3,
      label = "Odds Ratio", size = 3.8, fontface = "bold"
    ) +
    
    scale_x_continuous(
      limits = c(x_min_t - 0.8, x_limit_right),
      expand = c(0, 0),
      breaks = log2(tick_vals),
      labels = sprintf("%g", tick_vals)
    ) +
    scale_y_continuous(limits = c(-0.6, y_header + 1), expand = c(0, 0)) +
    scale_color_manual(
      values = c(
        "p < 0.001" = "#d73027",
        "p < 0.01"  = "#fc8d59",
        "p < 0.05"  = "#fee08b",
        "p ≥ 0.05"  = "#999999"
      ),
      name = "Significance"
    ) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position  = "bottom",
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 9),
      legend.margin    = margin(t = 10),
      plot.margin      = margin(20, 20, 20, 20),
      plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle    = element_text(hjust = 0.5, size = 11),
      plot.caption = element_text(hjust = 0)
    ) +
    labs(
      title = "Figure 2. Inverse-variance weighted Mendelian randomization estimates for genetically predicted effects of endometriosis liability on pregnancy outcomes",
      subtitle = "Two-sample Mendelian randomization (Inverse variance weighted method)",
      caption = paste(
        "Odds ratios (OR) with 95% confidence intervals were estimated using the IVW method\n",
        "and displayed on a log₂ scale for symmetric interpretation around the null (OR = 1).\n",
        "Point colours reflect nominal P-value thresholds (red: P < 0.001; orange: P < 0.01;\n",
        "yellow: P < 0.05; grey: P ≥ 0.05).\n",
        "Q values correspond to false discovery rate (FDR)–adjusted P values across all outcomes.\n",
        "Only the association shown in red reached FDR significance (Q < 0.05).\n",
        "The number of SNPs included as instruments for each outcome is provided in the rightmost column.",
        sep = ""
      )
    )
  
  ggsave(
    file.path(plots_dir, "Fig2forest_endoMR-PREG.png"),
    plot = p, width = 16, height = final_height,
    dpi = 300, bg = "white"
  )
  
  message("Saved: ", file.path(plots_dir, "Fig2forest_endoMR-PREG.png"))
}


###########################################
### FIGURE 3. MULTI-METHOD FOREST PLOTS ###
###########################################

message("Creating multi-method forest plots...")

method_levels <- c(
  "MR Egger",
  "Weighted median",
  "Inverse variance weighted",
  "Weighted mode"
)

# 0) Ensure filtered objects exist and carry a clean label
ivw_filtered <- ivw_results %>%
  dplyr::filter(method == "Inverse variance weighted") %>%
  dplyr::mutate(
    outcome_full = ifelse(is.na(outcome_full), outcome, outcome_full)
  )

all_filtered <- all_results %>%
  dplyr::filter(method %in% method_levels) %>%
  dplyr::mutate(
    outcome_full = ifelse(is.na(outcome_full), outcome, outcome_full)
  )

# --------------------------------------------------------------------------- #
# Helper: build data for multi-method forest plot
#   - Selection based on IVW P-value threshold (p_threshold)
#   - q-values kept for display and FDR colouring
# --------------------------------------------------------------------------- #

build_multi_method_df <- function(ivw_df,
                                  all_df,
                                  method_levels,
                                  top_n       = NULL,
                                  p_threshold = NULL) {
  # Require qval for FDR information
  if (!"qval" %in% names(ivw_df)) {
    stop("Column 'qval' not found in IVW results. FDR correction must be present.")
  }
  
  ivw_df <- ivw_df %>%
    dplyr::mutate(
      q_num = as.numeric(qval)
    )
  
  # --- Select focus outcomes using P-value threshold (IVW only) -------------
  if (!is.null(p_threshold)) {
    ivw_sig <- ivw_df %>%
      dplyr::filter(!is.na(pval), pval < p_threshold) %>%
      dplyr::arrange(pval)
    
    if (!is.null(top_n)) {
      ivw_sig <- dplyr::slice_head(ivw_sig, n = top_n)
    }
    
    focus_outcomes <- ivw_sig$outcome
  } else {
    # If no p-threshold is provided, keep all outcomes
    focus_outcomes <- unique(ivw_df$outcome)
  }
  
  focus_outcomes <- unique(focus_outcomes)
  if (!length(focus_outcomes)) {
    return(ivw_df[0, , drop = FALSE])
  }
  
  # --- Build multi-method dataset -------------------------------------------
  df_multi <- all_df %>%
    dplyr::filter(
      method %in% method_levels,
      outcome %in% focus_outcomes
    ) %>%
    dplyr::mutate(
      OR       = exp(as.numeric(b)),
      OR_lower = exp(as.numeric(b) - 1.96 * as.numeric(se)),
      OR_upper = exp(as.numeric(b) + 1.96 * as.numeric(se)),
      OR_text  = sprintf("%.2f", OR),
      CI_text  = sprintf("(%.2f–%.2f)", OR_lower, OR_upper),
      q_num    = as.numeric(qval),
      Q_text   = dplyr::case_when(
        is.na(q_num)           ~ "NA",
        q_num < 0.001          ~ "<0.001",
        TRUE                   ~ sprintf("%.3f", q_num)
      ),
      SNP_text = paste0("n=", nsnp),
      method   = factor(method, levels = method_levels),
      outcome_clean = ifelse(is.na(outcome_full), outcome, outcome_full),
      signif_group = dplyr::case_when(
        is.na(q_num)      ~ "q ≥ 0.05",
        q_num < 0.01      ~ "q < 0.01",
        q_num < 0.05      ~ "q < 0.05",
        TRUE              ~ "q ≥ 0.05"
      )
    )
  
  df_multi
}

# --------------------------------------------------------------------------- #
# Helper: draw multi-method forest plot
# --------------------------------------------------------------------------- #

make_multi_forest_plot <- function(df_multi,
                                   file_stub,
                                   subtitle_text) {
  if (!nrow(df_multi)) {
    warning("No data for multi-method forest plot: ", file_stub)
    return(invisible(NULL))
  }
  
  # Order outcomes alphabetically by label
  outcome_levels <- df_multi %>%
    dplyr::distinct(outcome_clean) %>%
    dplyr::arrange(outcome_clean) %>%
    dplyr::pull(outcome_clean)
  
  df_multi$outcome_clean <- factor(df_multi$outcome_clean, levels = outcome_levels)
  
  n_methods <- length(levels(df_multi$method))
  gap       <- 3
  
  base_pos <- data.frame(
    outcome_clean = outcome_levels,
    base_y        = cumsum(c(0, rep(n_methods + gap, length(outcome_levels) - 1)))
  )
  
  df_multi <- df_multi %>%
    dplyr::arrange(outcome_clean, method) %>%
    dplyr::group_by(outcome_clean) %>%
    dplyr::mutate(order_in_outcome = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(base_pos, by = "outcome_clean") %>%
    dplyr::mutate(
      y_pos        = base_y + (n_methods - order_in_outcome + 1),
      outcome_label = ifelse(order_in_outcome == 1,
                             as.character(outcome_clean), "")
    )
  
  # Visual scaling: compress OR-range on a log10 scale
  x_show_min    <- 0.1
  x_show_max    <- 6
  shrink_factor <- 0.25
  
  shrink_f <- function(x) 10^(log10(x) * shrink_factor)
  
  map_shrink <- function(x, band_left, band_right) {
    z  <- shrink_f(x)
    z0 <- shrink_f(x_show_min)
    z1 <- shrink_f(x_show_max)
    band_left + (z - z0) / (z1 - z0) * (band_right - band_left)
  }
  
  # Layout fractions
  left_frac         <- 0.25
  gap_label_to_band <- 0.06
  band_frac         <- 0.30
  gap_frac          <- 0.03
  col_step          <- 0.09
  col_pad           <- 0.010
  
  band_left   <- left_frac + gap_label_to_band
  band_right  <- band_left + band_frac
  right_left  <- band_right + gap_frac
  
  or_pos        <- right_left + col_pad
  ci_pos        <- or_pos + col_step
  q_pos         <- ci_pos + col_step
  snp_pos       <- q_pos + col_step
  x_right_limit <- snp_pos + 0.015
  
  label_x  <- left_frac * 0.60
  method_x <- left_frac * 0.72
  
  # Map OR, CI to transformed x
  df_multi <- df_multi %>%
    dplyr::mutate(
      OR_plot       = map_shrink(OR,       band_left, band_right),
      OR_lower_plot = map_shrink(OR_lower, band_left, band_right),
      OR_upper_plot = map_shrink(OR_upper, band_left, band_right)
    )
  
  y_min <- min(df_multi$y_pos) - 0.5
  y_max <- max(df_multi$y_pos) + 0.5
  
  panel_df <- data.frame(
    xmin = band_left,
    xmax = band_right,
    ymin = y_min,
    ymax = y_max
  )
  
  # Reference line at OR = 1
  line1_x <- map_shrink(1, band_left, band_right)
  line1_df <- data.frame(
    x    = line1_x,
    xend = line1_x,
    y    = y_min,
    yend = y_max
  )
  
  sig_palette <- c(
    "q < 0.01" = "#d73027",
    "q < 0.05" = "#fc8d59",
    "q ≥ 0.05" = "#999999"
  )
  
  or_ticks        <- c(0.1, 0.2, 0.5, 1, 2, 4, 6)
  x_ticks_plot    <- map_shrink(or_ticks, band_left, band_right)
  tick_baseline_y <- y_min - 2.00
  tick_text_y     <- y_min - 2.30
  
  p_multi <- ggplot(df_multi, aes(y = y_pos)) +
    # Background band
    geom_rect(
      data = panel_df, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "gray98", color = "gray85", alpha = 0.4
    ) +
    # Light vertical guide lines at OR ticks
    geom_segment(
      data = data.frame(x_plot = x_ticks_plot),
      inherit.aes = FALSE,
      aes(x = x_plot, xend = x_plot, y = y_min, yend = y_max),
      color = "gray90", linewidth = 0.4
    ) +
    # Reference line at OR = 1
    geom_segment(
      data = line1_df, inherit.aes = FALSE,
      aes(x = x, xend = xend, y = y, yend = yend),
      linetype = "dashed", color = "gray55", linewidth = 0.6
    ) +
    # Points + CI
    geom_pointrange(
      aes(x = OR_plot, xmin = OR_lower_plot, xmax = OR_upper_plot,
          colour = signif_group),
      size = 0.4, linewidth = 0.4
    ) +
    # Outcome (left)
    geom_text(
      aes(x = label_x, label = outcome_label),
      hjust = 1, vjust = 0.5, size = 3.6, fontface = "bold"
    ) +
    # Method
    geom_text(
      aes(x = method_x, label = method),
      hjust = 0, vjust = 0.5, size = 3.1
    ) +
    # Numeric columns
    geom_text(aes(x = or_pos,  label = OR_text),
              hjust = 0, vjust = 0.5, size = 3.0) +
    geom_text(aes(x = ci_pos,  label = CI_text),
              hjust = 0, vjust = 0.5, size = 3.0) +
    geom_text(aes(x = q_pos,   label = Q_text),
              hjust = 0, vjust = 0.5, size = 3.0) +
    geom_text(aes(x = snp_pos, label = SNP_text),
              hjust = 0, vjust = 0.5, size = 3.0) +
    # Column headers
    annotate(
      "text", x = label_x, y = y_max + 0.7,
      label = "Outcome", hjust = 1, size = 3.8, fontface = "bold"
    ) +
    annotate(
      "text", x = method_x, y = y_max + 0.7,
      label = "Method", hjust = 0, size = 3.8, fontface = "bold"
    ) +
    annotate(
      "text", x = or_pos, y = y_max + 0.7,
      label = "OR", hjust = 0, size = 3.6, fontface = "bold"
    ) +
    annotate(
      "text", x = ci_pos, y = y_max + 0.7,
      label = "95% CI", hjust = 0, size = 3.6, fontface = "bold"
    ) +
    annotate(
      "text", x = q_pos, y = y_max + 0.7,
      label = "q value (FDR)", hjust = 0, size = 3.6, fontface = "bold"
    ) +
    annotate(
      "text", x = snp_pos, y = y_max + 0.7,
      label = "SNPs", hjust = 0, size = 3.6, fontface = "bold"
    ) +
    # Compact OR axis underneath
    annotate(
      "segment",
      x = min(x_ticks_plot), xend = max(x_ticks_plot),
      y = tick_baseline_y, yend = tick_baseline_y,
      linewidth = 0.5, colour = "black"
    ) +
    geom_segment(
      data = data.frame(x_plot = x_ticks_plot, lab = or_ticks),
      inherit.aes = FALSE,
      aes(x = x_plot, xend = x_plot,
          y = tick_baseline_y - 0.06, yend = tick_baseline_y + 0.06),
      linewidth = 0.4, colour = "black"
    ) +
    geom_text(
      data = data.frame(x_plot = x_ticks_plot, lab = or_ticks),
      inherit.aes = FALSE,
      aes(x = x_plot, y = tick_text_y, label = lab),
      size = 3.2, hjust = 0.5
    ) +
    annotate(
      "text",
      x = mean(range(x_ticks_plot)), y = tick_text_y - 0.7,
      label = "Odds ratio", size = 3.8, fontface = "bold", hjust = 0.5
    ) +
    scale_x_continuous(limits = c(0, x_right_limit), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_min - 3.2, y_max + 1.2), expand = c(0, 0)) +
    scale_colour_manual(values = sig_palette, name = "FDR q-value") +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position  = "bottom",
      legend.margin    = margin(t = 40),
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 9),
      plot.margin      = margin(20, 35, 65, 45),
      plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle    = element_text(hjust = 0.5, size = 11)
    ) +
    labs(
      title ="Figure 3. Mendelian randomization estimates for genetically predicted endometriosis liability on pregnancy outcomes across different methods",
      subtitle = subtitle_text
    )
  
  out_file <- file.path(plots_dir, file_stub)
  ggsave(out_file, p_multi, width = 20, height = 30, dpi = 300, bg = "white")
  message("Saved multi-method forest plot: ", out_file)
}

# --------------------------------------------------------------------------- #
# FIGURE 3A – Top 8 IVW P-significant outcomes (P < 0.05)
# --------------------------------------------------------------------------- #

df_multi_top <- build_multi_method_df(
  ivw_df        = ivw_filtered,
  all_df        = all_filtered,
  method_levels = method_levels,
  top_n         = 8,
  p_threshold   = 0.05
)

if (nrow(df_multi_top)) {
  make_multi_forest_plot(
    df_multi_top,
    file_stub     = "Fig3a_forest_multiMethod_topP_endoMR-PREG.png",
    subtitle_text = "Two-sample MR – Multiple methods (top 8 outcomes with P < 0.05, IVW)"
  )
} else {
  warning("No outcomes with P < 0.05 for the top-P multi-method forest plot.")
}

# --------------------------------------------------------------------------- #
# FIGURE 3B – All outcomes (no P filter)
# --------------------------------------------------------------------------- #

df_multi_all <- build_multi_method_df(
  ivw_df        = ivw_filtered,
  all_df        = all_filtered,
  method_levels = method_levels,
  top_n         = NULL,
  p_threshold   = NULL
)

if (nrow(df_multi_all)) {
  make_multi_forest_plot(
    df_multi_all,
    file_stub     = "Fig3b_forest_multiMethod_allOutcomes_endoMR_PREG.png",
    subtitle_text = "Two-sample MR – Multiple methods (all pregnancy outcomes)"
  )
} else {
  warning("No outcomes available for the all-outcomes multi-method forest plot.")
}


##################################################################
### FIGURE 4. TRIOS FOREST PLOT: MATERNAL vs FETAL vs PATERNAL ###
##################################################################

message("Creating trios comparison forest plot (maternal vs fetal vs paternal)...")

if (!file.exists(trios_path)) {
  warning(
    "Trios comparison file not found: ", trios_path,
    "\nMaternal vs fetal vs paternal forest plot will be skipped."
  )
  
} else {
  
  ivw_all <- readr::read_csv(trios_path, show_col_types = FALSE)
  
  # Expected columns
  required_cols <- c("outcome_full", "OR", "LCL", "UCL", "origin")
  if (!all(required_cols %in% names(ivw_all))) {
    warning(
      "Trios file does not have the expected columns: ",
      paste(setdiff(required_cols, names(ivw_all)), collapse = ", "),
      "\nSkipping trios forest plot."
    )
    
  } else {
    
    # ------------------------------------------------------------------
    # 1) Define outcome order (by effect size) and spacing between rows
    # ------------------------------------------------------------------
    
    # Summarise effect size per outcome (max |log(OR)| across origins)
    outcome_order_df <- ivw_all %>%
      dplyr::group_by(outcome_full) %>%
      dplyr::summarise(
        order_val = max(abs(log(OR)), na.rm = TRUE),
        .groups   = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(order_val)) %>%
      dplyr::mutate(order_outcome = dplyr::row_number())
    
    # Join back to full data
    ivw_all <- ivw_all %>%
      dplyr::left_join(outcome_order_df, by = "outcome_full")
    
    # Order of outcomes on the y-axis
    outcome_levels <- outcome_order_df$outcome_full
    
    # Define a large vertical step between outcomes (controls gap between variables)
    group_step <- 2    # increase for more spacing (e.g. 2.0, 2.5, 3.0)
    
    base_pos <- outcome_order_df %>%
      dplyr::mutate(
        base_y = group_step * (order_outcome - 1L)  # 0, 2, 4, 6, ...
      ) %>%
      dplyr::select(outcome_full, base_y)
    
    ivw_all <- ivw_all %>%
      dplyr::left_join(base_pos, by = "outcome_full")
    
    # ------------------------------------------------------------------
    # 2) Within each outcome: equidistant Maternal / Fetal / Paternal
    # ------------------------------------------------------------------
    
    origin_levels <- c("Maternal", "Fetal", "Paternal")
    origin_offsets <- c(-0.3, 0, 0.3)   # vertical offsets around base_y
    
    ivw_all <- ivw_all %>%
      dplyr::mutate(
        origin = factor(origin, levels = origin_levels),
        origin_index = as.integer(origin),                 # 1, 2, 3
        y_pos = base_y + origin_offsets[origin_index]      # base_y ± offset
      )
    
    # ------------------------------------------------------------------
    # 3) Build the forest plot
    # ------------------------------------------------------------------
    
    p_trios <- ggplot(
      ivw_all,
      aes(
        x = OR,
        y = y_pos,
        colour = origin
      )
    ) +
      # Reference line OR = 1
      geom_vline(
        xintercept = 1,
        linetype   = "dashed",
        colour     = "grey60",
        linewidth  = 0.6
      ) +
      # Confidence intervals
      geom_errorbarh(
        aes(xmin = LCL, xmax = UCL),
        height    = 0.25,
        linewidth = 0.6
      ) +
      # Points
      geom_point(
        size = 2.8
      ) +
      # X-axis on log scale
      scale_x_log10(
        breaks = c(0.5, 0.7, 1, 1.4, 2),
        labels = c("0.5", "0.7", "1.0", "1.4", "2.0"),
        limits = c(0.45, 2.2)
      ) +
      
      # Y-axis: show only one label per outcome (at base_y)
      scale_y_continuous(
        breaks = base_pos$base_y,
        labels = base_pos$outcome_full,
        expand = expansion(mult = c(0.05, 0.05))
      ) +
      # Colours for maternal / fetal / paternal
      scale_colour_manual(
        values = c(
          "Maternal" = "#1f78b4",
          "Fetal"    = "#33a02c",
          "Paternal" = "#e31a1c"
        ),
        name = "Genetic effect"
      ) +
      labs(
        title = "Figure 5. Mendelian randomisation estimates for maternal blood pressure",
        subtitle = "Effects on primary pregnancy outcomes, adjusted for offspring genotype",
        x        = "Odds ratio (log scale, 95% CI)",
        y        = NULL,
        caption  = paste(
          "Odds ratios and 95% confidence intervals for the effect of genetically predicted endometriosis liability",
          "on pregnancy outcomes. Maternal, fetal, and paternal genetic effects were estimated using the same set",
          "of independent endometriosis-associated variants (Rahmioglu et al. 2023).",
          sep = "\n"
        )
      ) +
      theme_minimal(base_size = 11) +
      theme(
        # Titles
        plot.title    = element_text(size = 14, face = "bold", hjust = 0.5, lineheight = 1.1),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        
        # Axis text: we now label only outcomes at base_y
        axis.text.y = element_text(size = 9, margin = margin(t = 4, b = 4)),
        axis.text.x = element_text(size = 9),
        
        # Legend at bottom, left aligned
        legend.position      = "bottom",
        legend.justification = c(0, 0.5),
        legend.box.just      = "left",
        legend.title         = element_text(size = 10, face = "bold"),
        legend.text          = element_text(size = 9),
        
        # Caption
        plot.caption = element_text(
          hjust  = 0,
          size   = 8,
          colour = "grey40",
          margin = margin(t = 10)
        ),
        
        # Grid
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        
        # Margins
        plot.margin = margin(t = 15, r = 25, b = 20, l = 20)
      )
    
    # ------------------------------------------------------------------
    # 4) Save figure
    # ------------------------------------------------------------------
    
    out_png <- file.path(plots_dir, "Fig4_trios_maternal_fetal_paternal.png")
    
    ggsave(out_png, p_trios, width = 14, height = 12, dpi = 300, bg = "white")
    
    message("Saved trios comparison plot to:")
    message("  - ", out_png)
  }
}


### 6) SUMMARY ##################################################################

message("\n=== PLOT GENERATION COMPLETED ===")
message("Figures saved in: ", plots_dir)
message(" - Multi-method forest (top FDR): Fig3a_forest_multiMethod_topP_endoMR-PREG.png")
message(" - Multi-method forest (all outcomes): Fig3b_forest_multiMethod_allOutcomes_endoMR_PREG.png")
message(" - Trios comparison (if file present): Fig4_trios_maternal_fetal_paternal.png")
