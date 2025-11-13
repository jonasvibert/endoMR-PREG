#!/usr/bin/env Rscript
###############################################################################
# Visual MR Forest Plots for endoMR-PREG (pregnancy outcomes)
# Forest plot scripts for pregnancy outcomes
# - Creates IVW forest plots (single method, clean layout)
# - Creates multi-method forest plots (IVW, Egger, Weighted median, etc.)
# - Focuses on pregnancy, birth, and maternal health outcomes
# Authors: Claude Code
# Date: 2025-10-11
###############################################################################

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(readr)
  library(tidyr)
  library(stringr)
  library(forcats)
})

# Set working directory to endoMR-PREG
setwd("~/Desktop/GitHub/endoMR-PREG")

# ---------- Paths ----------
results_dir <- here("results")
plots_dir   <- here("plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Helpers ----------
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

message("=== FOREST PLOTS FOR endoMR-PREG ===")
message("Loading MR results and generating forest plots...")

# ---------- Load MR results ----------
ivw_results_path <- file.path(results_dir, "ivw_results.csv")
all_results_path <- file.path(results_dir, "all_mr_methods.csv")

if (!file.exists(ivw_results_path)) {
  stop("IVW results not found: ", ivw_results_path, "\nRun 4_MR_analysis_endoMR-PREG.R first.")
}
if (!file.exists(all_results_path)) {
  stop("All methods results not found: ", all_results_path, "\nRun 4_MR_analysis_endoMR-PREG.R first.")
}

ivw_results <- readr::read_csv(ivw_results_path, show_col_types = FALSE)
all_results <- readr::read_csv(all_results_path, show_col_types = FALSE)

message("Loaded IVW results: ", nrow(ivw_results), " outcomes")
message("Loaded all methods: ", nrow(all_results), " results")

# ---------- Outcome Categories and Labels ----------
# Define pregnancy/maternal health outcome categories and clean labels (filtered to 29 outcomes)
outcome_labels <- c(
  # Placental outcomes
  "finngen_R12_O15_PLAC_PRAEVIA"    = "Placenta praevia",
  "finngen_R12_O15_PLAC_DISORD"     = "Placental disorders",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR" = "Premature placental separation",
  
  # Cesarean delivery
  "el_cs"                  = "Elective caesarean section",
  "em_cs"                  = "Emergency caesarean section",
  "cs"                     = "Caesarean section",
  
  # Apgar scores
  "lowapgar1"              = "Low Apgar score at 1 min",
  "lowapgar5"              = "Low Apgar score at 5 min",
  
  # Labor and delivery complications
  "rup_memb"               = "Premature rupture of membranes",
  "induction"              = "Labour induction",
  
  # Gestational age and timing
  "ga_all"                 = "Gestational age",
  "pretb_all"              = "Preterm birth (any)",
  "vpretb_all"             = "Very preterm birth",
  "posttb_all"             = "Post-term birth",
  
  # Birth weight outcomes
  "hbw_all"                = "High birthweight (>4000g)",
  "lbw_all"                = "Low birthweight (<2500g)",
  "sga"                    = "Small for gestational age",
  
  # Maternal health
  "depr_subsamp"           = "Postpartum Depression",
  "anaemia_preg_all"       = "Pregnancy anemia",
  
  # Pregnancy complications
  "gdm_subsamp"            = "Gestational diabetes mellitus",
  "hdp_subsamp"            = "Hypertensive disorders of pregnancy",
  "gh_subsamp"             = "Gestational hypertension",
  "pe_subsamp"             = "Preeclampsia",
  
  # Neonatal outcomes
  "nicu"                   = "NICU admission",
  "sb_subsamp"             = "Stillbirth",
  
  # Hemorrhage and bleeding
  "Antepartum_bleeding_filtered"     = "Antepartum bleeding",
  "Postpartum_hemorrhage_filtered"   = "Postpartum hemorrhage",
  "Postpartum_hemorrhage_due_to_atony_filtered" = "PPH due to atony",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered" = "PPH due to retained placenta"
)

# ---------- Filter for pregnancy-relevant outcomes ----------
pregnancy_outcomes <- names(outcome_labels)
ivw_filtered <- ivw_results %>%
  filter(outcome %in% pregnancy_outcomes) %>%
  mutate(outcome_clean = dplyr::coalesce(outcome_labels[outcome], outcome_full, outcome))

all_filtered <- all_results %>%
  filter(outcome %in% pregnancy_outcomes) %>%
  mutate(outcome_clean = dplyr::coalesce(outcome_labels[outcome], outcome_full, outcome))

message("Filtered to ", nrow(ivw_filtered), " pregnancy-related outcomes for IVW")
message("Filtered to ", nrow(all_filtered), " results for multi-method analysis")

# ============================================================================ #
# FOREST PLOT 1: Clean IVW Forest Plot
# Clean forest plot style
# ============================================================================ #

message("Creating IVW forest plot...")

# Define specific outcomes (filtered to 29 outcomes)
target_outcomes <- c(
  "finngen_R12_O15_PLAC_PRAEVIA",    # Placenta praevia
  "finngen_R12_O15_PLAC_DISORD",     # Placental disorders  
  "finngen_R12_O15_PLAC_PREMAT_SEPAR", # Premature placental separation
  "el_cs",                           # Elective caesarean section
  "lowapgar1",                       # Low Apgar score at 1 min
  "em_cs",                           # Emergency caesarean section
  "lowapgar5",                       # Low Apgar score at 5 min
  "rup_memb",                        # Premature rupture of membranes
  "ga_all",                          # Gestational age
  "Antepartum_bleeding_filtered",    # Antepartum bleeding
  "hbw_all",                         # High birthweight (>4000g)
  "posttb_all",                      # Post-term birth
  "depr_subsamp",                    # Depression
  "cs",                              # Caesarean section
  "induction",                       # Labour induction
  "anaemia_preg_all",                # Pregnancy anemia
  "hdp_subsamp",                     # Hypertensive disorders of pregnancy
  "gdm_subsamp",                     # Gestational diabetes mellitus
  "Postpartum_hemorrhage_due_to_atony_filtered", # PPH due to atony
  "nicu",                            # NICU admission
  "gh_subsamp",                      # Gestational hypertension
  "Postpartum_hemorrhage_filtered",  # Postpartum hemorrhage
  "pe_subsamp",                      # Preeclampsia
  "sb_subsamp",                      # Stillbirth
  "lbw_all",                         # Low birthweight (<2500g)
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered", # PPH due to retained placenta
  "pretb_all",                       # Preterm birth (any)
  "sga",                             # Small for gestational age 
  "vpretb_all"                       # Very preterm birth
)

# Prepare IVW data with specific outcomes only
df_ivw <- ivw_filtered %>%
  filter(method == "Inverse variance weighted") %>%
  filter(outcome %in% target_outcomes) %>%  # Filter for specific outcomes only
  mutate(
    OR = exp(b),
    CI_low = exp(b - 1.96 * se),
    CI_high = exp(b + 1.96 * se),
    OR_txt = sprintf("%.2f", OR),
    CI_txt = sprintf("(%.2f, %.2f)", CI_low, CI_high),
    p_txt = case_when(
      pval < 0.001 ~ "<0.001",
      pval < 0.01 ~ sprintf("%.3f", pval),
      TRUE ~ sprintf("%.3f", pval)
    ),
    nsnp_txt = as.character(nsnp),
    # Significance colors
    sig_group = case_when(
      pval < 0.001 ~ "p < 0.001",
      pval < 0.01 ~ "p < 0.01", 
      pval < 0.05 ~ "p < 0.05",
      TRUE ~ "p ≥ 0.05"
    )
  ) %>%
  # Order by effect size (OR) for visual impact - descending order as in original image
  arrange(desc(OR))

if (nrow(df_ivw) == 0) {
  stop("No IVW pregnancy outcomes found in results.")
}

# Lock factor order for consistent plotting
df_ivw$outcome_clean <- factor(df_ivw$outcome_clean, levels = rev(unique(df_ivw$outcome_clean)))

# Axis limits and positioning (enhanced layout)
xmin <- 0.45
xmax_target <- 2.0
xmax_data <- max(df_ivw$CI_high, na.rm = TRUE) * 1.15
xmax <- max(xmax_target, xmax_data) * 2.1   # extra canvas; pushes forest band left

# Column positions (kept ≤ 1.0 * xmax so nothing is clipped)
x_or   <- 0.58 * xmax   # OR a bit right
x_ci   <- 0.76 * xmax   # CI to the right of OR
x_p    <- 0.90 * xmax   # p-value well separated
x_snp  <- 0.985 * xmax  # SNPs at the far right, still inside xlim

# Headers
y_hdr <- length(levels(df_ivw$outcome_clean)) + 1.0
hdr_outcome <- data.frame(x = xmin, y = y_hdr, lab = "Outcome")
hdr_cols <- data.frame(
  x = c(x_or, x_ci, x_p, x_snp),
  y = rep(y_hdr, 4),
  lab = c("OR", "95% CI", "p-value", "SNPs")
)

# Build forest plot
p_ivw <- ggplot(df_ivw, aes(y = outcome_clean)) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, color = sig_group),
                 height = 0.18, linewidth = 0.7) +
  geom_point(aes(x = OR, color = sig_group), size = 2.3) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5, color = "gray60") +
  
  # Text columns (with ample spacing; prevent overlaps)
  geom_text(aes(x = x_or,  label = OR_txt),   hjust = 0, family = "mono", size = 3.2, check_overlap = TRUE) +
  geom_text(aes(x = x_ci,  label = CI_txt),   hjust = 0, family = "mono", size = 3.2, check_overlap = TRUE) +
  geom_text(aes(x = x_p,   label = p_txt),    hjust = 0, family = "mono", size = 3.2, check_overlap = TRUE) +
  geom_text(aes(x = x_snp, label = nsnp_txt), hjust = 0, family = "mono", size = 3.2, check_overlap = TRUE) +
  
  # Headers
  geom_text(data = hdr_outcome, aes(x = x, y = y, label = lab),
            inherit.aes = FALSE, hjust = 0, fontface = 2, size = 4.6) +
  geom_text(data = hdr_cols, aes(x = x, y = y, label = lab),
            inherit.aes = FALSE, hjust = 0, fontface = 2, size = 4.0, family = "mono") +
  
  # Scales and theme
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2),
    limits = c(xmin, xmax),
    minor_breaks = NULL
  ) +
  scale_color_manual(values = c(
    "p < 0.001" = "#d73027",
    "p < 0.01"  = "#fc8d59",
    "p < 0.05"  = "#fee08b",
    "p ≥ 0.05"  = "#999999"
  ), name = "Significance") +
  labs(
    x = "Odds Ratio",
    y = NULL,
    title = "Association of Genetically Predicted Endometriosis with Pregnancy Outcomes",
    subtitle = "Two-sample Mendelian randomization (Inverse variance weighted method)"
  ) +
  coord_cartesian(xlim = c(xmin, xmax), clip = "off", expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 11, margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 8)),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(20, 360, 20, 20)   # generous right margin to show all text
  )


# Save IVW forest plot
ggsave(file.path(plots_dir, "forest_endoMR-PREG_IVW.png"),
       p_ivw, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: ", file.path(plots_dir, "forest_endoMR-PREG_IVW.png"))

# ============================================================================ #
# FOREST PLOT 2: Multi-Method Forest Plot  
# Multi-method forest plot approach
# ============================================================================ #

message("Creating multi-method forest plot...")

# Method levels for consistent ordering
method_levels <- c("MR Egger", "Weighted median", "Inverse variance weighted", "Weighted mode")

# Prepare multi-method data - focus on significant outcomes
significant_outcomes <- ivw_filtered %>%
  filter(pval < 0.05) %>%
  arrange(pval) %>%
  slice_head(n = 8) %>%  # Top 8 most significant
  pull(outcome)

df_multi <- all_filtered %>%
  filter(
    method %in% method_levels,
    outcome %in% significant_outcomes
  ) %>%
  mutate(
    OR = exp(as.numeric(b)),
    OR_lower = exp(as.numeric(b) - 1.96 * as.numeric(se)),
    OR_upper = exp(as.numeric(b) + 1.96 * as.numeric(se)),
    OR_text = sprintf("%.2f", OR),
    CI_text = sprintf("(%.2f–%.2f)", OR_lower, OR_upper),
    P_text = case_when(
      as.numeric(pval) < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", as.numeric(pval))
    ),
    SNP_text = paste0("n=", nsnp),
    method = factor(method, levels = method_levels),
    signif_group = case_when(
      as.numeric(pval) < 0.001 ~ "p < 0.001",
      as.numeric(pval) < 0.01 ~ "p < 0.01",
      as.numeric(pval) < 0.05 ~ "p < 0.05", 
      TRUE ~ "p ≥ 0.05"
    )
  )

if (nrow(df_multi) == 0) {
  warning("No significant outcomes found for multi-method plot")
} else {
  
  # Order outcomes alphabetically
  outcome_levels <- df_multi %>% distinct(outcome_clean) %>% arrange(outcome_clean) %>% pull(outcome_clean)
  df_multi$outcome_clean <- factor(df_multi$outcome_clean, levels = outcome_levels)
  
  # Vertical layout: 4 methods per outcome with gaps
  n_methods <- length(method_levels)
  gap <- 0.8
  base_pos <- data.frame(
    outcome_clean = outcome_levels,
    base_y = cumsum(c(0, rep(n_methods + gap, length(outcome_levels) - 1)))
  )
  
  df_multi <- df_multi %>%
    arrange(outcome_clean, method) %>%
    group_by(outcome_clean) %>%
    mutate(order_in_outcome = row_number()) %>%
    ungroup() %>%
    left_join(base_pos, by = "outcome_clean") %>%
    mutate(
      y_pos = base_y + (n_methods - order_in_outcome + 1),
      outcome_label = ifelse(order_in_outcome == 1, as.character(outcome_clean), "")
    )
  
  # Shrink transform for OR axis
  x_show_min <- 0.1
  x_show_max <- 6
  shrink_factor <- 0.25
  shrink_f <- function(x) 10^(log10(x) * shrink_factor)
  
  map_shrink <- function(x, band_left, band_right) {
    z <- shrink_f(x)
    z0 <- shrink_f(x_show_min)
    z1 <- shrink_f(x_show_max)
    band_left + (z - z0) / (z1 - z0) * (band_right - band_left)
  }
  
  # Layout fractions
  left_frac <- 0.25
  gap_label_to_band <- 0.06
  band_frac <- 0.30
  gap_frac <- 0.03
  
  col_step <- 0.09
  col_pad <- 0.010
  
  # Compute positions
  band_left <- left_frac + gap_label_to_band
  band_right <- band_left + band_frac
  right_left <- band_right + gap_frac
  
  or_pos <- right_left + col_pad
  ci_pos <- or_pos + col_step
  p_pos <- ci_pos + col_step
  snp_pos <- p_pos + col_step
  x_right_limit <- snp_pos + 0.015
  
  label_x <- left_frac * 0.60
  method_x <- left_frac * 0.72
  
  # Transform data to plotting coordinates
  df_multi <- df_multi %>%
    mutate(
      OR_plot = map_shrink(OR, band_left, band_right),
      OR_lower_plot = map_shrink(OR_lower, band_left, band_right),
      OR_upper_plot = map_shrink(OR_upper, band_left, band_right)
    )
  
  # Plot elements
  y_min <- min(df_multi$y_pos) - 0.5
  y_max <- max(df_multi$y_pos) + 0.5
  
  panel_df <- data.frame(xmin = band_left, xmax = band_right, ymin = y_min, ymax = y_max)
  line1_x <- map_shrink(1, band_left, band_right)
  line1_df <- data.frame(x = line1_x, xend = line1_x, y = y_min, yend = y_max)
  
  # Significance color palette
  sig_palette <- c("p < 0.001" = "#d73027", "p < 0.01" = "#fc8d59", 
                   "p < 0.05" = "#fee08b", "p ≥ 0.05" = "#999999")
  
  # Ticks for OR scale
  or_ticks <- c(0.1, 0.2, 0.5, 1, 2, 4, 6)
  x_ticks_plot <- map_shrink(or_ticks, band_left, band_right)
  tick_baseline_y <- y_min - 2.00
  tick_text_y <- y_min - 2.30
  
  # Build multi-method forest plot
  p_multi <- ggplot(df_multi, aes(y = y_pos)) +
    # Background
    geom_rect(data = panel_df, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "gray98", color = "gray85", alpha = 0.4) +
    
    # Vertical guides
    geom_segment(data = data.frame(x_plot = x_ticks_plot),
                 inherit.aes = FALSE,
                 aes(x = x_plot, xend = x_plot, y = y_min, yend = y_max),
                 color = "gray90", linewidth = 0.4) +
    
    # OR = 1 reference line
    geom_segment(data = line1_df, inherit.aes = FALSE,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 linetype = "dashed", color = "gray55", linewidth = 0.6) +
    
    # Data points and CIs
    geom_pointrange(aes(x = OR_plot, xmin = OR_lower_plot, xmax = OR_upper_plot, color = signif_group),
                    size = 0.4, linewidth = 0.4) +
    
    # Labels
    geom_text(aes(x = label_x, label = outcome_label),
              hjust = 1, vjust = 0.5, size = 3.6, fontface = "bold") +
    geom_text(aes(x = method_x, label = method),
              hjust = 0, vjust = 0.5, size = 3.1) +
    
    # Right columns
    geom_text(aes(x = or_pos, label = OR_text), hjust = 0, vjust = 0.5, size = 3.0) +
    geom_text(aes(x = ci_pos, label = CI_text), hjust = 0, vjust = 0.5, size = 3.0) +
    geom_text(aes(x = p_pos, label = P_text), hjust = 0, vjust = 0.5, size = 3.0) +
    geom_text(aes(x = snp_pos, label = SNP_text), hjust = 0, vjust = 0.5, size = 3.0) +
    
    # Headers
    annotate("text", x = label_x, y = y_max + 0.7, label = "Outcome",
             hjust = 1, size = 3.8, fontface = "bold") +
    annotate("text", x = method_x, y = y_max + 0.7, label = "Method", 
             hjust = 0, size = 3.8, fontface = "bold") +
    annotate("text", x = or_pos, y = y_max + 0.7, label = "OR",
             hjust = 0, size = 3.6, fontface = "bold") +
    annotate("text", x = ci_pos, y = y_max + 0.7, label = "95% CI",
             hjust = 0, size = 3.6, fontface = "bold") +
    annotate("text", x = p_pos, y = y_max + 0.7, label = "P value",
             hjust = 0, size = 3.6, fontface = "bold") +
    annotate("text", x = snp_pos, y = y_max + 0.7, label = "SNPs",
             hjust = 0, size = 3.6, fontface = "bold") +
    
    # Tick legend
    annotate("segment",
             x = min(x_ticks_plot), xend = max(x_ticks_plot),
             y = tick_baseline_y, yend = tick_baseline_y,
             linewidth = 0.5, colour = "black") +
    geom_segment(data = data.frame(x_plot = x_ticks_plot, lab = or_ticks),
                 inherit.aes = FALSE,
                 aes(x = x_plot, xend = x_plot, 
                     y = tick_baseline_y - 0.06, yend = tick_baseline_y + 0.06),
                 linewidth = 0.4, colour = "black") +
    geom_text(data = data.frame(x_plot = x_ticks_plot, lab = or_ticks),
              inherit.aes = FALSE,
              aes(x = x_plot, y = tick_text_y, label = lab),
              size = 3.2, hjust = 0.5) +
    annotate("text", x = mean(range(x_ticks_plot)), y = tick_text_y - 0.7,
             label = "Odds Ratio", size = 3.8, fontface = "bold", hjust = 0.5) +
    
    # Scales and theme
    scale_x_continuous(limits = c(0, x_right_limit), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_min - 3.2, y_max + 1.2), expand = c(0, 0)) +
    scale_color_manual(values = sig_palette, name = "Significance") +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.margin = margin(t = 40),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.margin = margin(20, 35, 65, 45),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    ) +
    labs(
      title = "Association of Genetically Predicted Endometriosis with Pregnancy Outcomes",
      subtitle = "Two-sample MR — Multiple methods comparison"
    )
  
  # Save multi-method forest plot
  ggsave(file.path(plots_dir, "forest_endoMR-PREG_multiMethod.png"),
         plot = p_multi, width = 13, height = 10, dpi = 300, bg = "white")
  message("Saved: ", file.path(plots_dir, "forest_endoMR-PREG_multiMethod.png"))
}

# ============================================================================ #
# Summary and completion
# ============================================================================ #

message("\n=== FOREST PLOTS COMPLETED ===")
message("Generated forest plots:")
message("1. IVW forest plot: forest_endoMR-PREG_IVW.png")
if (nrow(df_multi) > 0) {
  message("2. Multi-method forest plot: forest_endoMR-PREG_multiMethod.png")
}
message("All plots saved to: ", plots_dir)
message("Forest plot generation completed successfully!")

# Mark task as completed
message("Script 8 forest plots analysis completed.")

# ============================================================================ #
# FOREST PLOT 3: Improved Forest Plot Layout
# Better spacing, alignment, and label readability
# ============================================================================ #

message("Creating improved forest plot...")

# Prepare data using the same target_outcomes as defined earlier
df_improved_style <- ivw_filtered %>%
  filter(method == "Inverse variance weighted") %>%
  filter(outcome %in% target_outcomes) %>%  # Use the same target outcomes
  mutate(
    OR = exp(b),
    OR_lower = exp(b - 1.96 * se),
    OR_upper = exp(b + 1.96 * se),
    OR_text = sprintf("%.2f", OR),
    CI_text = sprintf("(%.2f-%.2f)", OR_lower, OR_upper), 
    P_text = ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval)),
    SNP_text = paste0("n=", nsnp),
    color_group = case_when(
      pval < 0.001 ~ "p < 0.001",
      pval < 0.01 ~ "p < 0.01", 
      pval < 0.05 ~ "p < 0.05",
      TRUE ~ "p ≥ 0.05"
    )
  ) %>%
  arrange(desc(OR)) %>%  # Sort by descending OR
  mutate(y_pos = n() - row_number() + 1)  # Largest OR gets highest y_pos (top)

if (nrow(df_improved_style) == 0) {
  message("No data available for improved forest plot")
} else {
  
  # Optional: wrap long outcome labels for better readability
  df_improved_style$outcome_wrap <- str_wrap(df_improved_style$outcome_clean, width = 32)
  
  # Apply spacing factor to y positions for better readability
  y_spacing_factor <- 1.3
  df_improved_style$y_pos <- df_improved_style$y_pos * y_spacing_factor
  
  # Dynamic limits for the odds ratio panel (with small padding)
  panel_pad_left  <- 0.05
  panel_pad_right <- 0.05
  true_min <- min(df_improved_style$OR_lower, na.rm = TRUE)
  true_max <- max(df_improved_style$OR_upper, na.rm = TRUE)
  x_min <- max(0.1, true_min - panel_pad_left)
  x_max <- true_max + panel_pad_right
  
  # Adjust label and text column positions
  x_label_pos <- x_min - 0.40   # More space for outcome labels
  text_start  <- x_max + 0.10
  text_spacing <- 0.35          # increased for better spacing
  
  # Text column positions
  or_pos  <- text_start
  ci_pos  <- text_start + text_spacing
  p_pos   <- text_start + 2 * text_spacing
  snp_pos <- text_start + 3 * text_spacing
  
  # Panel height and header with increased spacing
  y_bottom <- 0.5
  y_top <- max(df_improved_style$y_pos) + 0.5
  y_header <- max(df_improved_style$y_pos) + 1
  
  # Build improved forest plot
  p_forest <- ggplot(df_improved_style, aes(y = y_pos)) +
    
    # Background panel
    geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_bottom, ymax = y_top),
              fill = "gray98", color = "gray85", alpha = 0.3, inherit.aes = FALSE) +
    
    # Reference line (OR = 1)
    geom_segment(aes(x = 1, xend = 1, y = y_bottom, yend = y_top),
                 linetype = "dashed", color = "gray60", linewidth = 0.5, inherit.aes = FALSE) +
    
    # Points and confidence intervals
    geom_pointrange(aes(x = OR, xmin = OR_lower, xmax = OR_upper, color = color_group),
                    size = 0.8, linewidth = 0.6) +
    
    # Outcome labels (left)
    geom_text(aes(x = x_label_pos, label = outcome_wrap),
              hjust = 1, vjust = 0.5, size = 3.6, lineheight = 1.2) +
    
    # Text columns on the right
    geom_text(aes(x = or_pos,  label = OR_text),  hjust = 0, vjust = 0.5, size = 3.2) +
    geom_text(aes(x = ci_pos,  label = CI_text),  hjust = 0, vjust = 0.5, size = 3.2) +
    geom_text(aes(x = p_pos,   label = P_text),   hjust = 0, vjust = 0.5, size = 3.2) +
    geom_text(aes(x = snp_pos, label = SNP_text), hjust = 0, vjust = 0.5, size = 3.2) +
    
    # Column headers
    annotate("text", x = x_label_pos, y = y_header, label = "Outcome",
             hjust = 1, vjust = 0.5, size = 3.8, fontface = "bold") +
    annotate("text", x = or_pos,  y = y_header, label = "OR",
             hjust = 0, vjust = 0.5, size = 3.8, fontface = "bold") +
    annotate("text", x = ci_pos,  y = y_header, label = "95% CI",
             hjust = 0, vjust = 0.5, size = 3.8, fontface = "bold") +
    annotate("text", x = p_pos,   y = y_header, label = "P value",
             hjust = 0, vjust = 0.5, size = 3.8, fontface = "bold") +
    annotate("text", x = snp_pos, y = y_header, label = "SNPs",
             hjust = 0, vjust = 0.5, size = 3.8, fontface = "bold") +
    
    # Axis label for Odds Ratio
    annotate("text", x = (x_min + x_max) / 2, y = -0.3, label = "Odds Ratio",
             hjust = 0.5, vjust = 0.5, size = 3.8, fontface = "bold") +
    
    # OR scale ticks and labels
    annotate("segment", 
             x = c(0.5, 0.75, 1, 1.25, 1.5, 2.0), xend = c(0.5, 0.75, 1, 1.25, 1.5, 2.0),
             y = -0.8, yend = -0.6,
             linewidth = 0.4, colour = "black") +
    annotate("text", x = c(0.5, 0.75, 1, 1.25, 1.5, 2.0), y = -1.0, 
             label = c("0.5", "0.75", "1.0", "1.25", "1.5", "2.0"), 
             size = 3.2, hjust = 0.5) +
    
    # Scales and theme
    scale_x_continuous(limits = c(x_label_pos - 0.5, snp_pos + 0.4), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1.2, y_header + 1), expand = c(0, 0)) +
    scale_color_manual(values = c(
      "p < 0.001" = "#d73027",
      "p < 0.01"  = "#fc8d59",
      "p < 0.05"  = "#fee08b",
      "p ≥ 0.05"  = "#999999"
    ), name = "Significance") +
    
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.margin = margin(t = 10),
      plot.margin = margin(20, 30, 20, 60),  # Plus de marge à gauche
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 18)),
      plot.subtitle = element_text(hjust = 0.5, size = 11, margin = margin(b = 12))
    ) +
    labs(
      title = "Association of Genetically Predicted Endometriosis with Pregnancy Outcomes",
      subtitle = "Two-sample Mendelian randomization (Inverse variance weighted method)"
    )
  
  # Save improved forest plot
  ggsave(file.path(plots_dir, "forest_endoMR-PREG_improved.png"),
         plot = p_forest, width = 15, height = 10, dpi = 300, bg = "white")
  message("Saved: ", file.path(plots_dir, "forest_endoMR-PREG_improved.png"))
}

# ============================================================================ #
# FOREST PLOT 4: endoPAIN IVW Layout  (log2-transformed coordinates)
# Ensures equal visual distance 0.5→1 and 1→2, preserves visual layout
# ============================================================================ #

message("Creating endoPAIN layout forest plot...")

# Prepare data using the same target_outcomes as defined earlier
df_endopain_layout <- ivw_filtered %>%
  filter(method == "Inverse variance weighted") %>%
  filter(outcome %in% target_outcomes) %>%  # Use the same target outcomes
  mutate(
    OR = exp(b),
    OR_lower = exp(b - 1.96 * se),
    OR_upper = exp(b + 1.96 * se),
    
    # --- Transform to log2 space so 0.5→1 equals 1→2 visually ---
    OR_t       = log2(OR),
    OR_lower_t = log2(OR_lower),
    OR_upper_t = log2(OR_upper),
    
    OR_text = sprintf("%.2f", OR),
    CI_text = sprintf("(%.2f-%.2f)", OR_lower, OR_upper),
    P_text = ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval)),
    SNP_text = paste0("n=", nsnp),
    color_group = case_when(
      pval < 0.001 ~ "p < 0.001",
      pval < 0.01 ~ "p < 0.01", 
      pval < 0.05 ~ "p < 0.05",
      TRUE ~ "p ≥ 0.05"
    )
  ) %>%
  arrange(desc(OR)) %>%  # Sort by descending OR
  mutate(y_pos = n() - row_number() + 1)  # Largest OR gets highest y_pos (top)

# If outcome_clean is missing, fall back to outcome
if (!"outcome_clean" %in% names(df_endopain_layout)) {
  df_endopain_layout <- df_endopain_layout %>% mutate(outcome_clean = outcome)
}

if (nrow(df_endopain_layout) == 0) {
  message("No data available for endoPAIN layout forest plot")
} else {
  
  # Rename for endoPAIN code compatibility
  forest_data <- df_endopain_layout
  
  # Calculate dynamic height based on number of variables
  n_outcomes <- nrow(forest_data)
  base_height <- 8  # endoPAIN base height
  height_adjustment <- max(0, (n_outcomes - 20) * 0.3)  # 0.3 per additional variable
  final_height <- base_height + height_adjustment
  
  # 5) Create forest plot ------------------------------------------------------
  # Define x-axis limits and positions (in OR units)
  x_min <- 0.5
  x_max <- 2.5
  
  # --- Convert layout anchors to transformed (log2) space ---
  x_min_t <- log2(x_min)          # -1
  x_max_t <- log2(x_max)          # ~1.3219
  text_start_t <- x_max_t + 0.05  # Small margin after the forest panel
  text_spacing_t <- 0.25          # Spacing between text columns (in transformed units)
  
  # Text column positions (all in transformed units)
  or_pos_t  <- text_start_t
  ci_pos_t  <- text_start_t + text_spacing_t
  p_pos_t   <- text_start_t + 2 * text_spacing_t
  snp_pos_t <- text_start_t + 3 * text_spacing_t
  
  # ---- Precompute constants for guides/headers -------------------------------
  y_bottom <- 0.5
  y_top    <- max(forest_data$y_pos) + 0.5
  y_header <- max(forest_data$y_pos) + 1
  
  # Panel (now in transformed x)
  panel_df <- data.frame(
    xmin = x_min_t, xmax = x_max_t,
    ymin = y_bottom, ymax = y_top
  )
  
  # Reference line at OR = 1  -> log2(1) = 0
  ref1_df <- data.frame(
    x = 0, xend = 0,
    y = y_bottom, yend = y_top
  )
  
  # Compact tick baseline (transformed)
  axis_base_df <- data.frame(
    x = x_min_t, xend = x_max_t, y = 0.3, yend = 0.3
  )
  
  # Tick marks (choose any ORs you want labeled)
  tick_vals <- c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
  tick_df <- data.frame(x = log2(tick_vals))
  tick_lab_df <- transform(tick_df, lab = sprintf("%g", 2^x))
  
  header_df <- data.frame(
    x   = c(x_min_t - 0.15, or_pos_t,  ci_pos_t,   p_pos_t,    snp_pos_t),
    lab = c("Outcome",      "OR",      "95% CI",   "P value",  "SNPs")
  )
  
  # How far to extend x to include right-hand columns
  x_limit_right <- snp_pos_t + 0.30
  
  # ---- Plot -------------------------------------------------------------------
  p <- ggplot(forest_data, aes(y = y_pos)) +
    # Background panel
    geom_rect(data = panel_df, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "gray98", color = "gray85", alpha = 0.3) +
    
    # OR = 1 dashed reference
    geom_segment(data = ref1_df, inherit.aes = FALSE,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 linetype = "dashed", color = "gray60", linewidth = 0.5) +
    
    # CI ranges and points (NOTE: use transformed columns)
    geom_pointrange(
      aes(x = OR_t, xmin = OR_lower_t, xmax = OR_upper_t, color = color_group),
      size = 0.8, linewidth = 0.6
    ) +
    
    # Compact OR tick line underneath (manual axis, transformed coords)
    geom_segment(data = axis_base_df, inherit.aes = FALSE,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black", linewidth = 0.5) +
    geom_segment(data = transform(tick_df, xend = x, y = 0.25, yend = 0.35),
                 inherit.aes = FALSE,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black", linewidth = 0.4) +
    geom_text(data = transform(tick_lab_df, y = 0.1),
              inherit.aes = FALSE,
              aes(x = x, y = y, label = lab),
              hjust = 0.5, vjust = 0.5, size = 3.2) +
    
    # Outcome labels (left) - more space for long names
    geom_text(aes(x = x_min_t - 0.15, label = outcome_clean),
              hjust = 1, vjust = 0.5, size = 3.5) +
    
    # Right-side text columns (positions in transformed space)
    geom_text(aes(x = or_pos_t,  label = OR_text),  hjust = 0, vjust = 0.5, size = 3.2) +
    geom_text(aes(x = ci_pos_t,  label = CI_text),  hjust = 0, vjust = 0.5, size = 3.2) +
    geom_text(aes(x = p_pos_t,   label = P_text),   hjust = 0, vjust = 0.5, size = 3.2) +
    geom_text(aes(x = snp_pos_t, label = SNP_text), hjust = 0, vjust = 0.5, size = 3.2) +
    
    # Column headers (use the precomputed header y)
    geom_text(data = transform(header_df, y = y_header),
              inherit.aes = FALSE,
              aes(x = x, y = y, label = lab),
              hjust = c(1, 0, 0, 0, 0), vjust = 0.5, size = 3.8, fontface = "bold") +
    
    # "Odds Ratio" caption under the forest (centered in transformed space)
    annotate("text", x = (x_min_t + x_max_t) / 2, y = -0.3,
             label = "Odds Ratio", hjust = 0.5, vjust = 0.5, size = 3.8, fontface = "bold") +
    
    # Scales & theme (x-scale in transformed units; labels shown as OR)
    scale_x_continuous(
      limits = c(x_min_t - 0.8, x_limit_right),
      expand = c(0, 0),
      breaks = log2(tick_vals),
      labels = sprintf("%g", tick_vals)
    ) +
    scale_y_continuous(limits = c(-0.6, y_header + 1), expand = c(0, 0)) +
    scale_color_manual(values = c(
      "p < 0.001" = "#d73027",
      "p < 0.01"  = "#fc8d59",
      "p < 0.05"  = "#fee08b",
      "p ≥ 0.05"  = "#999999"
    ), name = "Significance") +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position  = "bottom",
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 9),
      legend.margin    = margin(t = 10),
      plot.margin      = margin(20, 20, 20, 20),
      plot.title       = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 20)),
      plot.subtitle    = element_text(hjust = 0.5, size = 11, margin = margin(b = 15))
    ) +
    labs(
      title = "Association of Genetically Predicted Endometriosis with Pregnancy Outcomes",
      subtitle = "Two-sample Mendelian randomization (Inverse variance weighted method)"
    )
  
  # Save endoPAIN layout forest plot
  ggsave(file.path(plots_dir, "forest_endoMR-PREG_endoPAIN_layout.png"),
         plot = p, width = 16, height = final_height, dpi = 300, bg = "white")
  message("Saved: ", file.path(plots_dir, "forest_endoMR-PREG_endoPAIN_layout.png"))
}


message("\n=== FOREST PLOTS COMPLETED (including endoPAIN layout) ===")
message("Generated forest plots:")
message("1. IVW forest plot: forest_endoMR-PREG_IVW.png")
if (nrow(df_multi) > 0) {
  message("2. Multi-method forest plot: forest_endoMR-PREG_multiMethod.png")
}
if (nrow(df_improved_style) > 0) {
  message("3. Improved forest plot: forest_endoMR-PREG_improved.png")
}
if (nrow(df_endopain_layout) > 0) {
  message("4. endoPAIN layout forest plot: forest_endoMR-PREG_endoPAIN_layout.png")
}
message("All plots saved to: ", plots_dir)