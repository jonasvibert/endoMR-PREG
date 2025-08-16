##################
# Visualizations #
##################

## A) Bar plot: SNP count per chromosome
chr_summary <- expdat_ready %>%
  filter(SNP %in% clumped$SNP) %>%    # keep only the 41 clumped SNPs
  mutate(chr = sub(":.*", "", chrpos)) %>%
  mutate(chr = factor(chr, levels = paste0("chr", 1:22))) %>%
  count(chr) %>%
  filter(n > 0)

ggplot(chr_summary, aes(x = chr, y = n)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Number of SNPs per Chromosome",
    x = "Chromosome",
    y = "Count of SNPs"
  ) +
  theme_minimal()

## B) Histogram of effect allele frequency (EAF)
ggplot(
  clumped2 %>% filter(SNP %in% clumped$SNP),
  aes(x = eaf.exposure)
) +
  geom_histogram(bins = 20, fill = "lightgray", color = "black") +
  labs(
    title = "Distribution of Effect Allele Frequency (EAF) for Clumped SNPs",
    x     = "EAF",
    y     = "Number of SNPs"
  ) +
  theme_minimal()

## C) Forest plot of β estimates with 95% CI
forest_df <- clumped2 %>%
  mutate(
    lower = beta.exposure - 1.96 * se.exposure,
    upper = beta.exposure + 1.96 * se.exposure
  ) %>%
  arrange(beta.exposure) %>%
  mutate(SNP = factor(SNP, levels = SNP))

ggplot(forest_df, aes(x = beta.exposure, y = SNP)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
  labs(
    title = "Forest Plot of SNP Effect Estimates",
    x = expression(beta~"(95% CI)"),
    y = "SNP"
  ) +
  theme_minimal()

## D) Volcano plot: β vs –log10(p‐value)
volc_df <- clumped2 %>%
  mutate(logp = -log10(pval.exposure.x))

ggplot(volc_df, aes(x = beta.exposure, y = logp)) +
  geom_point() +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
  labs(
    title = "Volcano Plot of SNP Effects",
    x = expression(beta),
    y = expression(-log[10]~"(p-value)")
  ) +
  theme_minimal()

## E) Scatterplot: F_stat vs R²
ggplot(clumped2, aes(x = R2, y = F_simple)) +
  geom_point() +
  labs(
    title = expression("F"["simple"]~"vs"~R^2),
    x = expression(R^2),
    y = expression(F["simple"])
  ) +
  theme_minimal()


##################
# TABLES
##################
library(dplyr)
library(openxlsx)

# Helper: coalesce over whatever columns exist
coalesce_into <- function(df, new_name, candidates) {
  df %>%
    mutate(
      "{new_name}" := coalesce(!!!select(., any_of(candidates)))
    )
}

# 1) Build unified columns safely
clumped2 <- clumped2 %>%
  coalesce_into("beta.exposure",       c("beta.exposure",       "beta.exposure.y",       "beta.exposure.x")) %>%
  coalesce_into("se.exposure",         c("se.exposure",         "se.exposure.y",         "se.exposure.x")) %>%
  coalesce_into("eaf.exposure",        c("eaf.exposure",        "eaf.exposure.y",        "eaf.exposure.x")) %>%
  coalesce_into("pval.exposure",       c("pval.exposure",       "pval.exposure.y",       "pval.exposure.x")) %>%
  coalesce_into("samplesize.exposure", c("samplesize.exposure", "samplesize.exposure.y", "samplesize.exposure.x"))

# 2) Ensure per-SNP R² and F are present
if (!all(c("R2_i","F_i") %in% names(clumped2))) {
  clumped2 <- clumped2 %>%
    mutate(
      R2_i = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
      F_i  = R2_i * (samplesize.exposure - 2) / (1 - R2_i)
    )
}

# 3) Build Supplementary Table S1 and format p-values as scientific notation
supp_table <- clumped2 %>%
  select(
    SNP,
    Beta = beta.exposure,
    SE = se.exposure,
    EAF = eaf.exposure,
    `P-value` = pval.exposure,
    `R²` = R2_i,
    `F-statistic` = F_i
  ) %>%
  mutate(`P-value` = formatC(`P-value`, format = "e", digits = 2))

# 4) Export to Excel
write.xlsx(supp_table, "Supplementary_Table_S1_Endometriosis_Instruments.xlsx", rowNames = FALSE)
cat("Saved: Supplementary_Table_S1_Endometriosis_Instruments.xlsx\n")
