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
  mutate(logp = -log10(pval.exposure))

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
