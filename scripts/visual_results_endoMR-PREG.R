# 1. Extract and prepare SNP‐specific results for female infertility
data <- res_single %>%
  filter(id.outcome == "finngen_R12_N14_FEMALEINFERT_filtered") %>%
  mutate(
    OR    = exp(b),
    lower = exp(b - 1.96 * se),
    upper = exp(b + 1.96 * se)
  ) %>%
  arrange(OR) %>%
  mutate(SNP = factor(SNP, levels = SNP))

# 2. Forest plot
library(ggplot2)
ggplot(data, aes(x = SNP, y = OR, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "Forest plot of SNP effects on Female Infertility (FinnGen)",
    x = "SNP",
    y = "Odds ratio (95% CI)"
  ) +
  theme_minimal(base_size = 9)
