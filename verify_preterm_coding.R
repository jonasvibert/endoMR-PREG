#!/usr/bin/env Rscript
###############################################################################
# Vérification du codage des variables preterm birth
# Analyse de la cohérence clinique et directionnelle
###############################################################################

library(dplyr)
library(readr)

cat("=== VÉRIFICATION DU CODAGE PRETERM BIRTH ===\n\n")

# Charger les données harmonisées
harm_data <- read_csv("results/harmonised_rahmioglu_bpo.csv", show_col_types = FALSE)

# Charger les résultats MR
mr_results <- read_csv("results/mr_results_30_outcomes_fdr.csv", show_col_types = FALSE)

# 1. ANALYSE DES RÉSULTATS MR GLOBAUX
cat("1. RÉSULTATS MR POUR LES ISSUES DE PRÉMATURITÉ:\n")
cat("================================================\n")

# Filtrer les résultats preterm
preterm_outcomes <- mr_results %>%
  filter(grepl("pretb|preterm", outcome, ignore.case = TRUE)) %>%
  select(outcome, outcome_clean, OR, OR_lower, OR_upper, pval, b, se) %>%
  arrange(pval)

for (i in 1:nrow(preterm_outcomes)) {
  row <- preterm_outcomes[i,]
  direction <- ifelse(row$OR > 1, "RISQUE AUGMENTÉ", "EFFET PROTECTEUR")
  cat(sprintf("%s:\n", row$outcome_clean))
  cat(sprintf("  - OR = %.2f (%.2f - %.2f)\n", row$OR, row$OR_lower, row$OR_upper))
  cat(sprintf("  - Beta = %.4f ± %.4f\n", row$b, row$se))
  cat(sprintf("  - P = %.6f\n", row$pval))
  cat(sprintf("  - Direction: %s\n", direction))
  cat("\n")
}

# 2. ANALYSE DES DONNÉES INDIVIDUELLES POUR PRETB_ALL
cat("2. ANALYSE SNP-LEVEL POUR PRETERM BIRTH (pretb_all):\n")
cat("==================================================\n")

# Filtrer les données pour pretb_all
pretb_data <- harm_data %>%
  filter(outcome == "pretb_all") %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, 
         beta.exposure, beta.outcome, se.exposure, se.outcome, 
         pval.exposure, pval.outcome) %>%
  mutate(
    # Direction de l'effet pour endométriose
    endo_direction = ifelse(beta.exposure > 0, "↑ Endometriosis risk", "↓ Endometriosis risk"),
    # Direction de l'effet pour preterm birth  
    pretb_direction = ifelse(beta.outcome > 0, "↑ Preterm birth risk", "↓ Preterm birth risk"),
    # Concordance des directions
    concordant = sign(beta.exposure) == sign(beta.outcome)
  )

# Afficher les premiers SNPs
cat("Premiers 10 SNPs instrumentaux:\n")
head_data <- pretb_data %>% slice_head(n = 10)
for (i in 1:nrow(head_data)) {
  row <- head_data[i,]
  concordance <- ifelse(row$concordant, "CONCORDANT", "DISCORDANT")
  cat(sprintf("%s: beta_endo=%.4f, beta_pretb=%.4f - %s\n", 
              row$SNP, row$beta.exposure, row$beta.outcome, concordance))
}

cat(sprintf("\nConcordance globale: %d/%d SNPs concordants (%.1f%%)\n",
            sum(pretb_data$concordant), nrow(pretb_data), 
            100 * sum(pretb_data$concordant) / nrow(pretb_data)))

# 3. ANALYSE DE COHÉRENCE CLINIQUE
cat("\n3. ANALYSE DE COHÉRENCE CLINIQUE:\n")
cat("=================================\n")

# Comparer avec les résultats cohérents
coherent_results <- mr_results %>%
  filter(outcome %in% c("finngen_R12_O15_PLAC_PRAEVIA", "el_cs", "rup_memb")) %>%
  select(outcome_clean, OR, pval, b) %>%
  mutate(direction = ifelse(OR > 1, "Augmente risque", "Diminue risque"))

cat("Résultats cliniquement cohérents:\n")
for (i in 1:nrow(coherent_results)) {
  row <- coherent_results[i,]
  cat(sprintf("- %s: OR=%.2f, %s\n", row$outcome_clean, row$OR, row$direction))
}

cat("\nLogique clinique attendue:\n")
cat("- Placenta praevia → césariennes précoces → PLUS de prématurité\n")
cat("- PROM → accouchement prématuré → PLUS de prématurité\n") 
cat("- Endométriose → complications → PLUS de prématurité\n")

cat("\nRésultat observé pour preterm birth:\n")
pretb_result <- mr_results %>% filter(outcome == "pretb_all")
cat(sprintf("- Preterm birth: OR=%.2f → MOINS de prématurité\n", pretb_result$OR))
cat("- ⚠️  INCOHÉRENCE CLINIQUE MAJEURE ⚠️\n")

# 4. HYPOTHÈSES EXPLICATIVES
cat("\n4. HYPOTHÈSES EXPLICATIVES POSSIBLES:\n")
cat("====================================\n")

cat("A) PROBLÈME DE CODAGE:\n")
cat("   - Variable 'preterm birth' pourrait être codée à l'envers\n")
cat("   - 1 = term birth, 0 = preterm birth (inverse de l'attendu)\n")
cat("   - Vérifier la documentation originale des GWAS\n\n")

cat("B) PROBLÈME DE SÉLECTION:\n") 
cat("   - Biais de survie (prématurés sévères exclus)\n")
cat("   - Population spécifique (femmes avec grossesses réussies)\n\n")

cat("C) CONFONDANT RÉSIDUEL:\n")
cat("   - Surveillance accrue → détection précoce → interventions\n")
cat("   - Endométriose → suivi médical → prévention prématurité\n\n")

cat("D) ERREUR D'ANALYSE:\n")
cat("   - Mauvaise harmonisation des allèles\n")
cat("   - Problème dans le calcul MR\n")

# 5. RECOMMANDATIONS
cat("\n5. RECOMMANDATIONS POUR VÉRIFICATION:\n")
cat("====================================\n")

cat("1. VÉRIFIER LE CODAGE ORIGINAL:\n")
cat("   - Consulter la publication EGG consortium\n")
cat("   - Vérifier si 1 = preterm ou 1 = term\n")
cat("   - Confirmer la direction des betas\n\n")

cat("2. ANALYSE SENSIBILITÉ:\n")
cat("   - Inverser le codage et refaire l'analyse\n")
cat("   - Comparer avec d'autres outcomes temporels\n")
cat("   - Vérifier la cohérence avec gestational age\n\n")

cat("3. VALIDATION EXTERNE:\n")
cat("   - Comparer avec littérature observationnelle\n")
cat("   - Vérifier d'autres cohortes (UK Biobank)\n")
cat("   - Contrôler avec outcomes similaires\n")

cat("\n=== CONCLUSION ===\n")
cat("Le résultat 'protecteur' pour preterm birth est CLINIQUEMENT INCOHÉRENT.\n")
cat("Il y a très probablement un problème de codage ou d'interprétation.\n")
cat("RECOMMANDATION: Vérifier le codage original avant publication.\n")