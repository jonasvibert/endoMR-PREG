###############################################
# Table 3. Primary Mendelian Randomization (IVW) Estimates
# for the Association of Genetically Predicted Endometriosis
# With Maternal and Fetal Pregnancy Outcomes
###############################################

library(dplyr)
library(openxlsx)

# Dictionnaire de noms complets pour les outcomes
outcome_labels <- c(
  pretb_all      = "Preterm birth (all)",
  el_cs          = "Elective caesarean section",
  rup_memb       = "Premature rupture of membranes",
  pretb_subsamp  = "Preterm birth (subsample)",
  apgar1         = "Apgar score at 1 min",
  vpretb_all     = "Very preterm birth (all)",
  em_cs          = "Emergency caesarean section",
  lowapgar1      = "Low Apgar score at 1 min",
  pe_subsamp     = "Preeclampsia (subsample)",
  lbw_all        = "Low birthweight (<2500 g)",
  ga_all         = "Gestational age (all)",
  gh_subsamp     = "Gestational hypertension (subsample)",
  ga_subsamp     = "Gestational age (subsample)",
  lga            = "Large for gestational age",
  zbw_all        = "Z-score birthweight (all)",
  nvp_sev_subsamp= "Severe nausea/vomiting (subsample)",
  nvp_sev_all    = "Severe nausea/vomiting (all)",
  sb_subsamp     = "Stillbirth (subsample)",
  gdm_subsamp    = "Gestational diabetes (subsample)",
  apgar5         = "Apgar score at 5 min",
  hdp_subsamp    = "Hypertensive disorders of pregnancy (subsample)",
  r_misc_subsamp = "Recurrent miscarriage (subsample)",
  bf_dur_4c      = "Breastfeeding ≥4 months",
  depr_subsamp   = "Postpartum depression (subsample)",
  hbw_all        = "High birthweight (>4000 g)",
  nicu           = "NICU admission",
  lowapgar5      = "Low Apgar score at 5 min",
  bf_sus         = "Breastfeeding cessation",
  posttb_all     = "Post-term birth",
  misc_subsamp   = "Miscarriage (subsample)",
  bf_ini         = "Breastfeeding initiation",
  bf_est         = "Exclusive breastfeeding",
  s_misc_subsamp = "Single miscarriage (subsample)",
  cs             = "Caesarean section (all)",
  hyp            = "Hyperemesis gravidarum",
  anaemia_preg_all = "Anaemia in pregnancy (all)",
  sga            = "Small for gestational age",
  induction      = "Induction of labour",
  finngen_R12_N14_FEMALEINFERT_filtered = "Female infertility",
  finngen_R12_O15_PLAC_PRAEVIA_filtered = "Placenta praevia",
  finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered = "Placental abruption",
  finngen_R12_O15_PLAC_DISORD_filtered = "Other placental disorders",
  finngen_R12_O15_PREG_ECTOP_filtered   = "Ectopic pregnancy",
  Early_bleeding_with_any_outcome_filtered = "Early bleeding (any outcome)",
  Postpartum_hemorrhage_due_to_retained_placenta_filtered = "PPH – retained placenta",
  Early_bleeding_ending_in_live_birth_filtered = "Early bleeding – live birth",
  Postpartum_hemorrhage_filtered = "Postpartum haemorrhage",
  Postpartum_hemorrhage_due_to_atony_filtered = "PPH – atony",
  Antepartum_bleeding_filtered = "Antepartum bleeding"
)

# Remplacer les codes par les noms complets
ivw_tab <- ivw_tab %>%
  mutate(Outcome = recode(Outcome, !!!outcome_labels))

# Renommer les colonnes pour plus de clarté
ivw_tab_pub <- ivw_tab %>%
  rename(
    `Data source`    = Source,
    `Outcome`        = Outcome,
    `No. SNPs`       = SNPs,
    `OR`             = OR_fmt,
    `95% CI`         = CI_fmt,
    `P`              = p_fmt,
    `P (Bonferroni)` = pBonf_fmt,
    `Q (FDR)`        = qFDR_fmt,
    `Significance`   = Signif
  )

# Afficher dans la console
print(ivw_tab_pub, n = Inf, width = Inf)

# Exporter vers Excel
write.xlsx(ivw_tab_pub, "ivw_primary_summary_named.xlsx", asTable = TRUE)
