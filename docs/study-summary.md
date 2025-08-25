# Study Summary: Association of Genetic Liability to Endometriosis With Pregnancy Outcomes

## Authors
Jonas Vibert, MD¹, Carolina Borges, PhD², Zoltán Kutalik, PhD³, David Baud, MD, PhD¹, Deborah A. Lawlor, PhD², Nicola Pluchino, MD, PhD¹

¹Department of Obstetrics and Gynecology, Lausanne University Hospital, Lausanne, Switzerland  
²MRC Integrative Epidemiology Unit at the University of Bristol, Bristol, UK  
³University Center for Primary Care and Public Health, Lausanne, Switzerland

## Study Design
**Two-sample Mendelian randomization study** investigating causal relationships between genetic liability to endometriosis and 49 maternal and fetal pregnancy outcomes.

## Key Methodology Details

### Exposure
- **Source**: Rahmioglu et al. (2023) Nature Genetics GWAS
- **Sample**: ~58,000 endometriosis cases, 733,000 controls (European ancestry)
- **Instruments**: 41 independent SNPs (genome-wide significant, P < 5×10⁻⁸)
- **Instrument strength**: Mean F-statistic = 279, multi-SNP F = 296
- **Explained variance**: ~5.6% of endometriosis liability

### Outcomes (49 total)
- **MR-PREG consortium**: Up to 678,001 women across multiple cohorts
- **FinnGen R12**: 176,899 participants (predominantly Finnish ancestry)
- **Westergaard et al.**: ~5,000 PPH cases, ~170,000 controls (European ancestry)

### Statistical Analysis
- **Primary method**: Inverse Variance Weighted (IVW)
- **Sensitivity analyses**: MR-Egger, Weighted Median, MR-PRESSO
- **Quality control**: Heterogeneity (Cochran's Q), pleiotropy testing, leave-one-out
- **Multiple testing correction**: Bonferroni (P < 0.001 for significance)
- **Software**: R v4.3.2, TwoSampleMR v0.5.6, MR-PRESSO v1.0.0

## Main Results

### Bonferroni-Corrected Significant Associations
1. **Placenta praevia**: OR 1.62 (95% CI: 1.33–1.97), P = 1.5×10⁻⁶
2. **Female infertility**: OR 1.63 (95% CI: 1.48–1.80), P = 6.1×10⁻²² *(internal validation)*

### Nominally Significant Associations (P < 0.05)
- Premature rupture of membranes: OR 1.12 (1.01–1.23), P = 0.025
- Elective cesarean delivery: OR 1.26 (1.04–1.53), P = 0.023
- Placental abruption: OR 1.36 (1.03–1.81), P = 0.031
- 1-minute Apgar score: OR 0.95 (0.91–0.99), P = 0.029

### No Significant Associations
- Preeclampsia
- Gestational diabetes
- Postpartum hemorrhage
- Stillbirth
- Fetal growth restriction
- Most other pregnancy complications

## Clinical Interpretation

### Primary Finding
**Placenta praevia** emerged as the only pregnancy outcome with strong causal evidence, supporting a 62% increased risk with genetic liability to endometriosis. This aligns with biological mechanisms involving:
- Impaired decidualization
- Defective spiral artery remodeling
- Altered trophoblast invasion
- Junctional zone abnormalities

### Clinical Implications
- Most previously reported associations between endometriosis and adverse pregnancy outcomes likely reflect **bias rather than causality**
- Potential sources of bias include:
  - Assisted reproductive technology use
  - Increased cesarean delivery rates
  - Maternal age effects
  - Intensified obstetric surveillance
  - Clinical misclassification (endometriosis vs. adenomyosis)

### Study Limitations
- Restricted to European ancestry populations
- Unable to assess heterogeneity by endometriosis subtype or severity
- Potential misclassification between endometriosis and adenomyosis
- Cannot exclude small effect sizes below detection threshold

## Conclusions

This comprehensive Mendelian randomization study provides evidence for a **causal association between genetic liability to endometriosis and placenta praevia**, but not with most other adverse pregnancy outcomes. These findings suggest that the clinical risks attributed to endometriosis in observational studies may be largely due to confounding factors rather than direct biological effects of the condition itself.

## Data Availability
- Complete analysis code: [GitHub repository](https://github.com/jonasvibert/endoMR-PREG)
- Summary statistics: Available from original GWAS sources
- Individual-level data: Not available (summary statistics only)