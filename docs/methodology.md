# Methodology

## Mendelian Randomization Analysis: Endometriosis and Pregnancy Outcomes

### Study Design

This study employs a two-sample Mendelian randomization (MR) approach to investigate causal relationships between genetic liability to endometriosis and pregnancy outcomes.

### Exposure Data

**Source**: Endometriosis GWAS
- **Sample size**: [Add sample size]
- **Population**: European ancestry
- **Instruments**: 41 genome-wide significant SNPs (P < 5×10⁻⁸)
- **Clumping**: r² < 0.001, 10,000 kb window using 1000 Genomes EUR reference

### Outcome Data

**Primary outcomes analyzed**:

1. **Fertility outcomes**
   - Female infertility (FinnGen R12)

2. **Placental disorders**
   - Placenta praevia (FinnGen R12)
   - Placental abruption (FinnGen R12)

3. **Preterm birth and complications**
   - Premature rupture of membranes
   - Preterm birth (various definitions)

4. **Delivery complications**
   - Caesarean section (elective/emergency)
   - Postpartum hemorrhage (overall and specific causes)

5. **Fetal growth and neonatal outcomes**
   - Birthweight categories
   - Apgar scores
   - NICU admission

### Statistical Methods

#### Primary Analysis
- **Inverse Variance Weighted (IVW)**: Main method assuming all variants are valid instruments
- **MR-Egger regression**: Allows for horizontal pleiotropy with intercept test
- **Weighted median**: Robust to up to 50% invalid instruments
- **Mode-based methods**: Simple and weighted mode for robustness

#### Sensitivity Analyses
1. **Heterogeneity assessment**: Cochran's Q test
2. **Pleiotropy testing**: MR-Egger intercept test
3. **Leave-one-out analysis**: Influence of individual SNPs
4. **Funnel plots**: Visual assessment of directional pleiotropy

#### Quality Control
- Minimum 3 SNPs required for MR-Egger analysis
- F-statistic > 10 for instrument strength
- Palindromic SNPs with intermediate allele frequencies excluded
- Harmonization of effect alleles between exposure and outcome

### Software

- **R version**: 4.3+
- **Primary package**: TwoSampleMR (version 0.5.7+)
- **Additional packages**: dplyr, ggplot2, data.table

### Statistical Significance

- **Primary outcomes**: P < 0.05 for IVW method
- **Multiple testing correction**: Bonferroni correction applied where appropriate
- **Effect sizes**: Reported as odds ratios (OR) with 95% confidence intervals

### Assumptions

MR analysis relies on three core assumptions:
1. **Relevance**: Genetic variants strongly associated with exposure
2. **Independence**: Genetic variants independent of confounders
3. **Exclusion restriction**: Genetic variants affect outcome only through exposure

### Limitations

- Limited to European ancestry populations
- Potential for population stratification
- Cannot detect non-linear causal relationships
- Assumes lifetime exposure effect

---

*For detailed statistical code, see `/scripts/` directory*