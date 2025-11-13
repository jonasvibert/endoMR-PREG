# Methodology

## Mendelian Randomization Analysis: Endometriosis and Pregnancy Outcomes

### Study Design

This study employs a two-sample Mendelian randomization (MR) approach to investigate causal relationships between genetic liability to endometriosis and pregnancy outcomes.

### Exposure Data

**Source**: Rahmioglu et al. (2023) Nature Genetics
- **Sample size**: ~58,000 cases and 733,000 controls
- **Population**: European ancestry
- **Instruments**: 41 independent genome-wide significant SNPs (P < 5×10⁻⁸)
- **Clumping**: r² < 0.001, 10,000 kb window using 1000 Genomes EUR reference
- **Instrument strength**: Mean F-statistic = 279, Multi-SNP F = 296
- **Variance explained**: ~5.6% of endometriosis liability

### Outcome Data

**29 specialized pregnancy and maternal outcomes analyzed** (focused on pregnancy/delivery complications, maternal health, and fetal development):

**Data sources**:
- **MR-PREG consortium**: Up to 934,566 women
- **FinnGen R12**: ~176,899 participants (Finnish ancestry)
- **Westergaard et al.**: ~5,000 PPH cases, ~170,000 controls

**Outcome categories** (29 specialized outcomes):
1. **Placental disorders**: Placenta praevia, premature placental separation, placental disorders
2. **Delivery complications**: Caesarean section (elective/emergency/any), labor induction
3. **Hypertensive disorders**: Preeclampsia, gestational hypertension, hypertensive disorders of pregnancy
4. **Metabolic**: Gestational diabetes mellitus, pregnancy anemia
5. **Preterm birth**: Any preterm birth, very preterm birth, premature rupture of membranes
6. **Birth outcomes**: High/low birthweight, small for gestational age, gestational age, post-term birth
7. **Neonatal**: Low Apgar scores (1 min/5 min), NICU admission, stillbirth
8. **Bleeding**: Postpartum hemorrhage (overall), PPH due to atony, PPH due to retained placenta, antepartum bleeding
9. **Maternal health**: Postpartum depression

### Statistical Methods

#### Primary Analysis
- **Inverse Variance Weighted (IVW)**: Main method assuming all variants are valid instruments
- **MR-Egger regression**: Allows for horizontal pleiotropy with intercept test
- **Weighted median**: Robust to up to 50% invalid instruments
- **Mode-based methods**: Simple and weighted mode for robustness

#### 🆕 Fetal Genetic Effect Analysis (TRIOS)
- **TRIOS methodology**: Uses family trio data (mother-father-offspring) to distinguish maternal vs fetal genetic effects
- **Maternal genetic effects**: Using maternal genotypes as instruments
- **Fetal genetic effects**: Using offspring genotypes as instruments  
- **Paternal genetic effects**: Using paternal genotypes as negative controls
- **Pathway identification**: Determines whether pregnancy outcomes are driven by maternal or fetal genetic pathways
- **Validation**: Paternal effects serve as negative controls to confirm pathway specificity

#### Sensitivity Analyses
1. **Heterogeneity assessment**: Cochran's Q test
2. **Pleiotropy testing**: MR-Egger intercept test
3. **Leave-one-out analysis**: Influence of individual SNPs
4. **Funnel plots**: Visual assessment of directional pleiotropy

#### Quality Control
- **Instrument strength**: F-statistic > 10 (all instruments met this threshold)
- **Independence**: LD clumping r² < 0.001, 10,000 kb window
- **Palindromic SNPs**: Excluded if MAF > 0.42 to avoid strand ambiguity
- **Harmonization**: Effect alleles aligned between exposure and outcome
- **Minimum variants**: ≥3 SNPs required for MR-Egger analysis
- **🆕 Maternal-fetal distinction**: TRIOS analysis applied to identify maternal vs fetal genetic pathways
- **TRIOS quality control**: Paternal genetic effects as negative controls, leave-one-out sensitivity

### Software

- **R version**: 4.3+
- **Primary package**: TwoSampleMR (version 0.5.7+)
- **Additional packages**: dplyr, ggplot2, data.table

### Statistical Significance

- **FDR threshold**: q < 0.05 (False Discovery Rate correction for 29 specialized pregnancy outcomes)
- **Suggestive evidence**: 0.05 < q < 0.2 (hypothesis-generating, requires replication)
- **Effect sizes**: Odds ratios (OR) with 95% confidence intervals per unit increase in genetic liability
- **Two-sided testing**: All P-values are two-sided

### Assumptions

MR analysis relies on three core assumptions:
1. **Relevance**: Genetic variants strongly associated with endometriosis (F > 10 ✓)
2. **Independence**: Genetic variants independent of confounders (tested via MR-Egger)
3. **Exclusion restriction**: Variants affect outcomes only through endometriosis (assessed via pleiotropy tests)

**Validation approaches**:
- Strong instruments (mean F = 279) satisfy relevance
- MR-Egger intercept tests assess independence
- Multiple methods provide robustness to violations

### Limitations

- **Population**: Limited to European ancestry (generalizability)
- **Disease classification**: Cannot distinguish endometriosis from adenomyosis
- **Heterogeneity**: Cannot assess by endometriosis subtype or severity
- **Linearity**: Assumes linear dose-response relationships
- **Timing**: Reflects lifetime genetic liability, not timing-specific effects
- **Power**: Some true small effects may be missed
- **Pleiotropy**: Possible horizontal pleiotropy for some variants
- **🆕 TRIOS limitations**: Requires family trio data availability, assumption of linear fetal genetic effects
- **Pathway interpretation**: Binary maternal vs fetal classification may oversimplify complex gene-environment interactions

---

*For detailed statistical code, see `/scripts/` directory*