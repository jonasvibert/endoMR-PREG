# 📊 Supplementary Methods

This two-sample Mendelian randomization (MR) study was designed to evaluate whether genetic liability to endometriosis has a causal effect on a range of adverse pregnancy and perinatal outcomes. 

## 🔬 Analysis Pipeline

The analysis pipeline was implemented in R using a set of modular scripts stored in the [endoMR-PREG GitHub repository](https://github.com/jonasvibert/endoMR-PREG/tree/3c6357c4e89748dc67e96dd292d10615feffbb93/scripts):

```
📁 scripts/
├── 01_select_instruments_endoMR-PREG.R
├── 02_prepare_outcomes_endoMR-PREG.R
├── 03_harmonise_data_endoMR-PREG.R
├── 04_main_analyses_endoMR-PREG.R
├── 04.2_fetal_effect_endoMR-PREG.R
├── 05_sensitivity_analyses_endoMR-PREG.R
├── 06_tables_endoMR-PREG.R
└── 07_plots_endoMR-PREG.R
```

Together, these scripts reproduce the full workflow, from instrument selection to figure generation.

## 🧬 Exposure: Genetic Liability to Endometriosis

Genetic instruments for endometriosis were derived from the largest published genome-wide association study (GWAS) meta-analysis of clinically confirmed endometriosis by **Rahmioglu et al.** *(Nature Genetics 2023)* <sup>[1]</sup>. 

### 📈 Study Characteristics
- **Cases**: ~58,000 surgically and/or clinically diagnosed
- **Controls**: ~733,000 female controls  
- **Ancestry**: Predominantly European
- **Contributing cohorts**: International Endometriosis Genetics Consortium, 23andMe, FinnGen and others

In `01_select_instruments_endoMR-PREG.R`, we imported the European-ancestry summary statistics from the meta-analysis that included 23andMe, and restricted all subsequent analyses to these data to maintain ancestry compatibility with the outcome GWAS.

### 🎯 Instrument Selection Procedure

#### Step 1: Variant Extraction
We extracted all variants reaching conventional genome-wide significance for endometriosis:
- **Threshold**: P < 5 × 10⁻⁸
- **Data extracted**: rsIDs, effect/non-effect alleles, beta, SE, EAF, sample size (N)
- **Coordinate mapping**: chr:pos format variants mapped to rsIDs using 1000 Genomes Project Phase 3 European reference panel <sup>[2]</sup>
- **Quality control**: Variants with ambiguous mapping or unresolved alleles were excluded

#### Step 2: Linkage Disequilibrium (LD) Clumping
To ensure independence between instruments, we performed LD clumping using **PLINK 1.9**:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `clump_r2` | 0.001 | LD threshold |
| `clump_kb` | 10,000 | Physical window (kb) |
| `clump_p1` | 5 × 10⁻⁸ | Primary P-value threshold |
| `clump_p2` | 1 | Secondary threshold |
| **Reference** | 1000 Genomes Phase 3 Europeans | Population reference |

> **Result**: 41 independent SNPs retained as instruments for genetic liability to endometriosis

### 💪 Instrument Strength Evaluation

Instrument strength was systematically evaluated using multiple metrics:

#### F-Statistics
**Simple F-statistic** for each SNP:
```
F_simple = (beta / SE)²
```

**Variance explained** in endometriosis liability (R²):
```
R² = [ 2 × EAF × (1 − EAF) × beta² ] / [ 2 × EAF × (1 − EAF) × beta² + SE² × N ]
```

**Exact per-SNP F-statistic**:
```
F_exact = R² × (N − 2) / (1 − R²)
```

#### 📊 Strength Assessment Results
- **Total variance explained**: ~5.6% (summing R² across all independent SNPs)
- **Mean per-SNP F-statistic**: ~279
- **Weak instrument risk**: Very low (F >> 10 threshold)

> 📋 Full instrument characteristics (rsID, genomic position, alleles, beta, SE, EAF, R², F) are provided in **Supplementary Table S1**.

## 🎯 Outcome GWAS and Phenotype Definition

We evaluated **29 maternal and perinatal outcomes** selected based on:
- 📚 Prior epidemiological evidence linking endometriosis with obstetric complications
- 🏥 Clinical relevance for maternal–fetal health  
- 📊 Availability of adequately powered GWAS in European-ancestry populations

### 📋 Data Sources
We drew outcome summary statistics from three main sources:

| Source | Description | Reference |
|--------|-------------|-----------|
| **MR-PREG** | Mendelian Randomization in Pregnancy collaboration | <sup>[3]</sup> |
| **FinnGen** | Release 12 | <sup>[4]</sup> |
| **PPH meta-analysis** | Large postpartum haemorrhage study | <sup>[5]</sup> |

> 🔧 Extraction, cleaning and initial formatting performed in `02_prepare_outcomes_endoMR-PREG.R`
### 🤝 MR-PREG Collaboration
**MR-PREG** is an international consortium that harmonises and meta-analyses GWAS of pregnancy and perinatal outcomes across several large European birth cohorts:

**Contributing cohorts**: ALSPAC, Born in Bradford, MoBa, UK Biobank

**Standardised pipeline** <sup>[3]</sup>:
- **Sample-level QC**: Removal of individuals with low call rate, sex discordance, outlying heterozygosity or cryptic relatedness
- **Ancestry restriction**: European genetic ancestry participants only
- **Variant-level QC**: Call rate, MAF, HWE, imputation quality filters  
- **Association analyses**: Logistic/linear regression or REGENIE with appropriate covariates
- **Meta-analysis**: Fixed-effects inverse-variance weighting

#### 🏥 MR-PREG Outcomes (21 total)
**From consortium's meta-analysed GWAS:**
- Gestational age (full sample and genotyped subsample)
- Gestational diabetes, gestational hypertension, hypertensive disorders of pregnancy
- High and low birthweight, small-for-gestational-age, birthweight z-score, large-for-gestational-age
- Preterm birth, very preterm birth, post-term birth
- Labour induction, premature rupture of membranes
- Apgar score <7 at 1 and 5 minutes, NICU admission
- Pregnancy anaemia, preeclampsia

**Additional outcomes** (same QC framework):
- Caesarean section phenotypes (overall, elective, emergency)
- Stillbirth, postpartum/peripartum depression
### 🇫🇮 FinnGen Release 12
To complement these outcomes with detailed **placental phenotypes**, we used FinnGen release 12 <sup>[4]</sup>.

**FinnGen framework** (Kurki et al., 2023):
- Array-based genotyping
- Imputation to SISu v3 Finnish reference panel  
- Strict sample and variant filtering
- Association testing using REGENIE
- Registry-based ICD code definitions

**FinnGen outcomes extracted:**
- 🩸 Placenta praevia
- 🔬 Placental disorders  
- ⚠️ Placental abruption

### 🩸 Bleeding-Related Outcomes
From **Westergaard et al. (2024)** meta-analysis <sup>[5]</sup> (Nordic biobanks + UK Biobank):
- Antepartum bleeding
- Postpartum haemorrhage (PPH) overall
- PPH due to uterine atony
- PPH due to retained placenta

> **Note**: Odds ratios converted to log-odds (β = log(OR)); missing SEs derived from P-values: SE = |β| / z

### 📊 Clinical Domain Organization
The **29 outcomes** are grouped into **9 clinical domains**:

| Domain | Outcomes |
|--------|----------|
| 🩸 **Placental disorders** | Placenta praevia, placental disorders, abruption |
| 🩺 **Hypertensive disorders** | Hypertensive disorders overall, gestational hypertension, preeclampsia |
| ⏰ **Pregnancy timing** | Gestational age, preterm birth, very preterm birth, post-term birth |
| 👶 **Fetal growth/birthweight** | High/low birthweight, SGA, birthweight z-score, LGA |
| 🏥 **Labour/delivery complications** | Labour induction, PROM, caesarean section (overall, elective, emergency) |
| 💉 **Bleeding/haemorrhage** | Antepartum bleeding, PPH (overall, atony, retained placenta) |
| 🔬 **Maternal metabolic/haematologic** | Gestational diabetes, pregnancy anaemia |
| 🧠 **Maternal mental health** | Postpartum/peripartum depression |
| 👶 **Neonatal condition** | Low Apgar scores (1 & 5 min), NICU admission, stillbirth |

> 📋 Detailed definitions, case/control counts and contributing cohorts: **Table 1** and **Supplementary Table S2**

## 🔗 Harmonisation of Exposure and Outcome Data

Harmonisation of endometriosis and outcome GWAS was carried out in `03_harmonise_data_endoMR-PREG.R` using:
- **Package**: TwoSampleMR (version 0.5.6)
- **Additional**: Manual quality checks

### 🎯 Harmonisation Procedure

#### Step 1: SNP Matching
- Match 41 endometriosis instrument SNPs to outcome SNPs by **rsID**

#### Step 2: Allele Alignment  
- Align alleles: identical effect alleles in exposure and outcome datasets
- **Beta flipping**: Adjust outcome beta sign when necessary
- **Frequency adjustment**: Modify allele frequencies as required

#### Step 3: Palindromic SNP Handling
| SNP Type | MAF Range | Action | Rationale |
|----------|-----------|--------|-----------|
| A/T or C/G | ~0.5 (intermediate) | **Removed** | Ambiguous strand assignment |
| A/T or C/G | Clearly low/high | **Retained** | Strand reliably inferred |

#### Step 4: Quality Control
**Excluded SNPs with:**
- Unresolved strand issues
- Mismatched alleles  
- Missing beta/SE in either dataset

### 📊 Final Instrument Count
- **Range**: 29-41 instruments per MR analysis
- **Reason**: Not all instrument SNPs survived QC or were available in every outcome GWAS

> 📋 Exact SNP counts per outcome: **Supplementary Table S2** (includes exposure and outcome sample sizes)

## 🧮 Two-Sample MR Analyses

Primary MR analyses were implemented in `04_main_analyses_endoMR-PREG.R`.

### 📏 Inverse-Variance Weighted (IVW) Estimator  
**Main causal effect measure** <sup>[6]</sup>

For each SNP *i*:
- β<sub>Xi</sub> = SNP–endometriosis association  
- β<sub>Yi</sub> = SNP–outcome association
- var<sub>Yi</sub> = variance of β<sub>Yi</sub>

**IVW estimator formula:**
```
β_IVW = Σ(β_Xi × β_Yi / var_Yi) / Σ(β_Xi² / var_Yi)
```

> This is equivalent to a weighted regression of β<sub>Yi</sub> on β<sub>Xi</sub> with no intercept and weights = 1/var<sub>Yi</sub>

### 📊 Effect Reporting
| Outcome Type | Metric | Formula | CI |
|--------------|--------|---------|-----|
| **Binary** | Odds Ratio | OR = exp(β<sub>IVW</sub>) | 95% CI |
| **Continuous** | Beta coefficient | β<sub>IVW</sub> | 95% CI |

*Continuous outcomes: gestational age, birthweight z-score*

Sensitivity analyses and assessment of MR assumptions

Sensitivity analyses were conducted in 05_sensitivity_analyses_endoMR-PREG.R to evaluate the robustness of IVW estimates to potential violations of MR assumptions, particularly horizontal pleiotropy. First, we performed MR-Egger regression, modelling beta_Yi as:

beta_Yi = intercept + beta_MR_Egger × beta_Xi + error

The slope beta_MR_Egger provides a pleiotropy-adjusted causal estimate under the Instrument Strength Independent of Direct Effect (InSIDE) assumption, while the intercept tests for directional pleiotropy: a non-zero intercept suggests that instruments have, on average, direct effects on the outcome independent of endometriosis liability.

Second, we estimated causal effects using the weighted median estimator, which remains consistent if at least 50% of the total instrument weight is contributed by valid instruments. This estimator is based on the median of SNP-specific ratio estimates (beta_Yi / beta_Xi), weighted by the inverse variance of beta_Yi. Where informative, we additionally computed weighted mode estimates, which assume that the most common (modal) causal estimate across SNPs arises from valid instruments.

Third, we quantified heterogeneity among SNP-specific causal estimates using Cochran’s Q statistic. For each SNP we calculated its ratio estimate beta_i = beta_Yi / beta_Xi and weight w_i = beta_Xi^2 / var_Yi, and then computed:

Q = sum( w_i × (beta_i − beta_IVW)^2 )

A large Q statistic with a low P-value indicates heterogeneity and may suggest pleiotropy or violation of model assumptions.

Fourth, we applied the MR-PRESSO (Pleiotropy RESidual Sum and Outlier) method (MRPRESSO package v1.0.0). For each outcome, the global test examined whether the pattern of residuals was compatible with no horizontal pleiotropy. When the global test was significant, MR-PRESSO identified outlier SNPs and recomputed outlier-corrected IVW estimates after excluding them. We report global test P-values, outlier SNPs where present and corrected effect estimates.

Finally, we performed single-SNP MR and leave-one-SNP-out analyses. Single-SNP estimates examine the contribution of each instrument individually, while leave-one-out analyses re-run IVW MR after removing each SNP in turn to evaluate whether any single variant unduly drives results. Figures summarising these diagnostics (scatter plots, forest plots, funnel plots, leave-one-out plots) were generated in 07_plots_endoMR-PREG.R and are presented in Supplementary Figures S1–S7.

Maternal, fetal and paternal genetic effects

Some outcomes, such as hypertensive disorders of pregnancy, birthweight and preterm birth, may be influenced by both maternal and fetal genomes. Where available, we therefore leveraged trio-based genetic effect estimates from MR-PREG to distinguish maternal and fetal pathways. These analyses were implemented in 04.2_fetal_effect_endoMR-PREG.R.

MR-PREG provides effect estimates from models jointly including maternal, offspring (fetal) and paternal genotypes, derived from mother–father–child trios or extended parent–offspring structures, using regression models of the form 

expected outcome = beta_m × G_m + beta_o × G_o + beta_p × G_p + covariates

, where GmG_mGm​, GoG_oGo​ and GpG_pGp​ represent maternal, offspring and paternal genotypes. Maternal effects (βm\beta_mβm​) capture direct maternal genetic influence on the intrauterine environment, fetal effects (βo\beta_oβo​) reflect the fetus’s own genome acting on fetal growth and development, and paternal effects (βp\beta_pβp​) serve as a negative control, as paternal genotype should not directly affect maternal pregnancy complications.

Trio-based genetic associations were available for 21 of the 29 investigated outcomes, namely small-for-gestational-age, large-for-gestational-age, high birthweight (>4000 g), low birthweight (<2500 g), birthweight z-score, gestational age, post-term birth, preterm birth <37 weeks, preeclampsia, gestational hypertension, hypertensive disorders of pregnancy, gestational diabetes, pregnancy anemia, labour induction, emergency caesarean section, elective caesarean section, low Apgar score at 1 minute, low Apgar score at 5 minutes, neonatal intensive care unit admission, postpartum depression and premature rupture of membranes. 
For each of these outcomes, we performed three parallel MR analyses using the endometriosis SNP instruments: one using maternal effect estimates, one using fetal effect estimates and one using paternal effect estimates as a negative control. The same causal estimators were applied in each pathway (IVW, MR-Egger, weighted median, weighted mode), together with heterogeneity statistics, MR-Egger intercepts, MR-PRESSO global and outlier tests and leave-one-out analyses. Because not all trio-based SNP associations were available for all outcomes, the number of instruments contributing to each maternal, fetal, and paternal analysis varied accordingly, ranging from 27 to 34 SNPs (details in Supplementary Table S3). Results of these analyses, allowing direct comparison of maternal, fetal and paternal pathways, are reported in Supplementary Tables S3–S4 and Supplementary Figure S8.

Multiple testing correction

Because we evaluated causal estimates for 29 outcomes, we controlled for multiplicity using the false discovery rate (FDR). In 05_sensitivity_analyses_endoMR-PREG.R, we applied the Benjamini–Hochberg procedure to the two-sided IVW P-values across all outcomes using R’s p.adjust function in R with method = "fdr". The resulting FDR-adjusted P-values (q-values) reflect the expected proportion of false positives among declared significant findings. We considered associations with q < 0.05 as statistically significant after FDR correction, and we report both nominal P-values and q-values in Supplementary Table S2.

Tables, figures and reproducibility

The script 06_tables_endoMR-PREG.R assembled results from all previous steps into the tables presented in the main manuscript and supplementary material. It joins MR estimates, confidence intervals, heterogeneity statistics, MR-Egger intercepts, MR-PRESSO outputs, SNP counts and sample sizes into publication-ready tables (for example, for IVW estimates, sensitivity analyses and trio-based MR). 07_plots_endoMR-PREG.R uses the same underlying results to generate forest plots of IVW estimates across outcomes (including the main 29-outcome figure), scatter plots contrasting SNP–exposure and SNP–outcome associations, funnel plots for asymmetry, leave-one-out plots and trio comparison plots.

All analyses were conducted in R version 4.3.2. Key packages included TwoSampleMR (v0.5.6) for harmonisation and MR methods, MRPRESSO (v1.0.0) for pleiotropy assessment, data.table and dplyr for data handling, ggplot2 for visualisation, and readr, stringr, tidyr and here for import and project organisation. LD clumping was performed using PLINK 1.9. The full set of scripts (01–07) and a detailed README are available in the publicly accessible endoMR-PREG repository (https://github.com/jonasvibert/endoMR-PREG
), allowing complete reproducibility of the analysis pipeline.

Data availability and ethics

This study relied exclusively on de-identified summary-level GWAS data from previously approved studies and consortia. Endometriosis GWAS summary statistics were obtained from Rahmioglu et al. (Nature Genetics 2023) (1). MR-PREG outcome GWAS were accessed through the MR-PREG collaboration under its data-sharing policies (3). FinnGen summary statistics (release 12) are publicly available through the FinnGen results portal (4), and bleeding outcomes were obtained from Westergaard et al. (Nature Genetics 2024) (5). Because only aggregated summary statistics were used and no individual-level genetic or clinical data were accessed, the present analyses do not constitute human subjects research as defined by 45 CFR 46.102, and no additional institutional review board approval was required. All original contributing studies obtained appropriate ethical approvals and informed consent from participants.

References
1.	Rahmioglu N, Mortlock S, Ghiasi M, Møller PL, Stefansdottir L, Galarneau G, et al. The genetic basis of endometriosis and comorbidity with other pain and inflammatory conditions. Nat Genet. 2023 Mar;55(3):423–36. 
2.	The 1000 Genomes Project Consortium, Corresponding authors, Auton A, Abecasis GR, Steering committee, Altshuler DM, et al. A global reference for human genetic variation. Nature. 2015 Oct 1;526(7571):68–74. 
3.	McBride N, Clayton GL, Goncalves Soares A, Yang Q, Bond TA, Taylor A, et al. Cohort Profile: The Mendelian Randomization in Pregnancy (MR-PREG) collaboration - Improving evidence for prevention and treatment of adverse pregnancy and perinatal outcomes [Internet]. Epidemiology; 2025 [cited 2025 Aug 9]. Available from: http://medrxiv.org/lookup/doi/10.1101/2025.03.22.25324447 
4.	Kurki MI, Karjalainen J, Palta P, Sipilä TP, Kristiansson K, Donner KM, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature. 2023 Jan 19;613(7944):508–18. 
5.	Westergaard D, Steinthorsdottir V, Stefansdottir L, Rohde PD, Wu X, Geller F, et al. Genome-wide association meta-analysis identifies five loci associated with postpartum hemorrhage. Nat Genet. 2024 Aug;56(8):1597–603. 
6.	Burgess S, Davey Smith G, Davies NM, Dudbridge F, Gill D, Glymour MM, et al. Guidelines for performing Mendelian randomization investigations: update for summer 2023. Wellcome Open Res. 2023 Aug 4;4:186. 
