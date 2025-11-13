# Data Sources and References

## GWAS Data Sources

### Exposure Data: Endometriosis

**Primary Source**: Rahmioglu et al. (2023)
- **Publication**: "The genetic basis of endometriosis and comorbidity with other pain and inflammatory conditions"
- **Journal**: Nature Genetics
- **DOI**: 10.1038/s41588-023-01323-z
- **Sample size**: ~791,000 individuals (~58,000 cases, ~733,000 controls)
- **Population**: European ancestry
- **Genetic instruments**: 41 independent genome-wide significant SNPs
- **Explained variance**: ~5.6% of endometriosis liability
- **Mean F-statistic**: 279 (strong instruments)

**Supplementary Sources**:
- Additional endometriosis GWAS for validation
- Meta-analysis datasets when available

### Outcome Data: 29 Specialized Pregnancy and Maternal Outcomes

#### MR-PREG Consortium

**Primary Source**: McBride et al. (2025)
**Population**: Up to 678,001 women (European ancestry)
**Focus**: Comprehensive pregnancy outcome GWAS meta-analysis
**🆕 TRIOS data**: Family trio data (mother-father-offspring) for fetal genetic effect analysis

**Outcomes included**:
- Preterm birth (any), very preterm birth
- Birth weight (high >4000g, low <2500g), small for gestational age
- Gestational age, post-term birth
- Caesarean section (elective, emergency, any)
- Labor induction, premature rupture of membranes
- Low Apgar scores (1 min, 5 min), NICU admission
- Gestational diabetes mellitus, pregnancy anemia
- Hypertensive disorders of pregnancy, gestational hypertension, preeclampsia
- Postpartum depression, stillbirth

#### FinnGen Release 12

**Access**: https://www.finngen.fi/
**Population**: ~176,899 participants (Finnish ancestry)
**Sample sizes**: Variable by outcome (typically 1,000-10,000 cases)

**Specialized placental outcomes**:
- `finngen_R12_O15_PLAC_PRAEVIA` - Placenta praevia
- `finngen_R12_O15_PLAC_PREMAT_SEPAR` - Premature placental separation
- `finngen_R12_O15_PLAC_DISORD` - Placental disorders (other)

#### Westergaard et al. Postpartum Hemorrhage GWAS

**Source**: Westergaard et al. (2024)
**Publication**: "Genome-wide association meta-analysis identifies five loci associated with postpartum hemorrhage"
**Journal**: Nature Genetics
**Sample**: ~5,000 PPH cases, ~170,000 controls

**Bleeding outcomes**:
- Postpartum hemorrhage (overall)
- PPH due to uterine atony
- PPH due to retained placenta
- Antepartum bleeding

### Reference Data

#### Linkage Disequilibrium Reference

**1000 Genomes Project Phase 3**
- **Population**: European (EUR) superpopulation
- **Purpose**: LD clumping, harmonization
- **Access**: https://www.internationalgenome.org/

#### Genome Build

**Primary build**: GRCh37/hg19
**Alternative**: GRCh38/hg38 (with liftOver when necessary)

### Data Processing Steps

#### Quality Control Filters

**Exposure data (Endometriosis)**:
- P-value threshold: 5 × 10⁻⁸
- MAF threshold: > 0.01
- INFO score: > 0.8
- LD clumping: r² < 0.001, 10,000 kb window

**Outcome data**:
- SNP-level filters applied per study
- Population ancestry restrictions
- Case/control minimum thresholds

#### Harmonization Protocol

1. **Allele alignment**: Effect alleles matched between exposure and outcome
2. **Strand flipping**: Correction for strand differences
3. **Palindromic SNPs**: Excluded if MAF > 0.42
4. **Effect size orientation**: Standardized to same allele

#### 🆕 TRIOS Data Processing

**Family trio structure**: Mother-father-offspring genotypes for same SNP set
**Quality control**:
- Mendelian inheritance checks
- Family structure validation
- Population stratification control

**Analysis approach**:
- **Maternal genetic effects**: Using maternal genotypes as instruments
- **Fetal genetic effects**: Using offspring genotypes as instruments
- **Paternal genetic effects**: Using paternal genotypes as negative controls
- **Pathway validation**: Minimal paternal effects confirm approach validity

### Ethical Considerations

#### Data Usage Agreements
- All data used under appropriate academic licenses
- FinnGen: Academic research agreement
- UK Biobank: Application 12345 (example)
- No individual-level data shared

#### IRB Approval
- Summary statistics only - no additional IRB required
- Original studies obtained appropriate ethical approvals

### Data Availability Statement

**Analysis code**: Fully available in this repository
**Summary statistics**: Available from original sources (URLs provided)
**Processed datasets**: Available upon reasonable request
**Individual-level data**: Not available (summary statistics only)

### Acknowledgments

We thank the following groups for making data publicly available:
- FinnGen consortium and participants
- UK Biobank and participants  
- EGG consortium
- 1000 Genomes Project
- All original GWAS authors and participants

### Version Information

**Data versions used**:
- FinnGen: Release 12
- 1000 Genomes: Phase 3
- UK Biobank: [Add data freeze date]

**Analysis date**: November 2025
**Last data update**: [Add date]

### Citation Requirements

When using this analysis, please cite:
1. This repository/manuscript
2. Original GWAS publications (see references)
3. TwoSampleMR package
4. FinnGen (if using FinnGen outcomes)

### References

1. Rahmioglu, N., et al. (2023). The genetic basis of endometriosis and comorbidity with other pain and inflammatory conditions. Nature Genetics, 55(3):423-36. DOI: 10.1038/s41588-023-01323-z

2. McBride, N., et al. (2025). Cohort Profile: The Mendelian Randomization in Pregnancy (MR-PREG) collaboration. [In preparation]

3. Westergaard, D., et al. (2024). Genome-wide association meta-analysis identifies five loci associated with postpartum hemorrhage. Nature Genetics, 56(8):1597-603.

4. Kurki, M.I., et al. (2023). FinnGen provides genetic insights from a well-phenotyped isolated population. Nature, 613(7944):508-18.

5. Hemani, G., et al. (2018). The MR-Base platform supports systematic causal inference across the human phenome. eLife, 7:e34408.

6. Warrington, N.M., et al. (2019). Maternal and fetal genetic effects on birth weight and their relevance to cardio-metabolic risk factors. Nature Genetics, 51:804-814. [TRIOS methodology]

7. FinnGen consortium. (2023). FinnGen Documentation of R12 release. https://www.finngen.fi/

---

*For data access questions, contact: [Your contact information]*