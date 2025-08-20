# Data Sources and References

## GWAS Data Sources

### Exposure Data: Endometriosis

**Primary Source**: Rahmioglu et al. (2018)
- **Publication**: "Large-scale genome-wide association meta-analysis of endometriosis reveals 13 novel loci and a causal association with depression"
- **Journal**: Human Molecular Genetics
- **DOI**: [Add DOI]
- **Sample size**: [Add sample size]
- **Population**: European ancestry
- **Cases/Controls**: [Add numbers]

**Supplementary Sources**:
- Additional endometriosis GWAS for validation
- Meta-analysis datasets when available

### Outcome Data: Pregnancy Outcomes

#### FinnGen Release 12

**Access**: https://www.finngen.fi/
**Population**: Finnish
**Sample sizes**: Variable by outcome

**Specific outcomes**:
- `finngen_R12_N14_FEMALEINFERT` - Female infertility
- `finngen_R12_O15_PLAC_PRAEVIA` - Placenta praevia  
- `finngen_R12_O15_PLAC_PREMAT_SEPAR` - Placental abruption
- `finngen_R12_O15_PLAC_DISORD` - Placental disorders (general)
- `finngen_R12_O15_PREG_ECTOP` - Ectopic pregnancy

#### UK Biobank

**Access**: https://www.ukbiobank.ac.uk/
**Population**: UK (primarily European ancestry)

**Pregnancy outcomes**:
- Postpartum hemorrhage variants
- Preterm birth outcomes
- Caesarean section records

#### Early Growth Genetics (EGG) Consortium

**Website**: https://egg-consortium.org/
**Focus**: Birth outcomes and early growth

**Key outcomes**:
- Birthweight
- Gestational duration
- Birth length
- Small/large for gestational age

#### Other Consortiums

**GWAS Catalog**: https://www.ebi.ac.uk/gwas/
- Additional pregnancy outcome GWAS
- Validation datasets

**dbGaP**: https://www.ncbi.nlm.nih.gov/gap/
- US-based pregnancy cohorts
- Maternal and fetal genetic data

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

**Analysis date**: August 2025
**Last data update**: [Add date]

### Citation Requirements

When using this analysis, please cite:
1. This repository/manuscript
2. Original GWAS publications (see references)
3. TwoSampleMR package
4. FinnGen (if using FinnGen outcomes)

### References

1. Rahmioglu, N., et al. (2018). Large-scale genome-wide association meta-analysis of endometriosis reveals 13 novel loci and a causal association with depression. Human Molecular Genetics.

2. Hemani, G., et al. (2018). The MR-Base platform supports systematic causal inference across the human phenome. eLife.

3. FinnGen consortium. (2023). FinnGen Documentation of R12 release. https://www.finngen.fi/

4. Kurki, M.I., et al. (2023). FinnGen provides genetic insights from a well-phenotyped isolated population. Nature.

---

*For data access questions, contact: [Your contact information]*