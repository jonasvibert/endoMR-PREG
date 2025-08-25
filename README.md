# EndoMR-PREG: Two-Sample Mendelian Randomization Study

## Association of Genetic Liability to Endometriosis With Pregnancy Outcomes

[![GitHub](https://img.shields.io/github/license/jonasvibert/endoMR-PREG)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.3.2-blue)](https://www.r-project.org/)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](.)

**Authors**: Jonas Vibert, MD¹, Carolina Borges, PhD², Zoltán Kutalik, PhD³, David Baud, MD, PhD¹, Deborah A. Lawlor, PhD², Nicola Pluchino, MD, PhD¹

¹Department of Obstetrics and Gynecology, Lausanne University Hospital, Lausanne, Switzerland  
²MRC Integrative Epidemiology Unit at the University of Bristol, Bristol, UK  
³University Center for Primary Care and Public Health, Lausanne, Switzerland

This repository contains the complete analysis code and results for a two-sample Mendelian randomization study investigating causal relationships between genetic liability to endometriosis and **49 maternal and fetal pregnancy outcomes**.

## 🎯 Key Findings

**Bonferroni-corrected significant associations (P < 0.001)**:
- **Placenta praevia**: OR 1.62 (95% CI: 1.33–1.97, P = 1.5×10⁻⁶)
- **Female infertility**: OR 1.63 (95% CI: 1.48–1.80, P = 6.1×10⁻²²) *(internal validation)*

**Nominally significant associations (P < 0.05)**:
- Premature rupture of membranes: OR 1.12 (95% CI: 1.01–1.23, P = 0.025)
- Elective cesarean delivery: OR 1.26 (95% CI: 1.04–1.53, P = 0.023)  
- Placental abruption: OR 1.36 (95% CI: 1.03–1.81, P = 0.031)
- Reduced 1-minute Apgar score: OR 0.95 (95% CI: 0.91–0.99, P = 0.029)

**No significant associations** found with preeclampsia, gestational diabetes, postpartum hemorrhage, stillbirth, fetal growth restriction, or most other adverse pregnancy outcomes after correction for multiple testing.

## 📁 Repository Structure

```
endoMR-PREG/
├── 📖 docs/                     # Complete documentation
│   ├── methodology.md          # Two-sample MR methods
│   ├── results.md             # Main findings and interpretation
│   ├── code-documentation.md  # Analysis pipeline guide
│   ├── data-sources.md        # GWAS datasets information
│   └── supplementary.md       # Sensitivity analyses
├── 🔬 scripts/                 # Complete analysis pipeline
│   ├── 1_selectSNPs_Rahmiglu_endoMR-PREG.R      # SNP selection and clumping
│   ├── 2_prep_data_endoMR-PREG_110725.R         # Outcome data preparation  
│   ├── 4_MR_analysis_endoMR-PREG_150725.R       # Main MR analyses
│   ├── 5_Tables_MR_analysis_endoMR-PREG.R       # Results tables
│   ├── 6_PlotsVF_MR_endoMR-PREG.R              # Publication figures
│   └── 1.2_visualization_selectedSNPs.R        # SNP visualization
├── 📈 results/                 # Analysis outputs
│   ├── harmonised_rahmioglu_bpo.csv            # Harmonized data
│   ├── all_mr_methods.csv                      # Complete MR results
│   ├── ivw_results.csv                         # IVW estimates
│   ├── egger_results.csv                       # MR-Egger results
│   ├── heterogeneity_results.csv               # Heterogeneity tests
│   └── pleiotropy_results.csv                  # Pleiotropy tests
├── 🎨 plots/                   # Generated figures
└── 📋 LICENSE                  # MIT License
```

## 🚀 Quick Start

### Prerequisites

**R version 4.3.2 or higher** with the following packages:

```r
# Install required packages
pkgs <- c("TwoSampleMR", "MRPRESSO", "dplyr", "ggplot2", "here", 
          "readr", "data.table", "gridExtra", "cowplot")
install.packages(pkgs)

# For MR-PRESSO (if not available on CRAN)
devtools::install_github("rondolab/MR-PRESSO")
```

### Running the Analysis

1. **Clone the repository**:
   ```bash
   git clone https://github.com/jonasvibert/endoMR-PREG.git
   cd endoMR-PREG
   ```

2. **Run the complete analysis pipeline**:
   ```r
   # 1. Select genetic instruments (41 independent SNPs)
   source("scripts/1_selectSNPs_Rahmiglu_endoMR-PREG.R")
   
   # 2. Prepare outcome data (49 pregnancy outcomes)
   source("scripts/2_prep_data_endoMR-PREG_110725.R")
   
   # 3. Main MR analysis (IVW, Egger, Weighted Median, MR-PRESSO)
   source("scripts/4_MR_analysis_endoMR-PREG_150725.R")
   
   # 4. Generate results tables
   source("scripts/5_Tables_MR_analysis_endoMR-PREG.R")
   
   # 5. Create publication figures
   source("scripts/6_PlotsVF_MR_endoMR-PREG.R")
   ```

3. **View results**:
   - **Tables**: `results/*.csv`
   - **Figures**: `plots/*.png`

## 📖 Documentation

For detailed information, see the [`docs/`](docs/) directory:

- **[Methodology](docs/methodology.md)**: Statistical methods and study design
- **[Results](docs/results.md)**: Main findings and interpretations
- **[Code Guide](docs/code-documentation.md)**: Script documentation
- **[Data Sources](docs/data-sources.md)**: GWAS datasets used
- **[Supplementary](docs/supplementary.md)**: Additional analyses

## 🔬 Methods Overview

- **Study design**: Two-sample Mendelian randomization
- **Exposure GWAS**: Endometriosis (Rahmioglu et al., Nature Genetics 2023)
  - ~58,000 cases and 733,000 controls of European ancestry
  - 41 independent SNPs as genetic instruments (mean F-statistic = 279)
  - Explained variance: ~5.6% of endometriosis liability
- **Outcome sources**: 
  - MR-PREG consortium (up to 678,001 women)
  - FinnGen R12 (176,899 participants, predominantly Finnish)
  - Westergaard et al. for postpartum hemorrhage (~5,000 cases, ~170,000 controls)
- **Statistical methods**:
  - **Primary**: Inverse Variance Weighted (IVW)
  - **Sensitivity**: MR-Egger, Weighted Median, MR-PRESSO
  - **Quality control**: Heterogeneity tests, pleiotropy assessment, leave-one-out
  - **Multiple testing**: Bonferroni correction for 49 outcomes (P < 0.001)
- **Software**: TwoSampleMR v0.5.6, MR-PRESSO v1.0.0 in R v4.3.2

## 📊 Key Results

| Outcome | OR (95% CI) | P-value | Significance |
|---------|-------------|---------|--------------|
| **Placenta praevia** | **1.62 (1.33–1.97)** | **1.5×10⁻⁶** | **✓ Bonferroni** |
| **Female infertility** | **1.63 (1.48–1.80)** | **6.1×10⁻²²** | **✓ Bonferroni** |
| Premature rupture of membranes | 1.12 (1.01–1.23) | 0.025 | Nominal |
| Elective cesarean delivery | 1.26 (1.04–1.53) | 0.023 | Nominal |
| Placental abruption | 1.36 (1.03–1.81) | 0.031 | Nominal |
| 1-minute Apgar score | 0.95 (0.91–0.99) | 0.029 | Nominal |
| Preterm birth (all) | 0.82 (0.71–0.95) | 0.007 | Nominal |

### Clinical Interpretation

**Primary finding**: Genetic liability to endometriosis increases the risk of **placenta praevia** by 62%, supporting a causal relationship. This finding aligns with biological mechanisms involving impaired decidualization and abnormal placentation in women with endometriosis.

**No evidence** for causal effects on preeclampsia, gestational diabetes, postpartum hemorrhage, stillbirth, or fetal growth restriction, suggesting that previously reported associations may reflect bias rather than causality.

*Complete results available in [`results/all_mr_methods.csv`](results/all_mr_methods.csv)*

## 🎨 Visualizations

High-resolution publication-ready figures available in [`plots/`](plots/):

- **Summary forest plot**: `summary_forest_all_IVW.png` - Overview of all 49 outcomes
- **Scatter plots**: SNP effects for significant outcomes (placenta praevia, infertility)
- **Forest plots**: Individual SNP and combined estimates 
- **Leave-one-out plots**: Sensitivity analyses removing each SNP
- **Funnel plots**: Assessment of directional pleiotropy

All plots generated with `6_PlotsVF_MR_endoMR-PREG.R` script.

## 💾 Data Availability

- **Analysis code**: Fully open source (this repository)
- **Summary statistics**: Available from original GWAS sources
- **Processed datasets**: Available upon request
- **Individual-level data**: Not available (summary statistics only)

## 📜 Citation

If you use this work, please cite:

```
Vibert J, Borges C, Kutalik Z, Baud D, Lawlor DA, Pluchino N. 
Association of Genetic Liability to Endometriosis With Pregnancy Outcomes: 
A Two-Sample Mendelian Randomization. [Journal] [Year].

GitHub repository: https://github.com/jonasvibert/endoMR-PREG
```

**Original data sources**:
- Rahmioglu N, et al. The genetic basis of endometriosis and comorbidity with other pain and inflammatory conditions. *Nat Genet*. 2023;55(3):423-36.
- McBride N, et al. Cohort Profile: The Mendelian Randomization in Pregnancy (MR-PREG) collaboration. 2025.
- Kurki MI, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. *Nature*. 2023;613(7944):508-18.
- Westergaard D, et al. Genome-wide association meta-analysis identifies five loci associated with postpartum hemorrhage. *Nat Genet*. 2024;56(8):1597-603.
- Hemani G, et al. The MR-Base platform supports systematic causal inference across the human phenome. *eLife*. 2018;7:e34408.

## 🤝 Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## 📧 Contact

**Corresponding Author**: Jonas Vibert, MD  
**Email**: jonas.vibert@chuv.ch  
**ORCID**: https://orcid.org/0009-0000-2449-7734  
**Institution**: Lausanne University Hospital

## ⚖️ License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- MR-PREG consortium and participants
- FinnGen consortium and participants
- Westergaard et al. for postpartum hemorrhage GWAS
- Rahmioglu et al. for endometriosis GWAS
- All study participants and original GWAS consortiums
- TwoSampleMR development team

---

**Last updated**: August 2025