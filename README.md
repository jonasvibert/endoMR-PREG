# EndoMR-PREG: Two-Sample Mendelian Randomization Study

## Association of Genetic Liability to Endometriosis With Pregnancy Outcomes

[![GitHub](https://img.shields.io/github/license/jonasvibert/endoMR-PREG)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.3.2-blue)](https://www.r-project.org/)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](.)

**Authors**: Jonas Vibert, MD¹, Carolina Borges, PhD², Zoltán Kutalik, PhD³, David Baud, MD, PhD¹, Deborah A. Lawlor, PhD², Nicola Pluchino, MD, PhD¹

¹Department of Obstetrics and Gynecology, Lausanne University Hospital, Lausanne, Switzerland  
²MRC Integrative Epidemiology Unit at the University of Bristol, Bristol, UK  
³University Center for Primary Care and Public Health, Lausanne, Switzerland

This repository contains the complete analysis code and results for a two-sample Mendelian randomization study investigating causal relationships between genetic liability to endometriosis and **29 specialized pregnancy and maternal outcomes**, including novel **fetal genetic effect correction** analysis.

## 🎯 Key Findings

**FDR-corrected significant associations (q < 0.05)**:
- **Placenta praevia**: OR 1.62 (95% CI: 1.33–1.97, P = 1.5×10⁻⁶, **q = 3.6×10⁻⁵**)
- **Female infertility**: OR 1.63 (95% CI: 1.48–1.80, P = 6.1×10⁻²²) *(internal validation)*

**FDR suggestive associations (0.05 < q < 0.2)**:
- Premature rupture of membranes: OR 1.12 (95% CI: 1.01–1.23, P = 0.025, q = 0.155)
- Elective cesarean delivery: OR 1.26 (95% CI: 1.04–1.53, P = 0.023, q = 0.155)  
- Premature placental separation: OR 1.36 (95% CI: 1.03–1.81, P = 0.031, q = 0.155)
- Preterm birth: OR 0.83 (95% CI: 0.71–0.97, P = 0.019, q = 0.155)

**No significant associations** found with preeclampsia, gestational diabetes, postpartum hemorrhage, stillbirth, fetal growth restriction, or most other adverse pregnancy outcomes after FDR correction for multiple testing.

**Novel Fetal Genetic Effect Analysis (TRIOS)**:
- **Birth weight**: Direct fetal genetic effect (β = -0.045, P = 0.007) stronger than maternal effect (β = -0.021, P = 0.21)
- **Maternal vs Fetal pathways**: Most pregnancy complications driven by maternal genetic effects; birth weight primarily influenced by fetal genetics
- **Validation**: Paternal genetic effects served as negative controls (minimal associations), confirming pathway specificity

## 📁 Repository Structure

```
endoMR-PREG/
├── 📖 docs/                     # Complete documentation
│   ├── methodology.md          # Two-sample MR methods
│   ├── results.md             # Main findings and interpretation
│   ├── fetal-genetic-analysis.md # Fetal genetic effect methodology
│   ├── code-documentation.md  # Analysis pipeline guide
│   ├── data-sources.md        # GWAS datasets information
│   └── supplementary.md       # Sensitivity analyses
├── 🔬 scripts/                 # Complete analysis pipeline (29 outcomes)
│   ├── 01_select_instruments_endoMR-PREG.R     # SNP selection and clumping
│   ├── 02_prepare_outcomes_endoMR-PREG.R       # Outcome data preparation  
│   ├── 03_harmonise_data_endoMR-PREG.R         # Data harmonization
│   ├── 04_main_analyses_endoMR-PREG.R          # Main MR analyses (29 outcomes)
│   ├── 04.2_fetal_effect_endoMR-PREG.R         # 🆕 Fetal genetic effect analysis
│   ├── 05_sensitivity_analyses_endoMR-PREG.R   # Sensitivity analyses
│   ├── 06_tables_endoMR-PREG.R                 # Results tables
│   ├── 07_plots_endoMR-PREG.R                  # Publication figures
│   ├── 08_forest_plots_endoMR-PREG.R           # Forest plot generation
│   └── 09_supplementary_tables.R               # Supplementary tables
├── 📈 results/                 # Analysis outputs
│   ├── harmonised_rahmioglu_bpo.csv            # Harmonized data
│   ├── all_mr_methods.csv                      # Complete MR results (29 outcomes)
│   ├── ivw_results.csv                         # IVW estimates
│   ├── trios_maternal_fetal_paternal_summary.csv # 🆕 TRIOS analysis results
│   ├── heterogeneity_results.csv               # Heterogeneity tests
│   └── pleiotropy_results.csv                  # Pleiotropy tests
├── 🧬 fetal_genetic_analysis/   # 🆕 Fetal genetic effect outputs
│   ├── mat_fetal_paternal_comparison.png        # TRIOS comparison plot
│   └── sensitivity/                            # Leave-one-out plots
├── 🎨 plots/                   # Generated figures (enhanced forest plots)
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
   source("scripts/01_select_instruments_endoMR-PREG.R")
   
   # 2. Prepare outcome data (29 pregnancy outcomes)
   source("scripts/02_prepare_outcomes_endoMR-PREG.R")
   
   # 3. Harmonize exposure and outcome data
   source("scripts/03_harmonise_data_endoMR-PREG.R")
   
   # 4. Main MR analysis (IVW, Egger, Weighted Median)
   source("scripts/04_main_analyses_endoMR-PREG.R")
   
   # 4.2. 🆕 Fetal genetic effect analysis (TRIOS methodology)
   source("scripts/04.2_fetal_effect_endoMR-PREG.R")
   
   # 5. Sensitivity analyses
   source("scripts/05_sensitivity_analyses_endoMR-PREG.R")
   
   # 6. Generate results tables
   source("scripts/06_tables_endoMR-PREG.R")
   
   # 7. Create publication figures
   source("scripts/07_plots_endoMR-PREG.R")
   
   # 8. Generate forest plots
   source("scripts/08_forest_plots_endoMR-PREG.R")
   
   # 9. Create supplementary tables
   source("scripts/09_supplementary_tables.R")
   ```

3. **View results**:
   - **Tables**: `results/*.csv`
   - **Figures**: `plots/*.png`

## 📖 Documentation

For detailed information, see the [`docs/`](docs/) directory:

- **[Methodology](docs/methodology.md)**: Statistical methods and study design
- **[Results](docs/results.md)**: Main findings and interpretations
- **[Fetal Genetic Analysis](docs/fetal-genetic-analysis.md)**: 🆕 TRIOS methodology and maternal vs fetal effects
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
  - MR-PREG consortium (up to 678,001 women) - 29 specialized pregnancy outcomes
  - FinnGen R12 (176,899 participants, predominantly Finnish)
  - Westergaard et al. for postpartum hemorrhage (~5,000 cases, ~170,000 controls)
- **Statistical methods**:
  - **Primary**: Inverse Variance Weighted (IVW)
  - **Sensitivity**: MR-Egger, Weighted Median, Mode-based methods
  - **🆕 Fetal genetic effects**: TRIOS analysis with maternal/fetal/paternal genotypes
  - **Quality control**: Heterogeneity tests, pleiotropy assessment, leave-one-out
  - **Multiple testing**: FDR correction for 29 outcomes (q < 0.05)
- **Software**: TwoSampleMR v0.5.6, MR-PRESSO v1.0.0 in R v4.3.2

## 📊 Key Results

| Outcome | OR (95% CI) | P-value | Significance |
|---------|-------------|---------|--------------|
| **Placenta praevia** | **1.62 (1.33–1.97)** | **1.5×10⁻⁶ (q=3.6×10⁻⁵)** | **✓ FDR** |
| **Female infertility** | **1.63 (1.48–1.80)** | **6.1×10⁻²²** | **✓ Internal validation** |
| Premature rupture of membranes | 1.12 (1.01–1.23) | 0.025 | Nominal |
| Elective cesarean delivery | 1.26 (1.04–1.53) | 0.023 | Nominal |
| Placental abruption | 1.36 (1.03–1.81) | 0.031 | Nominal |
| Low Apgar score at 1 minute | 1.18 (0.98–1.41) | 0.073 | Nominal |
| **🆕 Birth weight (fetal effect)** | **β = -0.045** | **0.007** | **Fetal pathway** |
| **🆕 Birth weight (maternal effect)** | **β = -0.021** | **0.21** | **Maternal pathway** |

### Clinical Interpretation

**Primary finding**: Genetic liability to endometriosis increases the risk of **placenta praevia** by 62%, supporting a causal relationship. This finding aligns with biological mechanisms involving impaired decidualization and abnormal placentation in women with endometriosis.

**🆕 Fetal genetic pathway discovery**: Birth weight represents a unique exception where **direct fetal genetic effects** (β = -0.045, P = 0.007) are stronger than maternal effects (β = -0.021, P = 0.21), indicating that endometriosis genetic liability can directly affect fetal growth independent of maternal pathways.

**Pathway specificity**: Most pregnancy complications are primarily driven by **maternal genetic pathways**, while birth weight shows predominant **fetal genetic effects**, suggesting different biological mechanisms underlying various pregnancy outcomes.

**No evidence** for causal effects on preeclampsia, gestational diabetes, postpartum hemorrhage, stillbirth after multiple testing correction, suggesting that previously reported associations may reflect bias rather than causality.

*Complete results available in [`results/all_mr_methods.csv`](results/all_mr_methods.csv)*

## 🎨 Visualizations

High-resolution publication-ready figures available in [`plots/`](plots/) and [`fetal_genetic_analysis/`](fetal_genetic_analysis/):

**Main Analysis Plots**:
- **Enhanced forest plots**: `forest_endoMR-PREG_*.png` - Multiple layouts (standard, improved, multi-method)
- **Scatter plots**: SNP effects for significant outcomes (placenta praevia, infertility)
- **Leave-one-out plots**: Sensitivity analyses removing each SNP
- **Funnel plots**: Assessment of directional pleiotropy

**🆕 Fetal Genetic Analysis Plots**:
- **TRIOS comparison**: `mat_fetal_paternal_comparison.png` - Maternal vs fetal vs paternal effects
- **Sensitivity analysis**: `sensitivity/` - Leave-one-out plots by outcome for TRIOS analysis

All main plots generated with scripts `07_plots_endoMR-PREG.R` and `08_forest_plots_endoMR-PREG.R`.  
Fetal analysis plots generated with `04.2_fetal_effect_endoMR-PREG.R`.

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

**Last updated**: November 2025