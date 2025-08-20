# EndoMR-PREG: Mendelian Randomization Analysis

## Association of Genetic Liability to Endometriosis with Pregnancy Outcomes

[![GitHub](https://img.shields.io/github/license/jonasvibert/endoMR-PREG)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](.)

**Authors**: Jonas Vibert¹, Carolina Borges², Zoltán Kutalik³, David Baud¹, Deborah A. Lawlor², Nicola Pluchino¹

¹Department of Obstetrics and Gynecology, Lausanne University Hospital, Switzerland  
²MRC Integrative Epidemiology Unit, University of Bristol, UK  
³University Center for Primary Care and Public Health, Lausanne, Switzerland

A comprehensive two-sample Mendelian randomization study investigating causal relationships between genetic liability to endometriosis and **49 maternal and fetal outcomes**.

## 🎯 Key Findings

**Significant causal association (Bonferroni-corrected)**:
- **Placenta praevia**: OR 1.62 (95% CI: 1.33–1.97, P = 1.5×10⁻⁶)

**Internal validation**:
- **Female infertility**: OR 1.63 (95% CI: 1.48–1.80, P = 6.1×10⁻²²)

**No significant associations** with hypertensive disorders, gestational diabetes, postpartum hemorrhage, preterm birth, or most other pregnancy outcomes after correction.

## 📁 Repository Structure

```
endoMR-PREG/
├── 📖 docs/                     # Complete documentation
│   ├── README.md               # Documentation overview
│   ├── methodology.md          # Detailed methods
│   ├── results.md             # Main findings
│   ├── code-documentation.md  # Script guide
│   ├── data-sources.md        # GWAS data info
│   └── supplementary.md       # Additional analyses
├── 🔬 scripts/                 # Analysis pipeline (R)
│   ├── 1_selectSNPs_Rahmiglu_endoMR-PREG.R
│   ├── 2_prep_data_endoMR-PREG_110725.R
│   ├── 3_Harmonise data_endoMR-PREG_110725.R
│   ├── 4_MR_analysis_endoMR-PREG_150725.R
│   ├── 4.2_MR_fetaleffect_endoMR-PREG.R
│   ├── 5_Tables_MR_analysis_endoMR-PREG.R
│   ├── 6_Plots_MR_analysis_endoMR-PREG.R
│   ├── 7_visualization_selectedSNPs.R
│   └── mr_visualization_analysis.R
├── 📊 data/                    # Input GWAS data
├── 📈 results/                 # Analysis outputs
├── 🎨 plot/                    # Publication figures
└── 📋 Supplementary files
```

## 🚀 Quick Start

### Prerequisites

```r
# Required R packages
install.packages(c("TwoSampleMR", "dplyr", "ggplot2", "data.table"))
```

### Running the Analysis

1. **Clone the repository**:
   ```bash
   git clone https://github.com/jonasvibert/endoMR-PREG.git
   cd endoMR-PREG
   ```

2. **Run analysis pipeline** (in order):
   ```r
   # 1. Prepare instruments
   source("scripts/1_selectSNPs_Rahmiglu_endoMR-PREG.R")
   
   # 2. Prepare outcomes  
   source("scripts/2_prep_data_endoMR-PREG_110725.R")
   
   # 3. Harmonize data
   source("scripts/3_Harmonise data_endoMR-PREG_110725.R")
   
   # 4. Main MR analysis
   source("scripts/4_MR_analysis_endoMR-PREG_150725.R")
   
   # 5. Generate outputs
   source("scripts/5_Tables_MR_analysis_endoMR-PREG.R")
   source("scripts/6_Plots_MR_analysis_endoMR-PREG.R")
   ```

3. **Check results**:
   - Tables: `results/`
   - Figures: `plot/`

## 📖 Documentation

For detailed information, see the [`docs/`](docs/) directory:

- **[Methodology](docs/methodology.md)**: Statistical methods and study design
- **[Results](docs/results.md)**: Main findings and interpretations
- **[Code Guide](docs/code-documentation.md)**: Script documentation
- **[Data Sources](docs/data-sources.md)**: GWAS datasets used
- **[Supplementary](docs/supplementary.md)**: Additional analyses

## 🔬 Methods Overview

- **Design**: Two-sample Mendelian randomization
- **Exposure**: Endometriosis genetic liability (41 independent SNPs, mean F = 279)
- **Outcomes**: 49 maternal and fetal pregnancy outcomes
- **Population**: European ancestry (~58,000 cases, 733,000 controls)
- **Software**: TwoSampleMR package in R

**Statistical approaches**:
- Inverse Variance Weighted (primary)
- MR-Egger regression
- Weighted median
- Sensitivity analyses (heterogeneity, pleiotropy, leave-one-out)

## 📊 Key Results

| Outcome | OR (95% CI) | P-value | Method |
|---------|-------------|---------|---------|
| **Placenta praevia** | **1.62 (1.33–1.97)** | **1.5×10⁻⁶** | IVW |
| **Female infertility** | **1.63 (1.48–1.80)** | **6.1×10⁻²²** | IVW |
| Premature rupture of membranes | 1.12 (1.01–1.23) | 0.026 | IVW |
| Placental abruption | 1.36 (1.03–1.81) | 0.031 | IVW |
| Preterm birth (all) | 0.82 (0.71–0.95) | 0.007 | IVW |

*Full results available in [`results/`](results/) directory*

## 🎨 Visualizations

High-resolution publication-ready figures available in [`plot/`](plot/):

- **Forest plots**: Individual and summary effect estimates
- **Scatter plots**: Exposure-outcome relationships with multiple methods
- **Funnel plots**: Assessment of directional pleiotropy
- **Leave-one-out**: Sensitivity to individual SNPs

## 💾 Data Availability

- **Analysis code**: Fully open source (this repository)
- **Summary statistics**: Available from original GWAS sources
- **Processed datasets**: Available upon request
- **Individual-level data**: Not available (summary statistics only)

## 📜 Citation

If you use this work, please cite:

```
[Your Citation]
EndoMR-PREG: Mendelian Randomization Analysis of Endometriosis and Pregnancy Outcomes
GitHub: https://github.com/jonasvibert/endoMR-PREG
```

**Original data sources to cite**:
- Endometriosis GWAS: Rahmioglu et al. (2023) Nature Genetics
- MR-PREG consortium: McBride et al. (2025)
- FinnGen: Kurki et al. (2023) Nature
- Postpartum hemorrhage: Westergaard et al. (2024) Nature Genetics
- TwoSampleMR: Hemani et al. (2018) eLife

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