# EndoMR-PREG: Mendelian Randomization Analysis

## Endometriosis and Pregnancy Outcomes Study

[![GitHub](https://img.shields.io/github/license/jonasvibert/endoMR-PREG)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](.)

A comprehensive two-sample Mendelian randomization analysis investigating the causal relationship between genetic liability to endometriosis and various pregnancy outcomes.

## 🎯 Key Findings

**Significant causal associations identified**:
- **Female infertility** (FinnGen)
- **Premature rupture of membranes**
- **Placenta praevia** (FinnGen)
- **Placental abruption** (FinnGen)

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
- **Exposure**: Endometriosis genetic liability (41 SNPs)
- **Outcomes**: 25+ pregnancy-related traits
- **Population**: European ancestry
- **Software**: TwoSampleMR package in R

**Statistical approaches**:
- Inverse Variance Weighted (primary)
- MR-Egger regression
- Weighted median
- Sensitivity analyses (heterogeneity, pleiotropy, leave-one-out)

## 📊 Key Results

| Outcome | OR (95% CI) | P-value | Method |
|---------|-------------|---------|---------|
| Female infertility | [Add] | [Add] | IVW |
| PROM | [Add] | [Add] | IVW |
| Placenta praevia | [Add] | [Add] | IVW |
| Placental abruption | [Add] | [Add] | IVW |

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
- Endometriosis GWAS: Rahmioglu et al. (2018)
- FinnGen: Kurki et al. (2023)
- TwoSampleMR: Hemani et al. (2018)

## 🤝 Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## 📧 Contact

**Author**: [Your Name]  
**Email**: [Your Email]  
**Institution**: [Your Institution]

## ⚖️ License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- FinnGen consortium and participants
- UK Biobank and participants
- Original GWAS authors and consortiums
- TwoSampleMR development team

---

**Last updated**: August 2025