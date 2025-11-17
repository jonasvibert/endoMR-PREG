# EndoMR-PREG: Two-Sample Mendelian Randomization Study
## Association of Genetic Liability to Endometriosis With Pregnancy Outcomes

[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue)](https://github.com/jonasvibert/endoMR-PREG) [![R](https://img.shields.io/badge/R-v4.3.2+-blue)](https://www.r-project.org/) [![DOI](https://img.shields.io/badge/DOI-Pending-yellow)]()

## 👥 Authors

<table>
<tr>
<td align="left"><b>Author</b></td>
<td align="left"><b>Affiliation</b></td>
<td align="left"><b>Role</b></td>
</tr>
<tr>
<td><strong>Jonas Vibert, MD</strong></td>
<td>🏥 Department of Obstetrics and Gynecology<br>Lausanne University Hospital, Switzerland</td>
<td>📧 <em>Corresponding Author</em></td>
</tr>
<tr>
<td><strong>Carolina Borges, PhD</strong></td>
<td>🔬 MRC Integrative Epidemiology Unit<br>University of Bristol, UK</td>
<td>Co-investigator</td>
</tr>
<tr>
<td><strong>Zoltán Kutalik, PhD</strong></td>
<td>📊 University Center for Primary Care and Public Health<br>Lausanne, Switzerland</td>
<td>Senior Investigator</td>
</tr>
<tr>
<td><strong>David Baud, MD, PhD</strong></td>
<td>🏥 Department of Obstetrics and Gynecology<br>Lausanne University Hospital, Switzerland</td>
<td>Senior Author</td>
</tr>
<tr>
<td><strong>Deborah A. Lawlor, PhD</strong></td>
<td>🔬 MRC Integrative Epidemiology Unit<br>University of Bristol, UK</td>
<td>Senior Investigator</td>
</tr>
<tr>
<td><strong>Nicola Pluchino, MD, PhD</strong></td>
<td>🏥 Department of Obstetrics and Gynecology<br>Lausanne University Hospital, Switzerland</td>
<td>Principal Investigator</td>
</tr>
</table>

### 🏛️ Institutional Affiliations

**🇨🇭 Switzerland**
- 🏥 **Lausanne University Hospital (CHUV)** - Department of Obstetrics and Gynecology
- 📊 **University Center for Primary Care and Public Health (Unisanté)** - Lausanne

**🇬🇧 United Kingdom** 
- 🔬 **University of Bristol** - MRC Integrative Epidemiology Unit

---

This repository contains the complete analysis code and results for a comprehensive two-sample Mendelian randomization study investigating causal relationships between genetic liability to endometriosis and **29 specialized pregnancy and perinatal outcomes**.

## 🎯 Key Findings

### Primary Result (FDR-corrected significant, q < 0.05):
- **Placenta praevia:** OR 1.62 (95% CI: 1.33–1.97, P = 1.5×10⁻⁶, q < 0.001)

### Nominally Significant Associations (P < 0.05, but q > 0.05 after FDR correction):
- **Premature placental separation:** OR 1.36 (95% CI: 1.03–1.81, P = 0.031)
- **Elective cesarean delivery:** OR 1.26 (95% CI: 1.03–1.54, P = 0.024)
- **Premature rupture of membranes:** OR 1.12 (95% CI: 1.01–1.23, P = 0.025)
- **Preterm birth:** OR 0.83 (95% CI: 0.71–0.97, P = 0.019) *[Protective effect]*

### No Significant Associations Found With:
Preeclampsia, gestational diabetes, postpartum hemorrhage, gestational hypertension, fetal growth outcomes (SGA/LGA), NICU admission, low Apgar scores, or most other adverse pregnancy outcomes after FDR correction for multiple testing.

---

## 📁 Repository Structure

```
endoMR-PREG/
├── 📖 docs/                          # Documentation (planned)
├── 🔧 config/                        # Configuration files
│   ├── config.R                      # Main configuration & parameters  
│   └── utils.R                       # Utility functions
├── 🔬 scripts/                       # Complete analysis pipeline
│   ├── 01_select_instruments_endoMR-PREG.R      # SNP selection & clumping
│   ├── 02_prepare_outcomes_endoMR-PREG.R        # Outcome data preparation  
│   ├── 03_harmonise_data_endoMR-PREG.R          # Data harmonization
│   ├── 04_main_analyses_endoMR-PREG.R           # Main MR analyses
│   ├── 04.2_fetal_effect_endoMR-PREG.R          # Fetal genetic effect analysis†
│   ├── 05_sensitivity_analyses_endoMR-PREG.R    # Sensitivity analyses  
│   ├── 06_tables_endoMR-PREG.R                  # Results tables
│   ├── 07_plots_endoMR-PREG.R                   # Publication figures
│   └── 08_diagnostic_plots_endoMR-PREG.R        # Diagnostic plots
├── 📈 results/                       # Analysis outputs
│   ├── tables/                       # All result tables (CSV & HTML)
│   ├── plots/                        # Generated figures
│   └── plots instruments selection/   # Instrument characterization plots
└── 📋 LICENSE                        # MIT License
```

†*Requires additional trio-based GWAS data (implementation ready)*

---

## 🚀 Quick Start

### Prerequisites
R version 4.3.2+ with required packages:

```r
# Install required packages
pkgs <- c("TwoSampleMR", "dplyr", "ggplot2", "here", 
          "readr", "data.table", "gridExtra", "gt", "openxlsx")
install.packages(pkgs)

# Install TwoSampleMR if needed
devtools::install_github("MRCIEU/TwoSampleMR")
```

### Running the Analysis

```bash
# Clone the repository
git clone https://github.com/jonasvibert/endoMR-PREG.git
cd endoMR-PREG
```

Run the complete analysis pipeline in R:

```r
# 1. Select genetic instruments (41 independent endometriosis SNPs)
source("scripts/01_select_instruments_endoMR-PREG.R")

# 2. Prepare outcome data (29 pregnancy outcomes from 3 sources)
source("scripts/02_prepare_outcomes_endoMR-PREG.R")

# 3. Harmonize exposure and outcome data
source("scripts/03_harmonise_data_endoMR-PREG.R")

# 4. Main MR analysis (IVW, MR-Egger, Weighted Median)
source("scripts/04_main_analyses_endoMR-PREG.R")

# 5. Comprehensive sensitivity analyses
source("scripts/05_sensitivity_analyses_endoMR-PREG.R")

# 6. Generate publication-ready tables
source("scripts/06_tables_endoMR-PREG.R")

# 7. Create publication-quality figures
source("scripts/07_plots_endoMR-PREG.R")

# 8. Generate diagnostic plots for significant outcomes
source("scripts/08_diagnostic_plots_endoMR-PREG.R")
```

**View results:**
- Tables: `results/tables/*.csv`
- Main figures: `results/plots/*.png` 
- Diagnostic plots: `results/plots/diagnostics/*.png`

---

## 🔬 Methods Overview

### Study Design
**Two-sample Mendelian randomization** using the largest available GWAS datasets

### Exposure GWAS: Endometriosis
- **Source:** Rahmioglu et al., Nature Genetics 2023
- **Sample:** 60,674 cases and 701,926 controls  
- **Ancestry:** Predominantly European (~98%) + Japanese (~2%)
- **Genetic instruments:** 41 independent SNPs (P < 5×10⁻⁸, clumped at r² < 0.001)
- **Mean F-statistic:** 279 (strong instruments)
- **Explained variance:** ~5.6% of endometriosis liability

### Outcome GWAS Sources
1. **MR-PREG consortium** (McBride et al., 2025)
   - Up to 678,001 women (outcome-specific)  
   - Adverse pregnancy and perinatal outcomes
   - Sources: ALSPAC, BiB, MoBa, UKB, FinnGen

2. **FinnGen R12** (2024)  
   - Finnish population-based cohort
   - Pregnancy and fertility ICD-10 phenotypes
   - Placental complications (O43-O45)

3. **Westergaard et al.** (2024)
   - Nordic registry-based GWAS
   - Up to 331,792 women
   - Bleeding complications and postpartum hemorrhage

### Statistical Methods
- **Primary:** Inverse Variance Weighted (IVW)
- **Sensitivity:** MR-Egger, Weighted Median, Weighted Mode
- **Quality control:** Heterogeneity (Cochran's Q), pleiotropy (MR-Egger intercept), leave-one-out SNP
- **Multiple testing:** FDR correction across 29 outcomes (Benjamini-Hochberg, q < 0.05)
- **Software:** TwoSampleMR v0.5.6 in R v4.3.2

---

## 📊 Complete Results Summary

| Outcome | Data Source | OR (95% CI) | P-value | q (FDR) | Significance |
|---------|-------------|-------------|---------|---------|--------------|
| **Placenta praevia** | FinnGen R12 | **1.62 (1.33–1.97)** | **1.5×10⁻⁶** | **< 0.001** | **✓ FDR** |
| Premature placental separation | FinnGen R12 | 1.36 (1.03–1.81) | 0.031 | 0.180 | Nominal |
| Elective cesarean section | MR-PREG | 1.26 (1.03–1.54) | 0.024 | 0.180 | Nominal |
| Premature rupture of membranes | MR-PREG | 1.12 (1.01–1.23) | 0.025 | 0.180 | Nominal |
| Preterm birth <37 weeks | MR-PREG | 0.83 (0.71–0.97) | 0.019 | 0.180 | Nominal |
| Low Apgar score at 1 min | MR-PREG | 1.18 (0.99–1.41) | 0.073 | 0.265 | — |
| Emergency cesarean section | MR-PREG | 1.16 (0.99–1.37) | 0.069 | 0.265 | — |
| Very preterm birth <34 weeks | MR-PREG | 0.72 (0.51–1.00) | 0.052 | 0.232 | — |

*Complete results available in `results/tables/Table3_IVW_FDR_main_results.csv`*

---

## 🎨 Visualizations

### Available Figures:
1. **Instrument characterization** (`results/plots instruments selection/`)
   - Chromosomal distribution of endometriosis SNPs
   - Effect allele frequency distribution  
   - Forest plot of instrument effects
   - F-statistics vs R² scatter plot

2. **Main analysis figures** (`results/plots/`)
   - Forest plot: IVW estimates for all 29 outcomes
   - Multi-method forest plots (significant outcomes)
   - Enhanced visualization with log₂ transformation

3. **Diagnostic plots** (`results/plots/diagnostics/`)
   - Scatter plots (exposure vs outcome effects)
   - Single-SNP forest plots
   - Leave-one-out sensitivity plots  
   - Funnel plots (directional pleiotropy assessment)

*All plots generated automatically with publication-ready quality (300 DPI)*

---

## 🧬 Clinical Interpretation

### Primary Finding: Placenta Praevia
Genetic liability to endometriosis **increases the risk of placenta praevia by 62%** (OR 1.62), providing strong evidence for a causal relationship. This finding aligns with biological mechanisms involving:
- Impaired decidualization in women with endometriosis
- Abnormal placentation due to inflammatory environment
- Altered angiogenesis and tissue remodeling

### Secondary Observations
- **Protective effect on preterm birth:** Unexpected finding requiring further investigation
- **Increased caesarean delivery risk:** Consistent with known obstetric management patterns
- **No evidence for major complications:** Preeclampsia, gestational diabetes, and postpartum hemorrhage show no causal associations after multiple testing correction

### Clinical Implications
- Enhanced prenatal screening for placenta praevia in women with endometriosis history
- Potential for early intervention and specialized obstetric care planning
- Evidence-based counseling for pregnancy risks in endometriosis patients

---

## 💾 Data Availability

- **Analysis code:** Fully open source (this repository)
- **Summary statistics:** Available from original GWAS publications
- **Processed datasets:** Available upon reasonable request
- **Individual-level data:** Not available (summary statistics only)

### Data Sources:
- Endometriosis GWAS: [Rahmioglu et al. 2023](https://doi.org/10.1038/s41588-023-01323-z)
- MR-PREG outcomes: [McBride et al. 2025](https://doi.org/TBD)
- FinnGen R12: [Kurki et al. 2023](https://doi.org/10.1038/s41586-022-05473-8)
- PPH GWAS: [Westergaard et al. 2024](https://doi.org/10.1038/s41588-024-01801-1)

---

## 📜 Citation

If you use this work, please cite:

```
Vibert J, Borges C, Kutalik Z, Baud D, Lawlor DA, Pluchino N. 
Association of Genetic Liability to Endometriosis With Pregnancy Outcomes: 
A Two-Sample Mendelian Randomization. [Journal] [Year].

GitHub repository: https://github.com/jonasvibert/endoMR-PREG
```

---

## 🤝 Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch  
3. Submit a pull request

For major changes, please open an issue first to discuss proposed modifications.

---

## 📧 Contact

<table>
<tr>
<td rowspan="4" width="100px" align="center">
<strong>📧<br>Corresponding<br>Author</strong>
</td>
<td><strong>Dr. Jonas Vibert, MD</strong></td>
</tr>
<tr>
<td>✉️ <strong>Email:</strong> <a href="mailto:jonas.vibert@chuv.ch">jonas.vibert@chuv.ch</a></td>
</tr>
<tr>
<td>🆔 <strong>ORCID:</strong> <a href="https://orcid.org/0009-0000-2449-7734" target="_blank">0009-0000-2449-7734</a></td>
</tr>
<tr>
<td>🏥 <strong>Institution:</strong> Lausanne University Hospital (CHUV)<br>&nbsp;&nbsp;&nbsp;&nbsp;Department of Obstetrics and Gynecology</td>
</tr>
</table>

### 🤝 Research Collaboration Inquiries
For questions about methodology, data access, or collaboration opportunities, please contact the corresponding author.

---

## ⚖️ License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- MR-PREG consortium participants and investigators
- FinnGen consortium participants and investigators  
- Westergaard et al. for postpartum hemorrhage GWAS
- Rahmioglu et al. for endometriosis GWAS
- All study participants and GWAS consortiums
- TwoSampleMR development team (Hemani et al.)

---

*Last updated: November 2025*