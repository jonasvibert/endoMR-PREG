# Code Documentation

## Analysis Pipeline Guide

### Overview

The analysis consists of 9 main R scripts that should be run sequentially. Each script has a specific purpose and generates outputs for subsequent steps. The pipeline analyzes **29 specialized pregnancy and maternal outcomes** with novel **fetal genetic effect correction**.

### Script Execution Order

#### 1. Data Preparation Scripts

**Script 01: `01_select_instruments_endoMR-PREG.R`**
- **Purpose**: Select and format endometriosis genetic instruments
- **Input**: Raw endometriosis GWAS data (Rahmioglu et al., 2023)
- **Output**: 41 clumped, genome-wide significant SNPs
- **Key functions**: SNP selection, LD clumping, F-statistic calculation, instrument validation

**Script 02: `02_prepare_outcomes_endoMR-PREG.R`**
- **Purpose**: Prepare 29 pregnancy outcome datasets for MR analysis
- **Input**: Raw GWAS summary statistics from MR-PREG, FinnGen, Westergaard
- **Output**: Formatted outcome data files for 29 specialized outcomes
- **Key functions**: Data formatting, quality control, outcome standardization

**Script 03: `03_harmonise_data_endoMR-PREG.R`**
- **Purpose**: Harmonize exposure and outcome data
- **Input**: Exposure instruments + outcome data
- **Output**: Harmonized dataset ready for MR analysis
- **Key functions**: Allele alignment, strand flipping, palindromic SNP handling

#### 2. Main Analysis Scripts

**Script 04: `04_main_analyses_endoMR-PREG.R`**
- **Purpose**: Main Mendelian randomization analysis (29 outcomes)
- **Input**: Harmonized data
- **Output**: MR results, sensitivity analyses for all outcomes
- **Key functions**: 
  - Primary MR methods (IVW, MR-Egger, Weighted Median, Mode-based)
  - Heterogeneity testing (Cochran's Q)
  - Pleiotropy assessment (MR-Egger intercept)
  - Leave-one-out sensitivity analysis

**🆕 Script 04.2: `04.2_fetal_effect_endoMR-PREG.R`**
- **Purpose**: Novel fetal genetic effect analysis using TRIOS methodology
- **Input**: Harmonized data + family trio genotypes (maternal/fetal/paternal)
- **Output**: Fetal vs maternal vs paternal effect comparisons
- **Key functions**: 
  - TRIOS analysis implementation
  - Maternal/fetal/paternal genetic effect decomposition
  - Pathway identification (maternal vs fetal)
  - Negative control validation (paternal effects)

**Script 05: `05_sensitivity_analyses_endoMR-PREG.R`**
- **Purpose**: Extended sensitivity analyses for all 29 outcomes
- **Input**: MR analysis results
- **Output**: Comprehensive sensitivity testing results
- **Key functions**: Additional robustness checks, method comparisons

#### 3. Output Generation Scripts

**Script 06: `06_tables_endoMR-PREG.R`**
- **Purpose**: Generate summary tables and statistical outputs
- **Input**: MR analysis results (29 outcomes)
- **Output**: Publication-ready tables (CSV, Excel format)
- **Key functions**: Result formatting, statistical summaries, effect size calculations

**Script 07: `07_plots_endoMR-PREG.R`**
- **Purpose**: Create publication-quality visualizations
- **Input**: MR results for significant outcomes
- **Output**: High-resolution plots (PNG format)
- **Key functions**: 
  - Scatter plots for significant outcomes
  - Funnel plots for pleiotropy assessment
  - Leave-one-out plots for sensitivity
  - Enhanced forest plots

**Script 08: `08_forest_plots_endoMR-PREG.R`**
- **Purpose**: Generate comprehensive forest plots with multiple layouts
- **Input**: IVW and multi-method MR results
- **Output**: Multiple forest plot styles
- **Key functions**:
  - Standard IVW forest plots
  - Multi-method comparison forests
  - Enhanced layout options
  - endoPAIN-style layouts

**Script 09: `09_supplementary_tables.R`**
- **Purpose**: Create supplementary tables for sensitivity analysis
- **Input**: IVW, MR-Egger, Weighted Median, heterogeneity, pleiotropy results
- **Output**: Formatted supplementary tables (CSV, Word, HTML)
- **Key functions**: Sensitivity analysis summaries, method comparisons

### Key Data Files

#### Input Data
- `data/endo_dat.txt` - Endometriosis exposure data (Rahmioglu et al., 2023)
- `data/OUTCOME_MR-PREG/` - MR-PREG consortium outcome data
- `data/OUTCOME_FINNGEN/` - FinnGen R12 outcome data  
- `data/OUTCOME_PPH/` - Westergaard postpartum hemorrhage data
- `Supplementary_Table_S1_Endometriosis_Instruments.xlsx` - Instrument details

#### Main Output Data
- `results/harmonised_rahmioglu_bpo.csv` - Main harmonized dataset (29 outcomes)
- `results/all_mr_methods.csv` - Complete MR results across all methods
- `results/ivw_results.csv` - Primary IVW MR results
- `results/heterogeneity_results.csv` - Heterogeneity statistics (Cochran's Q)
- `results/pleiotropy_results.csv` - Pleiotropy test results (MR-Egger intercept)

#### 🆕 Fetal Genetic Analysis Output
- `fetal_genetic_analysis/trios_maternal_fetal_paternal_summary.csv` - TRIOS analysis results
- `fetal_genetic_analysis/mat_fetal_paternal_comparison.png` - TRIOS comparison plot
- `fetal_genetic_analysis/sensitivity/` - Leave-one-out sensitivity plots for TRIOS

#### Visualization Output
- `plots/forest_endoMR-PREG_*.png` - Multiple forest plot layouts
- `plots/scatter_*.png` - Scatter plots for significant outcomes
- `plots/loo_*.png` - Leave-one-out sensitivity plots
- `plots/funnel_*.png` - Funnel plots for pleiotropy assessment

### Software Requirements

#### R Packages
```r
# Core MR analysis
library(TwoSampleMR)
library(MendelianRandomization)

# Data manipulation
library(dplyr)
library(data.table)
library(tidyr)

# Visualization
library(ggplot2)
library(forestplot)
library(gridExtra)

# Output formatting
library(openxlsx)
library(knitr)
```

#### System Requirements
- **R version**: ≥ 4.0.0
- **Memory**: ≥ 8GB RAM recommended
- **Storage**: ~2GB for full analysis with outputs

### Troubleshooting Common Issues

#### 1. Missing Data Errors
```r
# Check if harmonized data exists
if (!exists("dat")) {
  stop("Harmonized data not found. Run script 3 first.")
}
```

#### 2. Insufficient SNPs
```r
# Check SNP count for MR-Egger
if (nrow(dat_subset) < 3) {
  message("Insufficient SNPs for MR-Egger analysis")
}
```

#### 3. Memory Issues
```r
# For large datasets, use data.table
library(data.table)
dat <- fread("large_file.csv")
```

### Best Practices

1. **Run scripts sequentially** - Each script depends on previous outputs
2. **Check intermediate outputs** - Verify each step before proceeding
3. **Save workspace** - Use `save.image()` after major steps
4. **Document parameters** - Record any parameter changes
5. **Version control** - Use git to track changes

### Quality Control Checks

#### Before Analysis
- [ ] Instruments pass F-statistic threshold (F > 10)
- [ ] No excessive LD between instruments (r² < 0.001)
- [ ] Outcome data properly formatted

#### During Analysis
- [ ] Harmonization successful for all 29 outcomes
- [ ] No excessive heterogeneity (Q test P > 0.05)
- [ ] MR-Egger intercept test non-significant (P > 0.05)
- [ ] 🆕 TRIOS analysis: Paternal effects minimal (negative control validation)

#### After Analysis
- [ ] Results consistent across methods
- [ ] Effect directions biologically plausible
- [ ] Confidence intervals appropriately wide
- [ ] 🆕 Fetal genetic effects stronger than maternal for birth weight
- [ ] 🆕 Pathway specificity validated (maternal vs fetal effects)

---

*For specific function documentation, see individual script headers*