# Code Documentation

## Analysis Pipeline Guide

### Overview

The analysis consists of 7 main R scripts that should be run sequentially. Each script has a specific purpose and generates outputs for subsequent steps.

### Script Execution Order

#### 1. Data Preparation Scripts

**Script 1: `1_selectSNPs_Rahmiglu_endoMR-PREG.R`**
- **Purpose**: Select and format endometriosis genetic instruments
- **Input**: Raw endometriosis GWAS data
- **Output**: Clumped, genome-wide significant SNPs
- **Key functions**: SNP selection, clumping, instrument validation

**Script 2: `2_prep_data_endoMR-PREG_110725.R`**
- **Purpose**: Prepare outcome datasets for MR analysis
- **Input**: Raw GWAS summary statistics for pregnancy outcomes
- **Output**: Formatted outcome data files
- **Key functions**: Data formatting, quality control, standardization

**Script 3: `3_Harmonise data_endoMR-PREG_110725.R`**
- **Purpose**: Harmonize exposure and outcome data
- **Input**: Exposure instruments + outcome data
- **Output**: Harmonized dataset ready for MR
- **Key functions**: Allele alignment, strand flipping, palindromic SNP handling

#### 2. Analysis Scripts

**Script 4: `4_MR_analysis_endoMR-PREG_150725.R`**
- **Purpose**: Main Mendelian randomization analysis
- **Input**: Harmonized data
- **Output**: MR results, sensitivity analyses
- **Key functions**: 
  - Primary MR methods (IVW, MR-Egger, Weighted Median)
  - Heterogeneity testing
  - Pleiotropy assessment
  - Leave-one-out analysis

**Script 4.2: `4.2_MR_fetaleffect_endoMR-PREG.R`**
- **Purpose**: Fetal genetic effect analysis
- **Input**: Harmonized data + fetal genetic instruments
- **Output**: Fetal vs maternal effect comparisons
- **Key functions**: Fetal MR analysis, maternal-fetal comparison

#### 3. Output Generation Scripts

**Script 5: `5_Tables_MR_analysis_endoMR-PREG.R`**
- **Purpose**: Generate summary tables and statistical outputs
- **Input**: MR analysis results
- **Output**: Publication-ready tables (CSV, Excel)
- **Key functions**: Result formatting, statistical summaries

**Script 6: `6_Plots_MR_analysis_endoMR-PREG.R`**
- **Purpose**: Create publication-quality visualizations
- **Input**: MR results
- **Output**: High-resolution plots (PNG, PDF, SVG)
- **Key functions**: 
  - Forest plots
  - Scatter plots
  - Funnel plots
  - Leave-one-out plots

**Script 7: `7_visualization_selectedSNPs.R`**
- **Purpose**: Visualize genetic instrument characteristics
- **Input**: Selected SNPs and their properties
- **Output**: SNP characteristic plots
- **Key functions**: Instrument visualization, effect size distributions

#### 4. Comprehensive Visualization

**Script: `mr_visualization_analysis.R`**
- **Purpose**: Comprehensive automated visualization pipeline
- **Input**: Harmonized data and MR results
- **Output**: Complete set of publication-ready figures
- **Features**:
  - Automated plot generation
  - Consistent styling
  - Multiple export formats
  - Error handling

### Key Data Files

#### Input Data
- `data/endo_dat.txt` - Endometriosis exposure data
- `data/*_filtered_for_TSMR.csv` - Pregnancy outcome data
- `Supplementary_Table_S1_Endometriosis_Instruments.xlsx` - Instrument details

#### Output Data
- `results/harmonised_rahmioglu_bpo.csv` - Main harmonized dataset
- `results/ivw_results.csv` - Primary MR results
- `results/heterogeneity_results.csv` - Heterogeneity statistics
- `results/pleiotropy_results.csv` - Pleiotropy test results

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
- [ ] Harmonization successful for all outcomes
- [ ] No excessive heterogeneity (Q test P > 0.05)
- [ ] MR-Egger intercept test non-significant (P > 0.05)

#### After Analysis
- [ ] Results consistent across methods
- [ ] Effect directions biologically plausible
- [ ] Confidence intervals appropriately wide

---

*For specific function documentation, see individual script headers*