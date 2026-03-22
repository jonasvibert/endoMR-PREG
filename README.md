# endoMR-PREG: Two-Sample Mendelian Randomization Study
## Genetic liability to endometriosis and pregnancy outcomes: a two-sample Mendelian randomization study with maternal–fetal effect decomposition

[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue)](https://github.com/jonasvibert/endoMR-PREG) [![R](https://img.shields.io/badge/R-v4.3.2+-blue)](https://www.r-project.org/) [![DOI](https://img.shields.io/badge/DOI-Pending-yellow)]()

**Authors:** Jonas Vibert¹\*, Tuck Seng Cheng², Maria Christine Magnus³, Lizzy Aiton²·⁴, Zoltán Kutalik⁵, David Baud¹, Deborah A. Lawlor²·⁴, Maria Carolina Borges²·⁴†, Nicola Pluchino¹†

¹Department of Obstetrics and Gynecology, Lausanne University Hospital, Lausanne, Switzerland
²MRC Integrative Epidemiology Unit, University of Bristol, Bristol, UK
³Centre for Fertility and Health, Norwegian Institute of Public Health, Oslo, Norway
⁴Population Health Sciences, Bristol Medical School, University of Bristol, Bristol, UK
⁵University Center for Primary Care and Public Health, Lausanne, Switzerland

\*Corresponding author — jonas.vibert@chuv.ch — [ORCID 0009-0000-2449-7734](https://orcid.org/0009-0000-2449-7734)
†These authors contributed equally (listed alphabetically)

---

This repository contains the complete analysis code for the study described below.

> *We applied Mendelian randomization (MR) to estimate the causal effects of genetic liability to endometriosis on a broad range of maternal and perinatal outcomes, with particular attention to placentation-related pathways.*

## Key Findings

> *Across 30 outcomes, only placenta praevia reached FDR-corrected significance, with a robust and consistent causal signal (IVW OR 1·62, 95% CI 1·33–1·97; q<0·001). Within the placental disorders domain, estimates for premature placental separation and the broader placental disorders phenotype were directionally concordant but imprecise. For premature rupture of membranes, estimates were concordant across three methods, though the association was sensitive to cohort exclusion and did not survive multiple testing correction and should be interpreted cautiously. By contrast, hypertensive disorders, gestational diabetes, postpartum haemorrhage, stillbirth, and most neonatal outcomes showed estimates consistently close to the null across all methods. Trio-based analyses suggested predominantly maternal genetic pathways for most outcomes; fetal genetic contributions were not significant after correction for multiple testing, with exploratory signals observed for birthweight-related outcomes requiring independent replication.*

---

## Repository Structure

```
endoMR-PREG/
├── docs/                          # Documentation
├── config/                        # Configuration files
│   ├── config.R
│   └── utils.R
├── scripts/                       # Analysis pipeline (scripts 01–07)
│   ├── 01_select_instruments_endoMR-PREG.R      # SNP selection & LD clumping
│   ├── 02_prepare_outcomes_endoMR-PREG.R        # Outcome data preparation
│   ├── 03_harmonise_data_endoMR-PREG.R          # Harmonisation (main)
│   ├── 03b_harmonise_with_proxies_endoMR-PREG.R # Harmonisation with proxies
│   ├── 04_main_analyses_endoMR-PREG.R           # IVW primary analysis
│   ├── 04.2_fetal_effect_endoMR-PREG.R          # Trio-based maternal/fetal/paternal analysis
│   ├── 05_sensitivity_analyses_endoMR-PREG.R    # MR-Egger, WM, WMode, PRESSO, LOO
│   ├── 06_tables_endoMR-PREG.R                  # Supplementary and main tables
│   └── 07_plots_endoMR-PREG.R                   # All manuscript figures (Fig 2–6, S1–S16)
└── data/                          # Input GWAS summary statistics (not tracked)
```

Generated outputs (`results/`) are not tracked in git and are recreated by running the pipeline.

---

## Quick Start

### Prerequisites

R 4.3.2+ with:

```r
pkgs <- c("TwoSampleMR", "MRPRESSO", "dplyr", "ggplot2", "patchwork",
          "here", "readr", "data.table", "openxlsx", "readxl")
install.packages(pkgs)
devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("rondolab/MR-PRESSO")
```

### Running the Pipeline

```r
source("scripts/01_select_instruments_endoMR-PREG.R")   # 41 endometriosis SNPs
source("scripts/02_prepare_outcomes_endoMR-PREG.R")     # 30 outcome GWAS
source("scripts/03_harmonise_data_endoMR-PREG.R")       # harmonise exposure/outcomes
source("scripts/04_main_analyses_endoMR-PREG.R")        # IVW primary analysis
source("scripts/04.2_fetal_effect_endoMR-PREG.R")       # trio-based decomposition
source("scripts/05_sensitivity_analyses_endoMR-PREG.R") # sensitivity & LOO analyses
source("scripts/06_tables_endoMR-PREG.R")               # Tables 2–3, S1–S5
source("scripts/07_plots_endoMR-PREG.R")                # Figures 2–6, S1–S16
```

---

## Methods Overview

### Exposure
- **Source:** Rahmioglu et al., *Nature Genetics* 2023
- **Sample:** 60,674 cases; 701,926 controls (predominantly European)
- **Instruments:** 41 independent SNPs (p < 5×10⁻⁸, r² < 0.001, 10,000 kb window, 1000G EUR)
- **Instrument strength:** Mean F-statistic = 279; variance explained ≈ 5.6%

### Outcomes (30 total across 9 clinical domains)
| Source | Outcomes | N |
|--------|----------|---|
| **MR-PREG** (McBride et al., 2025) | 23 maternal/perinatal outcomes | Up to 678,001 women |
| **FinnGen R12** (Kurki et al., 2023) | Placenta praevia, placental disorders, premature placental separation | Up to 223,001 women |
| **Westergaard et al.** (2024) | PPH (any), PPH atony, PPH retained placenta, antepartum bleeding | Up to 331,792 women |

### Statistical Methods
- **Primary:** Inverse Variance Weighted (IVW)
- **Sensitivity:** MR-Egger, weighted median, weighted mode, simple mode, MR-PRESSO
- **Pleiotropy:** Egger intercept, PRESSO global test, Cochran's Q
- **Leave-one-out:** SNP-level (S11–S16) and cohort-level (Supplementary Table S5)
- **Trio-based:** Maternal/fetal/paternal effect decomposition (21 outcomes; Figure 6, S10)
- **Multiple testing:** Benjamini–Hochberg FDR across 30 outcomes (q < 0.05)
- **Software:** TwoSampleMR v0.5.6, MRPRESSO v1.0.0, R v4.3.2

---

## Interpretation

> *The pattern of findings — a robust causal signal for placenta praevia alongside directionally consistent estimates across the placental disorders domain — suggests that mechanisms related to abnormal implantation and placentation may play a central role through which endometriosis liability influences pregnancy, rather than systemic maternal mechanisms. These results suggest that previously reported associations with broader obstetric outcomes may partly reflect confounding or clinical management patterns, and support targeted surveillance for abnormal placentation rather than a generalised elevation of obstetric risk.*

---

## Citation

```
Vibert J, Cheng TS, Magnus MC, Aiton L, Kutalik Z, Baud D, Lawlor DA,
Borges MC, Pluchino N. Genetic liability to endometriosis and pregnancy
outcomes: a two-sample Mendelian randomization study with maternal–fetal
effect decomposition. [Journal] [Year].

GitHub: https://github.com/jonasvibert/endoMR-PREG
```

---

## Data Sources

- Endometriosis GWAS: [Rahmioglu et al. 2023](https://doi.org/10.1038/s41588-023-01323-z)
- MR-PREG: [McBride et al. 2025](https://doi.org/10.1101/2025.03.22.25324447)
- FinnGen R12: [Kurki et al. 2023](https://doi.org/10.1038/s41586-022-05473-8)
- PPH GWAS: [Westergaard et al. 2024](https://doi.org/10.1038/s41588-024-01839-y)

---

## Contact

**Jonas Vibert, MD** — jonas.vibert@chuv.ch
Department of Obstetrics and Gynecology, Lausanne University Hospital (CHUV)
Rue du Bugnon 21, 1011 Lausanne, Switzerland

---

## Acknowledgements

We thank the MR-PREG consortium participants and investigators for providing access to summary-level GWAS data across pregnancy and perinatal outcomes.

## Funding

MCB, LA, and DAL are members of the MRC Integrative Epidemiology Unit at the University of Bristol (MC_UU_00032/5). No specific funding was received for this study.

---

## License

MIT License — see LICENSE file for details.

---

*Last updated: March 2026*
