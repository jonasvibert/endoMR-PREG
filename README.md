# endoMR-PREG: Two-Sample Mendelian Randomization Study
## Genetic liability to endometriosis and pregnancy outcomes: a two-sample Mendelian randomization study with maternal–fetal effect decomposition

[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue)](https://github.com/jonasvibert/endoMR-PREG) [![R](https://img.shields.io/badge/R-v4.3.2+-blue)](https://www.r-project.org/) [![DOI](https://img.shields.io/badge/DOI-Pending-yellow)](#)

**Authors:** Jonas Vibert¹\*, Tuck Seng Cheng²·⁴, Maria Christine Magnus³, Elisabeth Aiton²·⁴, Zoltán Kutalik⁵, David Baud¹, Deborah A. Lawlor²·⁴, Maria Carolina Borges²·⁴†, Nicola Pluchino¹†

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

> *Across 30 outcomes, only placenta praevia reached FDR-corrected significance, with consistent estimates across four of five sensitivity methods (IVW OR 1·62, 95% CI 1·33–1·97; q<0·001). Within the placental disorders domain, estimates for premature placental separation and the broader placental disorders phenotype were directionally concordant but imprecise. For premature rupture of membranes, estimates were concordant across three methods, though the association was sensitive to cohort exclusion and did not survive multiple testing correction. By contrast, hypertensive disorders, gestational diabetes, postpartum haemorrhage, stillbirth, and most neonatal outcomes showed estimates consistently close to the null across all methods. Trio-based analyses suggested predominantly maternal genetic pathways for most outcomes; fetal genetic contributions were not significant after correction for multiple testing.*

---

## Repository Structure

```
endoMR-PREG/
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
│   └── 07_plots_endoMR-PREG.R                   # All manuscript figures (Fig 2–3, S1–S15)
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
source("scripts/06_tables_endoMR-PREG.R")               # Tables 1–3, S1–S5
source("scripts/07_plots_endoMR-PREG.R")                # Figures 2–3, S1–S15
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
| **MR-PREG** (McBride et al., 2026) | 23 maternal/perinatal outcomes | Up to 678,001 women |
| **FinnGen R12** (Kurki et al., 2023) | Placenta praevia, placental disorders, premature placental separation | Up to 223,001 women |
| **Westergaard et al.** (2024) | PPH (any), PPH atony, PPH retained placenta, antepartum bleeding | Up to 331,792 women |

### Statistical Methods
- **Primary:** Inverse Variance Weighted (IVW)
- **Sensitivity:** MR-Egger, weighted median, weighted mode, simple mode, MR-PRESSO
- **Pleiotropy:** Egger intercept, PRESSO global test, Cochran's Q
- **Leave-one-out:** SNP-level (S10–S15) and cohort-level (Supplementary Table S5)
- **Trio-based:** Maternal/fetal/paternal effect decomposition (21 outcomes; Figure 3, Table S4)
- **Multiple testing:** Benjamini–Hochberg FDR across 30 outcomes (q < 0.05)
- **Software:** TwoSampleMR v0.5.6, MRPRESSO v1.0.0, R v4.3.2
- **Very preterm birth** defined as < 34 weeks of gestation

---

## Interpretation

> *A robust causal signal for placenta praevia alongside directionally consistent estimates across the placental disorders domain suggests that abnormal implantation and placentation may constitute the primary mechanism linking endometriosis liability to adverse pregnancy outcomes. The absence of robust causal signals for hypertensive disorders, gestational diabetes, postpartum haemorrhage, fetal growth restriction, and most other obstetric outcomes suggests that previously reported observational associations for these outcomes may largely reflect confounding or clinical management patterns rather than direct biological effects of endometriosis. These findings support targeted attention to placental localisation rather than a generalised intensification of obstetric surveillance, and highlight the need to better characterise biological pathways linking endometriosis to impaired implantation and decidualisation.*

---

## Citation

```
Vibert J, Cheng TS, Magnus MC, Aiton E, Kutalik Z, Baud D, Lawlor DA,
Borges MC, Pluchino N. Genetic liability to endometriosis and pregnancy
outcomes: a two-sample Mendelian randomization study with maternal–fetal
effect decomposition. [Journal] [Year].

GitHub: https://github.com/jonasvibert/endoMR-PREG
```

---

## Data Sources

- Endometriosis GWAS: [Rahmioglu et al. 2023](https://doi.org/10.1038/s41588-023-01323-z)
- MR-PREG: [McBride et al. 2026](https://doi.org/10.1136/bmjopen-2025-103753)
- FinnGen R12: [Kurki et al. 2023](https://doi.org/10.1038/s41586-022-05473-8)
- PPH GWAS: [Westergaard et al. 2024](https://doi.org/10.1038/s41588-024-01839-y)

---

## Data Access

Analytical R code is available in this repository. Access to the underlying GWAS summary statistics requires individual applications to each data source:

- **MR-PREG** — Data can only be used for research covered by data agreements with contributing studies. ALSPAC data are available on request under proposal number B3844 ([apply here](http://www.bristol.ac.uk/alspac/researchers/access/)). BiB data available upon request ([borninbradford.nhs.uk](https://borninbradford.nhs.uk/research/how-to-access-data/)). MoBa data via [Helsedata](https://helsedata.no). UK Biobank data via the [Access Management System](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access).
- **Endometriosis GWAS** (Rahmioglu et al.) — Summary statistics (excluding 23andMe) available from the [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/) (Accession GCST90205183). 23andMe data available upon request (dataset-request@23andme.com).
- **PPH GWAS** (Westergaard et al.) — Summary statistics deposited at [deCODE](https://www.decode.com/summarydata).
- **FinnGen R12** — Publicly available; see [FinnGen documentation](https://finngen.gitbook.io/documentation/).

---

## Contact

**Jonas Vibert, MD** — jonas.vibert@chuv.ch
Department of Obstetrics and Gynecology, Lausanne University Hospital (CHUV)
Rue du Bugnon 21, 1011 Lausanne, Switzerland

---

## Acknowledgements

We are extremely grateful to all participants and families who contributed to the cohort studies included in this work (ALSPAC, Born in Bradford, MoBa, UK Biobank, FinnGen), as well as to the investigators, clinicians, and data collection teams involved. We thank the MR-PREG consortium participants and investigators for providing access to summary-level GWAS data across pregnancy and perinatal outcomes.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors. TSC, MCB, EA, and DAL are members of the MRC Integrative Epidemiology Unit at the University of Bristol (MC_UU_00032/5). MCB is also supported by the Leducq Foundation. MCM is supported by the Research Council of Norway through its Centres of Excellence funding scheme (project No 262700) and the research project "Endometriosis and adenomyosis throughout the life-course" (project No 351058).

---

## Declaration of Interest

MCB and DAL received funding from Novartis for unrelated research. JV, TSC, MCM, EA, ZK, DB, and NP declare no competing interests.

---

## License

MIT License — see LICENSE file for details.

---

*Last updated: May 2026*
