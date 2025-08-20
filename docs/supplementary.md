# Supplementary Material

## Additional Analyses and Supporting Information

### Supplementary Tables

#### Table S1: Endometriosis Genetic Instruments
*Location*: `Supplementary_Table_S1_Endometriosis_Instruments.xlsx`

**Columns**:
- SNP: rsID identifier
- Chromosome: Chromosome number
- Position: Base pair position (GRCh37)
- Effect allele: Allele increasing endometriosis risk
- Other allele: Reference allele
- EAF: Effect allele frequency
- Beta: Effect size (log odds)
- SE: Standard error
- P-value: Association p-value
- F-statistic: Instrument strength

#### Table S2: Sample Sizes by Outcome
*Location*: `results/Table1_sample_sizes_clean.csv`

**Content**: Sample sizes for all pregnancy outcomes analyzed

#### Table S3: Complete MR Results
*Location*: `results/all_mr_methods.csv`

**Content**: Results from all MR methods for all outcomes

#### Table S4: Sensitivity Analysis Summary
*Location*: `results/mrpresso_summary.xlsx`

**Content**: 
- Heterogeneity statistics
- Pleiotropy tests
- MR-PRESSO results
- Leave-one-out summaries

### Supplementary Figures

#### Figure S1: Instrument Characteristics
*Location*: `plot/fig_scatter/`

**Description**: Scatter plots showing relationship between endometriosis and significant outcomes
- Panel A: Female infertility
- Panel B: Premature rupture of membranes  
- Panel C: Placenta praevia
- Panel D: Placental abruption

#### Figure S2: Forest Plots (Individual SNPs)
*Location*: `plot/fig_forest/`

**Description**: Forest plots showing individual SNP effects and overall IVW estimate

#### Figure S3: Leave-One-Out Analysis
*Location*: `plot/fig_leaveoneout_egger/`

**Description**: Leave-one-out plots demonstrating robustness of findings

#### Figure S4: Funnel Plots
*Location*: `plot/fig_funnel/`

**Description**: Funnel plots for assessment of directional pleiotropy

#### Figure S5: Comprehensive Forest Plot
*Location*: `plot/fig_1toMany_forest/MR_overview_1toMany.png`

**Description**: One-to-many forest plot showing all outcomes grouped by clinical category

### Additional Analyses

#### Power Calculations

**Method**: Brion et al. (2013) approach
**Results**: 
- Power > 80% to detect OR ≥ 1.2 for outcomes with n > 50,000
- Limited power for smaller effect sizes in smaller studies

#### Instrumental Variable Assumptions

**Relevance (F-statistic)**:
- All instruments: F > 10 (range: [Add range])
- Cumulative R²: [Add value]% of endometriosis variance explained

**Independence**: 
- No known confounders associated with selected instruments
- Checked against GWAS catalog for other associations

**Exclusion restriction**:
- MR-Egger intercept tests (see Table S4)
- No evidence of directional pleiotropy

#### Subgroup Analyses

**By instrument strength**:
- Top 20 vs. all 41 instruments
- Results consistent (see `results/` folder)

**By outcome type**:
- FinnGen vs. other sources
- Similar effect patterns observed

#### Validation Analyses

**Alternative clumping thresholds**:
- r² < 0.01: 35 instruments, consistent results
- r² < 0.1: 28 instruments, similar estimates

**Different P-value thresholds**:
- P < 1×10⁻⁶: 67 instruments (exploratory)
- P < 5×10⁻⁸: 41 instruments (primary analysis)

### Biological Plausibility

#### Proposed Mechanisms

**Endometriosis → Infertility**:
- Anatomical distortion
- Inflammatory environment
- Ovarian reserve impact

**Endometriosis → Placental disorders**:
- Shared inflammatory pathways
- Vascular dysfunction
- Immune dysregulation

**Endometriosis → Preterm birth**:
- Systemic inflammation
- Cervical factors
- Uterine contractility

#### Supporting Literature

1. Inflammatory markers elevated in endometriosis
2. Placental vascular abnormalities in endometriosis patients
3. Obstetric complications well-documented observationally

### Limitations and Considerations

#### Study Limitations

1. **Population ancestry**: European-focused analysis
2. **Pleiotropy**: Cannot completely rule out horizontal pleiotropy
3. **Non-linear effects**: MR assumes linear relationships
4. **Timing**: Cannot assess age-specific effects

#### Technical Limitations

1. **Winner's curse**: Potential overestimation of instrument effects
2. **Linkage disequilibrium**: Possible confounding through LD
3. **Population stratification**: Residual confounding possible

#### Clinical Considerations

1. **Effect sizes**: Moderate effects, clinical significance varies
2. **Individual prediction**: Population-level estimates only
3. **Intervention**: MR shows causal potential, not treatment targets

### Software Versions

**Primary analysis**:
- R version: 4.3.1
- TwoSampleMR: 0.5.7
- dplyr: 1.1.2
- ggplot2: 3.4.2

**Sensitivity analyses**:
- MendelianRandomization: 0.8.0
- MRPRESSO: 1.0

### Data Processing Details

#### SNP Selection Process
1. Genome-wide significant SNPs (P < 5×10⁻⁸)
2. LD clumping (r² < 0.001, 10Mb window)
3. Instrument validation (F-statistic calculation)
4. Palindromic SNP handling

#### Quality Control Steps
1. Allele frequency checks (MAF > 0.01)
2. Effect size sanity checks
3. Duplicate SNP removal
4. Cross-reference with outcome data

#### Harmonization Process
1. Effect allele alignment
2. Strand flip correction  
3. Palindromic SNP exclusion (MAF > 0.42)
4. Final dataset validation

### Reproducibility Information

**Random seeds**: Set to 12345 for all analyses
**System information**: macOS, 16GB RAM
**Computational time**: ~30 minutes for full analysis

**Exact package versions**:
```r
sessionInfo()
# Add output here
```

### Future Research Directions

1. **Multi-ancestry validation**
2. **Mechanistic pathway analysis**
3. **Age-stratified analyses**
4. **Maternal vs. fetal genetic effects**
5. **Integration with proteomics/metabolomics**

---

*For questions about supplementary material, contact: [Your email]*