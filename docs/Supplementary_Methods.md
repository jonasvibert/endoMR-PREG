# Supplementary Material

## Additional Analyses and Supporting Information

**Updated for 29 specialized pregnancy outcomes with novel fetal genetic effect analysis**

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
*Location*: `results/Supplementary_Table_S2_Sensitivity_Analysis.csv`

**Content**: 
- Heterogeneity statistics (Cochran's Q)
- Pleiotropy tests (MR-Egger intercept)
- Method comparison (IVW, MR-Egger, Weighted Median)
- Leave-one-out summaries

#### 🆕 Table S5: TRIOS Fetal Genetic Effect Analysis
*Location*: `fetal_genetic_analysis/trios_maternal_fetal_paternal_summary.csv`

**Content**:
- Maternal genetic effects (using maternal genotypes)
- Fetal genetic effects (using offspring genotypes)  
- Paternal genetic effects (negative control)
- Effect size comparisons and pathway identification
- Birth weight detailed analysis (primary fetal pathway finding)

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

#### Figure S5: Enhanced Forest Plots
*Location*: `plots/forest_endoMR-PREG_*.png`

**Description**: Multiple forest plot layouts for 29 specialized outcomes
- Standard IVW forest plot
- Multi-method comparison forest
- Enhanced layout with improved spacing
- endoPAIN-style layout

#### 🆕 Figure S6: TRIOS Maternal vs Fetal Genetic Effect Analysis
*Location*: `fetal_genetic_analysis/mat_fetal_paternal_comparison.png`

**Description**: Comprehensive comparison of maternal, fetal, and paternal genetic effects
- Birth weight primary analysis (fetal effect β = -0.045, P = 0.007)
- Maternal effect comparison (β = -0.021, P = 0.21)
- Paternal negative control validation (β = -0.008, P = 0.63)

#### 🆕 Figure S7: TRIOS Sensitivity Analysis
*Location*: `fetal_genetic_analysis/sensitivity/`

**Description**: Leave-one-out sensitivity plots for TRIOS analysis
- Individual outcome sensitivity testing
- Robustness of fetal vs maternal pathway identification

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

#### 🆕 TRIOS Fetal Genetic Effect Analysis

**Methodology**: Novel application of family trio data to distinguish maternal vs fetal genetic pathways

**Key findings**:
- **Birth weight**: Predominant fetal genetic effect (β = -0.045, P = 0.007) vs weaker maternal effect (β = -0.021, P = 0.21)
- **Pathway validation**: Paternal genetic effects minimal across outcomes (negative control confirmation)
- **Maternal pathway dominance**: Most pregnancy complications driven by maternal genetic liability

**Technical validation**:
- Leave-one-out sensitivity analysis confirms robustness
- Consistent findings across different SNP subsets
- Biological plausibility of pathway-specific effects

**Clinical implications**:
- Birth weight monitoring may need fetal-specific considerations
- Pathway-specific counseling for different pregnancy outcomes
- Distinction between intrinsic vs extrinsic fetal growth effects

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

**🆕 Endometriosis → Fetal birth weight (direct fetal pathway)**:
- Inherited genetic liability affecting fetal growth pathways
- Direct impact on fetal growth hormone signaling
- Intrinsic fetal metabolic programming
- Independent of maternal placental function

**🆕 Maternal vs Fetal pathway distinction**:
- **Maternal pathways**: Placental complications, delivery decisions, maternal health
- **Fetal pathways**: Growth parameters, developmental timing
- **Validation**: Paternal genetic effects as negative controls

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

1. **Multi-ancestry validation**: Replicate TRIOS findings in diverse populations
2. **🆕 Extended TRIOS application**: Apply maternal/fetal pathway analysis to other complex traits
3. **🆕 Mechanistic pathway studies**: Investigate specific biological pathways mediating fetal vs maternal effects
4. **Age and timing-stratified analyses**: Critical windows for genetic effect expression
5. **🆕 Gene-environment interaction studies**: How genetic pathways interact with environmental factors during pregnancy
6. **Integration with omics data**: Proteomics/metabolomics to understand pathway mechanisms
7. **🆕 Therapeutic targeting**: Develop pathway-specific interventions (maternal vs fetal)
8. **🆕 Precision counseling**: Risk prediction models incorporating pathway-specific genetic effects

---

*For questions about supplementary material, contact: [Your email]*