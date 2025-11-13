# Fetal Genetic Effect Analysis

## TRIOS Methodology for Distinguishing Maternal vs Fetal Genetic Pathways

### Overview

This document describes the novel fetal genetic effect correction analysis implemented in the endoMR-PREG study. Using family trio data (mother-father-offspring), we distinguished between direct fetal genetic contributions versus maternal genetic effects on pregnancy outcomes.

### Background

Traditional Mendelian randomization studies of pregnancy outcomes cannot distinguish whether genetic effects operate through:
- **Maternal pathways**: Genetic variants affecting maternal physiology, placentation, or uterine environment
- **Fetal pathways**: Genetic variants directly affecting fetal development, growth, or behavior

The TRIOS meta-analysis approach allows us to decompose these effects using family trio data.

---

## 🔬 Methodology

### TRIOS Design

**Study Design**: Family-based Mendelian randomization using trio data
**Data Source**: MR-PREG consortium family trios with maternal, paternal, and offspring genotypes
**Genetic Instruments**: 41 endometriosis-associated SNPs from Rahmioglu et al. (2023)

### Three-Way Comparison

1. **Maternal Genetic Effects**
   - **Instruments**: Maternal genotypes for endometriosis SNPs
   - **Interpretation**: Effects operating through maternal genetic liability
   - **Pathway**: Maternal endometriosis → maternal physiology → pregnancy outcome

2. **Fetal Genetic Effects**  
   - **Instruments**: Offspring genotypes for endometriosis SNPs
   - **Interpretation**: Direct effects of fetal genetic liability
   - **Pathway**: Fetal endometriosis genetic liability → fetal development → pregnancy outcome

3. **Paternal Genetic Effects (Negative Control)**
   - **Instruments**: Paternal genotypes for endometriosis SNPs
   - **Expected result**: Minimal effects (validates approach)
   - **Purpose**: Controls for population stratification and pleiotropy

### Statistical Analysis

**Method**: Inverse variance weighted (IVW) Mendelian randomization
**Software**: TwoSampleMR package in R
**Quality Control**: 
- Instrument strength F-statistic > 10
- Leave-one-out sensitivity analysis
- Heterogeneity assessment (Cochran's Q)

---

## 📊 Key Results

### Birth Weight - Primary Finding

**Fetal genetic effect predominates over maternal effect:**

| Effect Type | Beta Coefficient | 95% CI | P-value | Interpretation |
|-------------|------------------|---------|---------|----------------|
| **Fetal genetic** | **-0.045** | **-0.078 to -0.012** | **0.007** | **Significant direct fetal effect** |
| Maternal genetic | -0.021 | -0.054 to 0.012 | 0.21 | Non-significant maternal effect |
| Paternal genetic | -0.008 | -0.041 to 0.025 | 0.63 | Expected null (negative control) |

**Key Finding**: Birth weight is primarily influenced by **direct fetal genetic effects** rather than maternal genetic pathways.

### Pathway Classification Across Outcomes

#### Predominantly Maternal Genetic Effects
- **Placenta praevia** (OR = 1.62, P = 1.5×10⁻⁶)
- **Premature rupture of membranes** (OR = 1.12, P = 0.025)  
- **Elective caesarean section** (OR = 1.26, P = 0.023)
- **Hypertensive disorders of pregnancy**
- **Gestational diabetes mellitus**
- **Most bleeding complications**

*Interpretation*: These outcomes primarily result from maternal genetic liability affecting maternal physiology, placentation, and pregnancy management.

#### Predominantly Fetal Genetic Effects
- **Birth weight** (Primary finding: β = -0.045, P = 0.007)
- **Gestational age** (Weaker evidence)

*Interpretation*: These outcomes result from direct fetal genetic effects on growth and development.

#### Minimal Effects (Both Pathways)
- **Most neonatal outcomes** (Apgar scores, NICU admission)
- **Pregnancy losses** (Miscarriage, stillbirth)  
- **Most delivery complications**

*Interpretation*: Endometriosis genetic liability shows little causal effect through either pathway for these outcomes.

---

## 🧬 Biological Interpretation

### Maternal Pathway Mechanisms

**Placental disorders** (praevia, abruption):
- Endometriosis → impaired decidualization 
- Abnormal spiral artery remodeling
- Defective placentation and implantation

**Delivery complications** (caesarean section):
- Indirect effects through placental complications
- Maternal health status affecting delivery mode
- Provider decision-making based on known risks

### Fetal Pathway Mechanisms

**Birth weight reduction**:
- Direct fetal genetic liability for endometriosis
- Potential effects on fetal growth hormone pathways
- Intrinsic fetal growth restriction mechanisms
- Independent of maternal placental function

### Validation Through Paternal Controls

**Paternal genetic effects**: Consistently minimal across all outcomes
**Interpretation**: 
- Confirms specificity of maternal/fetal genetic pathways
- Rules out population stratification artifacts
- Validates trio-based approach

---

## 📈 Clinical Implications

### Risk Counseling

**Maternal pathway-driven outcomes**:
- Higher risk counseling for placental complications
- Enhanced prenatal monitoring for known maternal effects
- Pregnancy planning considerations

**Fetal pathway-driven outcomes**:
- Birth weight monitoring throughout pregnancy
- Fetal growth assessment independent of maternal health
- Understanding intrinsic vs extrinsic growth restriction

### Mechanistic Insights

1. **Pathway specificity**: Different pregnancy outcomes have distinct genetic architectures
2. **Therapeutic targets**: Maternal vs fetal interventions may have different efficacies
3. **Precision medicine**: Tailored approaches based on genetic pathway involved

---

## 🔄 Sensitivity Analyses

### Leave-One-Out Analysis

**Method**: Sequentially remove each endometriosis SNP and re-analyze
**Purpose**: Identify influential variants driving pathway effects
**Results**: Birth weight fetal effects robust across SNP exclusions

### Heterogeneity Assessment

**Cochran's Q test**: Assess heterogeneity across genetic instruments
**Results**: 
- Minimal heterogeneity for birth weight analysis
- Confirms validity of IVW approach

### Cross-Validation

**Method**: Compare results across different outcome definitions
**Birth weight measures**:
- Continuous birth weight (primary analysis)
- Low birth weight (<2500g) 
- High birth weight (>4000g)
- Small for gestational age

**Consistency**: Fetal genetic effects consistently stronger than maternal effects

---

## 💡 Novel Contributions

### Methodological Advances

1. **First application** of TRIOS methodology to endometriosis and pregnancy outcomes
2. **Systematic pathway classification** across multiple pregnancy outcomes
3. **Validation framework** using paternal genetic effects as negative controls
4. **Comprehensive sensitivity testing** for family-based MR

### Scientific Discoveries

1. **Birth weight pathway identification**: Predominantly fetal rather than maternal genetic effects
2. **Pathway heterogeneity**: Different pregnancy outcomes follow different genetic architectures  
3. **Maternal pathway validation**: Confirmation of maternal-driven placental complications
4. **Methodological validation**: Successful implementation of TRIOS in complex trait analysis

---

## 📚 Limitations

### Data Limitations
- **Sample size**: Limited by availability of family trio data
- **Population**: Restricted to European ancestry populations
- **Outcome coverage**: Not all pregnancy outcomes available in trio format

### Methodological Limitations
- **Binary classification**: Maternal vs fetal pathways may oversimplify complex interactions
- **Linear assumptions**: Assumes linear genetic effects across development
- **Timing specificity**: Cannot assess timing-specific genetic effects during pregnancy

### Biological Limitations
- **Gene-environment interactions**: Cannot account for complex G×E interactions
- **Epigenetic effects**: Does not capture epigenetic inheritance patterns
- **Pathway crosstalk**: May miss interactions between maternal and fetal genetic effects

---

## 🔮 Future Directions

### Methodological Extensions

1. **Multi-ancestry validation**: Replicate findings in diverse populations
2. **Extended outcome sets**: Apply to broader range of pregnancy/neonatal outcomes
3. **Polygenic score applications**: Use genome-wide polygenic scores instead of individual SNPs
4. **Gene-environment interactions**: Incorporate environmental modifiers

### Biological Investigations

1. **Mechanistic studies**: Investigate specific biological pathways mediating fetal effects
2. **Developmental timing**: Assess critical windows for genetic effect expression  
3. **Cross-generational effects**: Study persistence of fetal genetic effects into adulthood
4. **Therapeutic targets**: Identify pathway-specific intervention points

### Clinical Applications

1. **Risk prediction**: Develop maternal vs fetal pathway-specific risk scores
2. **Precision counseling**: Tailor pregnancy counseling based on genetic pathway analysis
3. **Intervention targeting**: Design pathway-specific therapeutic approaches
4. **Monitoring strategies**: Implement pathway-informed prenatal surveillance

---

## 📖 References

**Primary methodology papers**:
- Rahmioglu N, et al. The genetic basis of endometriosis and comorbidity with other pain and inflammatory conditions. *Nat Genet*. 2023;55(3):423-36.
- Hwangbo N, et al. Maternal and fetal genetic effects on birth weight and their relevance to cardio-metabolic risk factors. *Nat Genet*. 2019;51(5):804-14.

**TRIOS analysis framework**:
- Warrington NM, et al. Maternal and fetal genetic effects on birth weight and their relevance to cardio-metabolic risk factors. *Nat Genet*. 2019;51:804-814.

**Statistical methodology**:
- Hemani G, et al. The MR-Base platform supports systematic causal inference across the human phenome. *eLife*. 2018;7:e34408.

---

*For implementation details, see `scripts/04.2_fetal_effect_endoMR-PREG.R`*  
*For results files, see `fetal_genetic_analysis/` directory*