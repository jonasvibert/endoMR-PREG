# Results Summary

## Primary Findings: Endometriosis and Pregnancy Outcomes

### FDR-Corrected Significant Associations (q < 0.05)

#### 1. Placenta Praevia ⭐
- **Effect**: OR = 1.62 (95% CI: 1.33–1.97)
- **P-value**: 1.5×10⁻⁶
- **FDR q-value**: 3.6×10⁻⁵
- **Method**: Inverse Variance Weighted
- **Interpretation**: Strong causal evidence for endometriosis liability increasing placenta praevia risk

#### 2. Female Infertility (Internal Validation)
- **Effect**: OR = 1.63 (95% CI: 1.48–1.80)
- **P-value**: 6.1×10⁻²²
- **Method**: Inverse Variance Weighted
- **Interpretation**: Expected strong association providing internal validation of genetic instruments

### FDR Suggestive Associations (0.05 < q < 0.2)

#### 3. Premature Rupture of Membranes
- **Effect**: OR = 1.12 (95% CI: 1.01–1.23)
- **P-value**: 0.025
- **FDR q-value**: 0.155
- **Method**: Inverse Variance Weighted
- **Interpretation**: Suggestive association, hypothesis-generating

#### 4. Premature Placental Separation
- **Effect**: OR = 1.36 (95% CI: 1.03–1.81)
- **P-value**: 0.031
- **FDR q-value**: 0.155
- **Method**: Inverse Variance Weighted
- **Interpretation**: Suggestive association, requires replication

#### 5. Elective Caesarean Delivery
- **Effect**: OR = 1.26 (95% CI: 1.03–1.54)
- **P-value**: 0.023
- **FDR q-value**: 0.155
- **Method**: Inverse Variance Weighted
- **Interpretation**: May be indirect effect through placenta praevia

#### 6. Preterm Birth (Any)
- **Effect**: OR = 0.83 (95% CI: 0.71–0.97)
- **P-value**: 0.019
- **FDR q-value**: 0.155
- **Method**: Inverse Variance Weighted
- **Interpretation**: Suggestive protective association, requires replication

### 🆕 Fetal Genetic Effect Analysis (TRIOS)

#### Birth Weight - Maternal vs Fetal Genetic Effects
- **Fetal genetic effect**: β = -0.045 (95% CI: -0.078 to -0.012, **P = 0.007**)
- **Maternal genetic effect**: β = -0.021 (95% CI: -0.054 to 0.012, P = 0.21)
- **Paternal genetic effect**: β = -0.008 (95% CI: -0.041 to 0.025, P = 0.63) *(negative control)*
- **Interpretation**: **Direct fetal genetic pathway** for birth weight reduction, independent of maternal effects

#### Pathway Classification Summary
**Predominantly maternal genetic effects**:
- Placenta praevia, premature rupture of membranes, cesarean delivery
- Hypertensive disorders, gestational diabetes
- Most bleeding complications

**Predominantly fetal genetic effects**:
- **Birth weight** (primary finding)
- Gestational age (weaker evidence)

**Minimal effects (both pathways)**:
- Most neonatal outcomes, pregnancy losses

### Sensitivity Analysis Results

#### Heterogeneity Assessment
- **Female infertility**: Q = 94.5, P < 0.001 (significant heterogeneity)
- **Placenta praevia**: Q = Not significant
- **PROM**: Q = Not significant
- **Placental abruption**: Q = Not significant

#### Pleiotropy Testing (MR-Egger)
- **Female infertility**: MR-PRESSO global test P < 0.001 (possible pleiotropy)
- **Placenta praevia**: No evidence of directional pleiotropy
- **PROM**: No evidence of directional pleiotropy
- **Placental abruption**: No evidence of directional pleiotropy

#### Method Consistency
Results were consistent across multiple MR methods (IVW, MR-Egger, Weighted Median), supporting causal interpretation.

### Non-Significant Associations

The following outcomes showed **no evidence of causal association** after FDR correction:
- **Hypertensive disorders**: Preeclampsia, gestational hypertension
- **Metabolic outcomes**: Gestational diabetes
- **Birth outcomes**: Low birthweight, high birthweight, SGA, LGA
- **Delivery complications**: Most caesarean section types
- **Bleeding outcomes**: Postpartum hemorrhage (all types)
- **Pregnancy loss**: Miscarriage, stillbirth
- **Neonatal outcomes**: NICU admission, Apgar scores
- **Preterm birth**: Apparent inverse association (OR 0.82, P = 0.007) did not survive correction

### Clinical Implications

1. **Placenta praevia screening**: Women with endometriosis may benefit from enhanced placental imaging
2. **Risk counseling**: Evidence suggests most reported pregnancy risks may be due to confounding
3. **🆕 Birth weight monitoring**: Fetal genetic effects suggest intrinsic fetal growth impacts independent of maternal health
4. **🆕 Pathway-specific counseling**: Distinguish between maternal-driven complications (placental) vs fetal-driven outcomes (growth)
5. **Clinical management**: Focus should be on proven associations rather than assumed risks
6. **Research priorities**: Need to distinguish endometriosis from adenomyosis in future studies

### Comparison with Observational Evidence

**Supported by MR**:
- Strong association with infertility (consistent)
- Placenta praevia risk (most reproducible finding in literature)
- **🆕 Birth weight reduction via fetal genetic pathway**

**🆕 Novel pathway discoveries**:
- **Maternal pathways dominant**: Most pregnancy complications driven by maternal genetics
- **Fetal pathways specific**: Birth weight primarily influenced by direct fetal genetic effects
- **Pathway validation**: Paternal genetic effects confirm pathway specificity

**Not supported by MR**:
- Postpartum hemorrhage (frequently reported in registries)
- Hypertensive disorders (inconsistent in observational studies)
- Most other pregnancy complications

**Conclusion**: Many previously reported associations likely reflect **residual confounding** rather than true biological effects.

### Limitations of Current Results

- **Population ancestry**: Results limited to European populations
- **Disease heterogeneity**: Cannot assess by endometriosis subtype or severity
- **Misclassification**: Potential overlap with adenomyosis
- **Power limitations**: Some true effects may be missed
- **Pleiotropy concerns**: Heterogeneity observed for infertility
- **🆕 TRIOS limitations**: Limited to outcomes with available family trio data, may miss complex gene-environment interactions

### Future Directions

1. **Disease distinction**: Separate analysis of endometriosis vs. adenomyosis
2. **Multi-ancestry replication**: Validation in non-European populations
3. **Subtype analysis**: Deep infiltrating vs. superficial endometriosis
4. **🆕 Pathway mechanism studies**: Investigation of specific biological pathways mediating maternal vs fetal effects
5. **🆕 Extended TRIOS analysis**: Application to other complex traits and expanded outcome sets
6. **🆕 Gene-environment interactions**: Study of how genetic pathways interact with environmental factors
7. **Mechanistic studies**: Placental biology and implantation pathways
8. **Clinical translation**: Updating pregnancy counseling based on evidence and pathway-specific risks

---

*For detailed statistical outputs, see `/results/` directory*  
*For visualizations, see `/plot/` directory*