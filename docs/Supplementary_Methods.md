Supplementary Methods

This two-sample Mendelian randomization (MR) study was designed to evaluate whether genetic liability to endometriosis has a causal effect on a range of adverse pregnancy and perinatal outcomes. The analysis pipeline was implemented in R using a set of modular scripts stored in the endoMR-PREG GitHub repository (https://github.com/jonasvibert/endoMR-PREG/tree/3c6357c4e89748dc67e96dd292d10615feffbb93/scripts):
01_select_instruments_endoMR-PREG.R
02_prepare_outcomes_endoMR-PREG.R
03_harmonise_data_endoMR-PREG.R
04_main_analyses_endoMR-PREG.R
04.2_fetal_effect_endoMR-PREG.R
05_sensitivity_analyses_endoMR-PREG.R
06_tables_endoMR-PREG.R
07_plots_endoMR-PREG.R)
Together, these scripts reproduce the full workflow, from instrument selection to figure generation.

Exposure: genetic liability to endometriosis

Genetic instruments for endometriosis were derived from the largest published genome-wide association study (GWAS) meta-analysis of clinically confirmed endometriosis by Rahmioglu et al. (Nature Genetics 2023) (1). This meta-analysis included approximately 58,000 surgically and/or clinically diagnosed cases and about 733,000 female controls, predominantly of European ancestry, contributed by the International Endometriosis Genetics Consortium, 23andMe, FinnGen and other cohorts. In 01_select_instruments_endoMR-PREG.R, we imported the European-ancestry summary statistics from the meta-analysis that included 23andMe, and restricted all subsequent analyses to these data to maintain ancestry compatibility with the outcome GWAS.

Instrument selection followed a two-step procedure. First, we extracted all variants reaching conventional genome-wide significance for endometriosis (P < 5 × 10^-8) together with their rsIDs, effect and non-effect alleles, effect sizes (beta), standard errors (SE), effect allele frequencies (EAF) and, where available, sample size (N). Some variants were reported by Rahmioglu et al. in chromosome:position (chr:pos) format rather than rsID. For these, 01_select_instruments_endoMR-PREG.R mapped coordinates to rsIDs using the 1000 Genomes Project Phase 3 European reference panel (2). When multiple rsIDs were present at a given position, we retained the allele-specific identifier consistent with the effect allele reported in the GWAS. Variants with ambiguous mapping or unresolved alleles were excluded at this stage.

To ensure independence between instruments, we performed linkage disequilibrium (LD) clumping using PLINK 1.9 with 1000 Genomes Phase 3 Europeans as the reference population. The script 01_select_instruments_endoMR-PREG.R applied an r^2 threshold of 0.001 within a 10,000 kb physical window (parameters clump_r2 = 0.001, clump_kb = 10,000), using 5 × 10^-8 as the primary P-value threshold (clump_p1) and 1 as the secondary threshold (clump_p2). After clumping, 41 independent single nucleotide polymorphisms (SNPs) were retained as instruments for genetic liability to endometriosis.

Instrument strength was systematically evaluated in the same script. For each SNP, we calculated a simple F-statistic as:

F_simple = (beta / SE)^2

We then estimated the variance explained in endometriosis liability (R^2) using the standard formula:

R^2 = [ 2 × EAF × (1 − EAF) × beta^2 ] / [ 2 × EAF × (1 − EAF) × beta^2 + SE^2 × N ]

From this, we derived an exact per-SNP F-statistic:

F_exact = R^2 × (N − 2) / (1 − R^2)

Summing R^2 across all independent SNPs (assuming negligible LD after clumping) gave a total variance explained of approximately 5.6%. The mean per-SNP F-statistic was around 279, substantially above the conventional threshold of 10 and indicating a very low risk of weak instrument bias. Full instrument characteristics (rsID, genomic position, alleles, beta, SE, EAF, R^2, F) are provided in Supplementary Table S1.

Outcome GWAS and phenotype definition

We evaluated 29 maternal and perinatal outcomes. These outcomes were chosen based on prior epidemiological evidence linking endometriosis with obstetric complications, clinical relevance for maternal–fetal health, and the availability of adequately powered GWAS in European-ancestry populations. We drew outcome summary statistics from three sources: the Mendelian Randomization in Pregnancy (MR-PREG) collaboration (3), FinnGen release 12 (4), and a large postpartum haemorrhage (PPH) meta-analysis (5). Extraction, cleaning and initial formatting of these GWAS were performed in 02_prepare_outcomes_endoMR-PREG.R.
MR-PREG is an international consortium that harmonises and meta-analyses GWAS of pregnancy and perinatal outcomes across several large European birth cohorts, including ALSPAC, Born in Bradford, MoBa and UK Biobank, following a standardised pipeline described in detail in McBride et al. (2025) (3). Briefly, each cohort applied standard sample-level quality control (removal of individuals with low call rate, sex discordance, outlying heterozygosity or cryptic relatedness) and restricted analyses to participants of European genetic ancestry. Variant-level quality control included filters for call rate, minor allele frequency, Hardy–Weinberg equilibrium and imputation quality, as described in McBride et al. (3). Association analyses were performed using logistic or linear regression, or REGENIE for binary and continuous traits, with appropriate adjustment for age, ancestry principal components, batch effects and cohort-specific covariates, and rare outcomes were handled using Firth correction. Cohort-level results were meta-analysed using fixed-effects inverse-variance weighting.
From MR-PREG we used 21 outcomes provided in the consortium’s meta-analysed GWAS: gestational age (full sample and genotyped subsample), gestational diabetes, gestational hypertension, hypertensive disorders of pregnancy, high and low birthweight, small-for-gestational-age, birthweight z-score, large-for-gestational-age, preterm birth, very preterm birth, post-term birth, labour induction, premature rupture of membranes, Apgar score <7 at 1 and 5 minutes, NICU admission, pregnancy anaemia and preeclampsia. Caesarean section phenotypes (overall, elective and emergency), stillbirth and postpartum/peripartum depression were obtained through MR-PREG or collaborating consortia using the same QC and meta-analytic framework. For each outcome GWAS we extracted SNP-level beta, SE, effect allele, allele frequency and sample size directly from the published summary statistics.
To complement these outcomes with more detailed placental phenotypes, we used FinnGen release 12 (4). FinnGen GWAS follow a standardised QC and analysis framework described in Kurki et al. (2023), including array-based genotyping, imputation to the SISu v3 Finnish reference panel, strict sample and variant filtering, and association testing using REGENIE. From FinnGen, we extracted placenta praevia, placental disorders and placental abruption, defined through registry-based ICD codes.
Finally, bleeding-related outcomes were obtained from the GWAS meta-analysis by Westergaard et al. (2024) (5), which aggregates data from multiple Nordic biobanks and UK Biobank. We included antepartum bleeding, postpartum haemorrhage (PPH) overall, PPH due to uterine atony and PPH due to retained placenta as defined in that study. Reported odds ratios were converted to log-odds (beta = log(OR)); when standard errors were not available, they were derived from the reported P-values using SE = |beta| / z, where z is the standard normal quantile corresponding to P/2.
Overall, the 29 outcomes can be grouped into nine clinical domains: placental disorders (placenta praevia, placental disorders, abruption); hypertensive disorders of pregnancy (hypertensive disorders overall, gestational hypertension, preeclampsia); pregnancy timing (gestational age, preterm birth, very preterm birth, post-term birth); fetal growth and birthweight (high birthweight, low birthweight, small-for-gestational-age, birthweight z-score, large-for-gestational-age); labour and delivery complications (labour induction, premature rupture of membranes, caesarean section overall, elective and emergency caesarean); bleeding and haemorrhage (antepartum bleeding, PPH overall, PPH due to atony, PPH due to retained placenta); maternal metabolic/haematologic complications (gestational diabetes, pregnancy anaemia); maternal mental health (postpartum/peripartum depression); and neonatal condition at birth (low Apgar scores at 1 and 5 minutes, NICU admission, stillbirth). Detailed definitions, case/control counts and contributing cohorts are presented in Table 1 and Supplementary Table S2.

Harmonisation of exposure and outcome data

Harmonisation of endometriosis and outcome GWAS was carried out in 03_harmonise_data_endoMR-PREG.R using the TwoSampleMR package (version 0.5.6) with additional manual checks. For each outcome, we first matched the 41 endometriosis instrument SNPs to outcome SNPs by rsID. We then aligned alleles so that the effect allele was identical in exposure and outcome datasets, flipping the sign of the outcome beta and adjusting allele frequencies where necessary. Palindromic SNPs with intermediate allele frequencies (A/T or C/G variants with minor allele frequency around 0.5) were removed using the harmonise_data function with action = 2, because strand assignment is ambiguous in this range. Palindromic SNPs with clearly low or high MAF, for which strand could be reliably inferred, were retained. SNPs with unresolved strand issues, mismatched alleles or missing beta/SE in either dataset were excluded.

Because not all instrument SNPs survived quality control or were available in every outcome GWAS, the effective number of instruments contributing to each MR analysis ranged from 29 to 41. The exact SNP count per outcome is reported in Supplementary Table S2, alongside exposure and outcome sample sizes.

Two-sample MR analyses

Primary MR analyses were implemented in 04_main_analyses_endoMR-PREG.R. For each outcome, we used the inverse-variance weighted (IVW) estimator as the main causal effect measure (6). For each SNP i, we denote the SNP–endometriosis association as beta_Xi, the SNP–outcome association as beta_Yi, and the variance of beta_Yi as var_Yi. The IVW estimator can be written as:

beta_IVW = sum( beta_Xi × beta_Yi / var_Yi ) / sum( beta_Xi^2 / var_Yi )

This is equivalent to a weighted regression of beta_Yi on beta_Xi with no intercept and weights equal to 1 / var_Yi. For binary outcomes we exponentiated beta_IVW to obtain odds ratios (OR = exp(beta_IVW)) with 95% confidence intervals, whereas for continuous outcomes (gestational age, birthweight z-score) we reported the beta coefficient and corresponding 95% confidence intervals.

Sensitivity analyses and assessment of MR assumptions

Sensitivity analyses were conducted in 05_sensitivity_analyses_endoMR-PREG.R to evaluate the robustness of IVW estimates to potential violations of MR assumptions, particularly horizontal pleiotropy. First, we performed MR-Egger regression, modelling beta_Yi as:

beta_Yi = intercept + beta_MR_Egger × beta_Xi + error

The slope beta_MR_Egger provides a pleiotropy-adjusted causal estimate under the Instrument Strength Independent of Direct Effect (InSIDE) assumption, while the intercept tests for directional pleiotropy: a non-zero intercept suggests that instruments have, on average, direct effects on the outcome independent of endometriosis liability.

Second, we estimated causal effects using the weighted median estimator, which remains consistent if at least 50% of the total instrument weight is contributed by valid instruments. This estimator is based on the median of SNP-specific ratio estimates (beta_Yi / beta_Xi), weighted by the inverse variance of beta_Yi. Where informative, we additionally computed weighted mode estimates, which assume that the most common (modal) causal estimate across SNPs arises from valid instruments.

Third, we quantified heterogeneity among SNP-specific causal estimates using Cochran’s Q statistic. For each SNP we calculated its ratio estimate beta_i = beta_Yi / beta_Xi and weight w_i = beta_Xi^2 / var_Yi, and then computed:

Q = sum( w_i × (beta_i − beta_IVW)^2 )

A large Q statistic with a low P-value indicates heterogeneity and may suggest pleiotropy or violation of model assumptions.

Fourth, we applied the MR-PRESSO (Pleiotropy RESidual Sum and Outlier) method (MRPRESSO package v1.0.0). For each outcome, the global test examined whether the pattern of residuals was compatible with no horizontal pleiotropy. When the global test was significant, MR-PRESSO identified outlier SNPs and recomputed outlier-corrected IVW estimates after excluding them. We report global test P-values, outlier SNPs where present and corrected effect estimates.

Finally, we performed single-SNP MR and leave-one-SNP-out analyses. Single-SNP estimates examine the contribution of each instrument individually, while leave-one-out analyses re-run IVW MR after removing each SNP in turn to evaluate whether any single variant unduly drives results. Figures summarising these diagnostics (scatter plots, forest plots, funnel plots, leave-one-out plots) were generated in 07_plots_endoMR-PREG.R and are presented in Supplementary Figures S1–S7.

Maternal, fetal and paternal genetic effects

Some outcomes, such as hypertensive disorders of pregnancy, birthweight and preterm birth, may be influenced by both maternal and fetal genomes. Where available, we therefore leveraged trio-based genetic effect estimates from MR-PREG to distinguish maternal and fetal pathways. These analyses were implemented in 04.2_fetal_effect_endoMR-PREG.R.

MR-PREG provides effect estimates from models jointly including maternal, offspring (fetal) and paternal genotypes, derived from mother–father–child trios or extended parent–offspring structures, using regression models of the form 

expected outcome = beta_m × G_m + beta_o × G_o + beta_p × G_p + covariates

, where GmG_mGm​, GoG_oGo​ and GpG_pGp​ represent maternal, offspring and paternal genotypes. Maternal effects (βm\beta_mβm​) capture direct maternal genetic influence on the intrauterine environment, fetal effects (βo\beta_oβo​) reflect the fetus’s own genome acting on fetal growth and development, and paternal effects (βp\beta_pβp​) serve as a negative control, as paternal genotype should not directly affect maternal pregnancy complications.

Trio-based genetic associations were available for 21 of the 29 investigated outcomes, namely small-for-gestational-age, large-for-gestational-age, high birthweight (>4000 g), low birthweight (<2500 g), birthweight z-score, gestational age, post-term birth, preterm birth <37 weeks, preeclampsia, gestational hypertension, hypertensive disorders of pregnancy, gestational diabetes, pregnancy anemia, labour induction, emergency caesarean section, elective caesarean section, low Apgar score at 1 minute, low Apgar score at 5 minutes, neonatal intensive care unit admission, postpartum depression and premature rupture of membranes. 
For each of these outcomes, we performed three parallel MR analyses using the endometriosis SNP instruments: one using maternal effect estimates, one using fetal effect estimates and one using paternal effect estimates as a negative control. The same causal estimators were applied in each pathway (IVW, MR-Egger, weighted median, weighted mode), together with heterogeneity statistics, MR-Egger intercepts, MR-PRESSO global and outlier tests and leave-one-out analyses. Because not all trio-based SNP associations were available for all outcomes, the number of instruments contributing to each maternal, fetal, and paternal analysis varied accordingly, ranging from 27 to 34 SNPs (details in Supplementary Table S3). Results of these analyses, allowing direct comparison of maternal, fetal and paternal pathways, are reported in Supplementary Tables S3–S4 and Supplementary Figure S8.

Multiple testing correction

Because we evaluated causal estimates for 29 outcomes, we controlled for multiplicity using the false discovery rate (FDR). In 05_sensitivity_analyses_endoMR-PREG.R, we applied the Benjamini–Hochberg procedure to the two-sided IVW P-values across all outcomes using R’s p.adjust function in R with method = "fdr". The resulting FDR-adjusted P-values (q-values) reflect the expected proportion of false positives among declared significant findings. We considered associations with q < 0.05 as statistically significant after FDR correction, and we report both nominal P-values and q-values in Supplementary Table S2.

Tables, figures and reproducibility

The script 06_tables_endoMR-PREG.R assembled results from all previous steps into the tables presented in the main manuscript and supplementary material. It joins MR estimates, confidence intervals, heterogeneity statistics, MR-Egger intercepts, MR-PRESSO outputs, SNP counts and sample sizes into publication-ready tables (for example, for IVW estimates, sensitivity analyses and trio-based MR). 07_plots_endoMR-PREG.R uses the same underlying results to generate forest plots of IVW estimates across outcomes (including the main 29-outcome figure), scatter plots contrasting SNP–exposure and SNP–outcome associations, funnel plots for asymmetry, leave-one-out plots and trio comparison plots.

All analyses were conducted in R version 4.3.2. Key packages included TwoSampleMR (v0.5.6) for harmonisation and MR methods, MRPRESSO (v1.0.0) for pleiotropy assessment, data.table and dplyr for data handling, ggplot2 for visualisation, and readr, stringr, tidyr and here for import and project organisation. LD clumping was performed using PLINK 1.9. The full set of scripts (01–07) and a detailed README are available in the publicly accessible endoMR-PREG repository (https://github.com/jonasvibert/endoMR-PREG
), allowing complete reproducibility of the analysis pipeline.

Data availability and ethics

This study relied exclusively on de-identified summary-level GWAS data from previously approved studies and consortia. Endometriosis GWAS summary statistics were obtained from Rahmioglu et al. (Nature Genetics 2023) (1). MR-PREG outcome GWAS were accessed through the MR-PREG collaboration under its data-sharing policies (3). FinnGen summary statistics (release 12) are publicly available through the FinnGen results portal (4), and bleeding outcomes were obtained from Westergaard et al. (Nature Genetics 2024) (5). Because only aggregated summary statistics were used and no individual-level genetic or clinical data were accessed, the present analyses do not constitute human subjects research as defined by 45 CFR 46.102, and no additional institutional review board approval was required. All original contributing studies obtained appropriate ethical approvals and informed consent from participants.

References
1.	Rahmioglu N, Mortlock S, Ghiasi M, Møller PL, Stefansdottir L, Galarneau G, et al. The genetic basis of endometriosis and comorbidity with other pain and inflammatory conditions. Nat Genet. 2023 Mar;55(3):423–36. 
2.	The 1000 Genomes Project Consortium, Corresponding authors, Auton A, Abecasis GR, Steering committee, Altshuler DM, et al. A global reference for human genetic variation. Nature. 2015 Oct 1;526(7571):68–74. 
3.	McBride N, Clayton GL, Goncalves Soares A, Yang Q, Bond TA, Taylor A, et al. Cohort Profile: The Mendelian Randomization in Pregnancy (MR-PREG) collaboration - Improving evidence for prevention and treatment of adverse pregnancy and perinatal outcomes [Internet]. Epidemiology; 2025 [cited 2025 Aug 9]. Available from: http://medrxiv.org/lookup/doi/10.1101/2025.03.22.25324447 
4.	Kurki MI, Karjalainen J, Palta P, Sipilä TP, Kristiansson K, Donner KM, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature. 2023 Jan 19;613(7944):508–18. 
5.	Westergaard D, Steinthorsdottir V, Stefansdottir L, Rohde PD, Wu X, Geller F, et al. Genome-wide association meta-analysis identifies five loci associated with postpartum hemorrhage. Nat Genet. 2024 Aug;56(8):1597–603. 
6.	Burgess S, Davey Smith G, Davies NM, Dudbridge F, Gill D, Glymour MM, et al. Guidelines for performing Mendelian randomization investigations: update for summer 2023. Wellcome Open Res. 2023 Aug 4;4:186. 
