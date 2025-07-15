### 3. PLOTS  ###################################################
### Scatter plot
res <- mr(dat)
p1 <- mr_scatter_plot(res, dat)
length(p1)
ggsave(p1[[1]], file = "SCATTERPLOTRESULTS.pdf", width = 7, height = 7)
ggsave(p1[[2]], file = "SCATTERPLOTRESULTS2.pdf", width = 7, height = 7)
ggsave(p1[[33]], file = "SCATTERPLOTRESULTS2.pdf", width = 7, height = 7)
print(p1[[1]])

### Forest plot
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
p2[[37]]

### Leave-one-out plot
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

### Funnel plot
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p4[[1]]





### 7. PLOTTING ###################################################################
# 7.1 Scatter plot (IVW vs Egger vs modes)
p_scatter <- mr_scatter_plot(all_res, dat)
ggsave(plot = p_scatter[[1]], filename = file.path(results_dir, "scatter_ivw_egger.png"),
       width = 7, height = 7)

# 7.2 Forest plot for single-SNP
p_forest <- mr_forest_plot(single_res)
ggsave(plot = p_forest[[1]], filename = file.path(results_dir, "forest_singlesnp.png"),
       width = 7, height = 7)

# 7.3 Funnel plot for single-SNP
p_funnel <- mr_funnel_plot(single_res)
ggsave(plot = p_funnel[[1]], filename = file.path(results_dir, "funnel_singlesnp.png"),
       width = 7, height = 7)

# 7.4 Leave-one-out SNP plot
p_loo <- mr_leaveoneout_plot(loo_snp_res)
ggsave(plot = p_loo[[1]], filename = file.path(results_dir, "leaveoneout_plot.png"),
       width = 7, height = 7)

### 8. LEAVE-ONE-OUT BY COHORT ###################################################
# Requires a file 'stu_out_dat.txt' with columns: SNP, beta, se, effect_allele, other_allele, eaf, Phenotype, study
if (file.exists("stu_out_dat.txt")) {
  mr_data <- readr::read_delim("stu_out_dat.txt", delim = "\t", col_types = cols())
  cohorts  <- unique(mr_data$study)
  
  loo_cohort <- lapply(cohorts, function(coh) {
    df   <- dplyr::filter(mr_data, study != coh)
    exp  <- df %>% rename(beta.exposure = beta, se.exposure = se) %>% mutate(id.exposure = unique(df$Phenotype))
    out  <- df %>% rename(beta.outcome = beta, se.outcome = se, id.outcome = Phenotype)
    hmd  <- harmonise_data(exp, out, action = 2)
    res  <- mr(hmd, method_list = "mr_ivw")
    res$left_out_cohort <- coh
    res
  })
  loo_cohort_df <- dplyr::bind_rows(loo_cohort)
  export_csv(loo_cohort_df, "leaveoneout_by_cohort")
}

message("✅ Script finished. All results are saved in ", results_dir, "/")
