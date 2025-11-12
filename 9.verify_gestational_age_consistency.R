> #!/usr/bin/env Rscript
  > 
  > suppressPackageStartupMessages({
    +     library(dplyr)
    +     library(readr)
    +     library(stringr)
    +     library(TwoSampleMR)
    +     library(ggplot2)
    + })
> 
  > cat("=== SYSTEMATIC CHECK: preterm birth vs gestational age (and post-term) ===\n\n")
=== SYSTEMATIC CHECK: preterm birth vs gestational age (and post-term) ===
  
  > 
  > # Inputs
  > mr_results_path <- "results/mr_results_30_outcomes_fdr.csv"
  > harm_path       <- "results/harmonised_rahmioglu_bpo.csv"
  > out_csv         <- "results/check_preterm_vs_ga_report.csv"
  > 
    > stopifnot(file.exists(mr_results_path), file.exists(harm_path))
  > mr_results <- read_csv(mr_results_path, show_col_types = FALSE)
  > harm_data  <- read_csv(harm_path,       show_col_types = FALSE)                                                              
  >                                                                                                                              
    > # 1) Global MR comparison (add post-term if present)
    > temporal_keep <- c("pretb_all", "pretb_subsamp", "ga_all", "posttb_all")
    > temporal_outcomes <- mr_results %>%
      +     filter(outcome %in% temporal_keep) %>%
      +     select(outcome, outcome_clean, OR, OR_lower, OR_upper, pval, b, se) %>%
      +     arrange(outcome)
    > 
      > cat("Global MR results (temporal outcomes):\n\n")
    Global MR results (temporal outcomes):
      
      > apply(temporal_outcomes, 1, function(r){
        +     cat(sprintf("%s:\n", r[["outcome_clean"]]))
        +     cat(sprintf("  - OR = %.3f (%.3f–%.3f)\n", as.numeric(r[["OR"]]),
                          +                 as.numeric(r[["OR_lower"]]), as.numeric(r[["OR_upper"]])))
        +     cat(sprintf("  - Beta = %.4f ± %.4f\n", as.numeric(r[["b"]]), as.numeric(r[["se"]])))
        +     cat(sprintf("  - P = %.6f\n\n", as.numeric(r[["pval"]])))
        + })
    Gestational age (all):
      - OR = 1.108 (0.964–1.274)
    - Beta = 0.1026 ± 0.0713
    - P = 0.150148
    
    Post-term birth:
      - OR = 1.038 (0.882–1.220)
    - Beta = 0.0370 ± 0.0828
    - P = 0.654867
    
    Preterm birth (all):
      - OR = 0.822 (0.712–0.948)
    - Beta = -0.1965 ± 0.0732
    - P = 0.007311
    
    Preterm birth (subsample):
      - OR = 0.839 (0.717–0.982)
    - Beta = -0.1757 ± 0.0803
    - P = 0.028716
    
    NULL
    > 
      > pretb_result <- mr_results %>% filter(outcome == "pretb_all") %>% slice(1)
    > ga_result    <- mr_results %>% filter(outcome == "ga_all")    %>% slice(1)
    > posttb_result<- mr_results %>% filter(outcome == "posttb_all")%>% slice(1)
    > 
      > logic_ok <- NA
    > if (nrow(pretb_result) && nrow(ga_result)) {
      +     pretb_dir <- ifelse(pretb_result$OR > 1, "increase", "decrease")
      +     ga_dir    <- ifelse(ga_result$OR > 1, "increase", "decrease")
      +     logic_ok  <- (pretb_dir == "decrease" && ga_dir == "increase") ||
        +         (pretb_dir == "increase" && ga_dir == "decrease")
      +     cat("EXPECTED LOGIC:\n- endometriosis ↑ → preterm ↑ & GA ↓ OR endometriosis ↓ → preterm ↓ & GA ↑\n\n")
      +     cat("OBSERVED:\n")
      +     cat(sprintf("- Preterm birth (all): OR = %.3f → %s prematurity\n",
                        +                 pretb_result$OR, ifelse(pretb_result$OR > 1, "MORE", "LESS")))
      +     cat(sprintf("- Gestational age (all): OR = %.3f → %s gestational age\n\n",
                        +                 ga_result$OR, ifelse(ga_result$OR > 1, "LONGER", "SHORTER")))
      +     cat(sprintf("COHERENCE: %s\n\n", ifelse(logic_ok, "✅ CONSISTENT", "❌ INCONSISTENT")))
      + }
    EXPECTED LOGIC:
      - endometriosis ↑ → preterm ↑ & GA ↓ OR endometriosis ↓ → preterm ↓ & GA ↑
    
    OBSERVED:
      - Preterm birth (all): OR = 0.822 → LESS prematurity
    - Gestational age (all): OR = 1.108 → LONGER gestational age
    
    COHERENCE: ✅ CONSISTENT
    
    > 
      > check_sign_consistency <- function(or, b) {
        +     if (is.na(or) || is.na(b)) return(NA)
        +     (or > 1 & b > 0) || (or < 1 & b < 0)
        + }
    > pretb_sign_ok <- if (nrow(pretb_result)) check_sign_consistency(pretb_result$OR, pretb_result$b) else NA
    > ga_sign_ok    <- if (nrow(ga_result))    check_sign_consistency(ga_result$OR,    ga_result$b)    else NA
    > cat(sprintf("Sign consistency — Preterm: %s | GA: %s\n\n",
                  +             ifelse(isTRUE(pretb_sign_ok), "OK", "Check"),
                  +             ifelse(isTRUE(ga_sign_ok),    "OK", "Check")))
    Sign consistency — Preterm: OK | GA: OK
    
    > 
      > # 2) SNP-level comparison (preterm vs GA)
      > pretb_snps <- harm_data %>%
      +     filter(outcome == "pretb_all") %>%
      +     select(SNP, beta.exposure, beta.outcome) %>%
      +     rename(beta_pretb = beta.outcome)
    > 
      > ga_snps <- harm_data %>%
      +     filter(outcome == "ga_all") %>%
      +     select(SNP, beta.exposure, beta.outcome) %>%
      +     rename(beta_ga = beta.outcome)
    > 
      > combined_snps <- inner_join(pretb_snps, ga_snps, by = c("SNP","beta.exposure")) %>%
      +     mutate(
        +         opposite_dir     = sign(beta_pretb) != sign(beta_ga),
        +         pretb_concordant = sign(beta.exposure) == sign(beta_pretb),
        +         ga_concordant    = sign(beta.exposure) == sign(beta_ga)
        +     )
    > 
      > if (nrow(combined_snps)) {
        +     cat("SNP-LEVEL COMPARISON (preterm vs GA):\n")
        +     cat(sprintf("- Common SNPs: %d\n", nrow(combined_snps)))
        +     cat(sprintf("- Opposite directions (pretb vs GA): %d/%d (%.1f%%)\n",
                          +                 sum(combined_snps$opposite_dir), nrow(combined_snps),
                          +                 100 * mean(combined_snps$opposite_dir)))
        +     r_pg <- suppressWarnings(cor(combined_snps$beta_pretb, combined_snps$beta_ga, use = "complete.obs"))
        +     cat(sprintf("- Cor(beta_pretb, beta_ga) = %.3f\n\n", r_pg))
        +     
          +     snp_compare_csv <- "results/check_preterm_vs_ga_snp_table.csv"
          +     readr::write_csv(combined_snps, snp_compare_csv)
          +     cat("Saved SNP table to:", snp_compare_csv, "\n")
          +     
            +     offenders <- combined_snps %>% filter(!opposite_dir)
            +     cat(sprintf("SNPs not showing opposite directions (pretb vs GA): %d\n",
                              +                 nrow(offenders)))
            +     if (nrow(offenders)) {
              +         cat(paste(head(offenders$SNP, 10), collapse = ", "),
                            +             if (nrow(offenders) > 10) " ..." else "", "\n", sep = "")
              +         offenders_csv <- "results/check_preterm_vs_ga_offenders.csv"
              +         readr::write_csv(offenders, offenders_csv)
              +         cat("Saved offenders list to:", offenders_csv, "\n")
              +     }
            +     
              +     missing_in_ga    <- setdiff(pretb_snps$SNP, ga_snps$SNP)
              +     missing_in_pretb <- setdiff(ga_snps$SNP,    pretb_snps$SNP)
              +     cat(sprintf("Overlap: %d common | %d only in pretb | %d only in GA\n",
                                +                 length(intersect(pretb_snps$SNP, ga_snps$SNP)),
                                +                 length(missing_in_ga), length(missing_in_pretb)))
              +     
                +     p_scatter <- ggplot(combined_snps, aes(x = beta_pretb, y = beta_ga)) +
                  +         geom_point(alpha = 0.8) +
                  +         geom_abline(slope = -1, intercept = 0, linetype = "dashed") +
                  +         labs(title = "SNP effects: preterm vs gestational age",
                                 +              x = expression(beta["preterm"]), y = expression(beta["GA"])) +
                  +         theme_minimal()
                +     scatter_png <- "results/check_preterm_vs_ga_scatter.png"
                +     ggsave(scatter_png, p_scatter, width = 6, height = 5, dpi = 300, bg = "white")
                +     cat("Saved scatter to:", scatter_png, "\n\n")
                + }
    SNP-LEVEL COMPARISON (preterm vs GA):
      - Common SNPs: 34
    - Opposite directions (pretb vs GA): 26/34 (76.5%)
    - Cor(beta_pretb, beta_ga) = -0.837
    
    Saved SNP table to: results/check_preterm_vs_ga_snp_table.csv                                                                
    SNPs not showing opposite directions (pretb vs GA): 8
    rs10983311, rs12030576, rs1352889, rs2421985, rs3803042, rs56090796, rs66683298, rs7924571
    Saved offenders list to: results/check_preterm_vs_ga_offenders.csv                                                           
    Overlap: 34 common | 0 only in pretb | 0 only in GA
    Saved scatter to: results/check_preterm_vs_ga_scatter.png 
    
    > 
      > # 2b) SNP-level comparison (preterm vs post-term), if post-term exists
      > postterm_available <- nrow(posttb_result) == 1
    > posttb_snps <- NULL
    > combined_pt <- NULL
    > if (postterm_available) {
      +     posttb_snps <- harm_data %>%
        +         filter(outcome == "posttb_all") %>%
        +         select(SNP, beta.exposure, beta.outcome) %>%
        +         rename(beta_posttb = beta.outcome)
      +     
        +     combined_pt <- inner_join(pretb_snps, posttb_snps, by = c("SNP","beta.exposure")) %>%
          +         mutate(opposite_dir = sign(beta_pretb) != sign(beta_posttb))
        +     
          +     if (nrow(combined_pt)) {
            +         cat("SNP-LEVEL COMPARISON (preterm vs post-term):\n")
            +         cat(sprintf("- Common SNPs: %d\n", nrow(combined_pt)))
            +         cat(sprintf("- Opposite directions (pretb vs post-term): %d/%d (%.1f%%)\n",
                                  +                     sum(combined_pt$opposite_dir), nrow(combined_pt),
                                  +                     100 * mean(combined_pt$opposite_dir)))
            +         r_pp <- suppressWarnings(cor(combined_pt$beta_pretb, combined_pt$beta_posttb, use = "complete.obs"))
            +         cat(sprintf("- Cor(beta_pretb, beta_posttb) = %.3f\n\n", r_pp))
            +         pt_csv <- "results/check_preterm_vs_postterm_snp_table.csv"
            +         readr::write_csv(combined_pt, pt_csv)
            +         cat("Saved SNP table to:", pt_csv, "\n\n")
            +     }
        + }
    SNP-LEVEL COMPARISON (preterm vs post-term):
      - Common SNPs: 34
    - Opposite directions (pretb vs post-term): 24/34 (70.6%)
    - Cor(beta_pretb, beta_posttb) = -0.711
    
    Saved SNP table to: results/check_preterm_vs_postterm_snp_table.csv                                                          
    
    > 
      > # 3) MR diagnostics on pretb_all (if columns available)
      > has_cols <- all(c("se.exposure","se.outcome","pval.outcome") %in% names(harm_data))
    > diag_summary <- tibble()
    > 
      > .p_from_b_se <- function(b, se) { z <- abs(b / se); 2 * stats::pnorm(z, lower.tail = FALSE) }
    > .std_names   <- function(df) { names(df) <- tolower(names(df)); df }
    > 
      > if (has_cols) {
        +     dat <- harm_data %>%
          +         filter(outcome == "pretb_all") %>%
          +         mutate(across(c(beta.exposure, se.exposure, beta.outcome, se.outcome, pval.outcome), as.numeric)) %>%
          +         filter(is.finite(beta.exposure), is.finite(se.exposure),
                           +                is.finite(beta.outcome),  is.finite(se.outcome))
        +     
          +     if (nrow(dat)) {
            +         F_i <- (dat$beta.exposure^2) / (dat$se.exposure^2)
            +         F_median <- median(F_i, na.rm = TRUE); F_min <- min(F_i, na.rm = TRUE)
            +         
              +         dat_s <- steiger_filtering(dat)
              +         st_tab <- table(dat_s$steiger_dir, useNA = "ifany")
              +         st_true  <- ifelse(!is.na(st_tab["TRUE"]),  st_tab["TRUE"], 0)
              +         st_false <- ifelse(!is.na(st_tab["FALSE"]), st_tab["FALSE"], 0)
              +         
                +         res <- mr(dat_s) %>%
                  +             mutate(OR = exp(b),
                                       +                    CI_low = exp(b - 1.96*se),
                                       +                    CI_high = exp(b + 1.96*se))
                +         
                  +         het   <- mr_heterogeneity(dat_s)
                  +         pleio <- tryCatch(mr_pleiotropy_test(dat_s), error = function(e) NULL)
                  +         
                    +         loo <- mr_leaveoneout(dat_s) %>% .std_names()
                    +         if (!"p" %in% names(loo)) loo$p <- .p_from_b_se(loo$b, loo$se)
                    +         snp_col    <- if ("snp" %in% names(loo)) loo$snp else if ("id" %in% names(loo)) loo$id else rep(NA_character_, nrow(loo))
                    +         method_col <- if ("method" %in% names(loo)) loo$method else rep(NA_character_, nrow(loo))
                    +         
                      +         loo_sum <- tibble(
                        +             SNP    = snp_col,
                        +             method = method_col,
                        +             OR     = exp(loo$b),
                        +             CI_low = exp(loo$b - 1.96 * loo$se),
                        +             CI_high= exp(loo$b + 1.96 * loo$se),
                        +             p      = loo$p
                        +         )
                      +         
                        +         ivw_row <- res %>% filter(method == "Inverse variance weighted") %>%
                          +             transmute(method, OR, OR_low = CI_low, OR_high = CI_high, pval)
                        +         
                          +         diag_summary <- tibble(
                            +             outcome = "pretb_all",
                            +             F_median = round(F_median, 1),
                            +             F_min = round(F_min, 1),
                            +             Steiger_TRUE = as.integer(st_true),
                            +             Steiger_FALSE = as.integer(st_false),
                            +             IVW_OR = round(ivw_row$OR, 3),
                            +             IVW_CI = sprintf("(%.3f–%.3f)", ivw_row$OR_low, ivw_row$OR_high),
                            +             IVW_p = signif(ivw_row$pval, 3),
                            +             Egger_intercept_p = if (!is.null(pleio)) signif(pleio$pval, 3) else NA
                            +         )
                          +         
                            +         cat("MR DIAGNOSTICS (pretb_all):\n"); print(diag_summary); cat("\n")
                          +         
                            +         loo_csv <- "results/check_preterm_leaveoneout_detail.csv"
                            +         tryCatch(readr::write_csv(loo_sum, loo_csv), error = function(e) NULL)
                            +         
                              +         dat_flip <- dat_s; dat_flip$beta.outcome <- -dat_flip$beta.outcome
                              +         res_flip <- mr(dat_flip) %>% dplyr::filter(method == "Inverse variance weighted") %>%
                                +             mutate(OR = exp(b), CI_low = exp(b - 1.96*se), CI_high = exp(b + 1.96*se)) %>%
                                +             select(method, OR, CI_low, CI_high, pval)
                              +         cat("Sanity flip (pretb_all, IVW):\n"); print(res_flip); cat("\n")
                              +     }
        + }
    Analysing 'pvh6wX' on 'WdIiHx'
    MR DIAGNOSTICS (pretb_all):
      # A tibble: 1 × 9
      outcome   F_median F_min Steiger_TRUE Steiger_FALSE IVW_OR IVW_CI          IVW_p Egger_intercept_p
    <chr>        <dbl> <dbl>        <int>         <int>  <dbl> <chr>           <dbl>             <dbl>
      1 pretb_all     37.8  29.7           33             1  0.822 (0.712–0.948) 0.00731             0.719
    
    Analysing 'pvh6wX' on 'WdIiHx'
    Sanity flip (pretb_all, IVW):
      method          OR      CI_low     CI_high           pval
    1 Inverse variance weighted 1.217090144 1.054330258 1.404975726 0.007311188425
    
    > 
      > # 4) GA unit/sense checks
      > ga_rows <- harm_data %>% filter(outcome == "ga_all")
    > ga_meta_cols <- intersect(names(ga_rows), c("units.outcome","unit.outcome","units","outcome_units"))
    > ga_units <- if (length(ga_meta_cols)) unique(na.omit(unlist(ga_rows[ga_meta_cols]))) else character(0)
    > cat("GA unit metadata (if any): ", if (length(ga_units)) paste(ga_units, collapse=" | ") else "none found", "\n", sep = "")
    GA unit metadata (if any): none found
    > 
      > if (nrow(ga_rows)) {
        +     frac_pos_beta <- mean(sign(ga_rows$beta.outcome) > 0, na.rm = TRUE)
        +     cat(sprintf("GA sign check: fraction of SNPs with beta.outcome > 0 = %.1f%%\n", 100*frac_pos_beta))
        +     if (frac_pos_beta > 0.5) {
          +         cat("Interpretation hint: majority β>0 → GA tends to increase (consistent with OR>1 if OR reported).\n\n")
          +     } else {
            +         cat("Interpretation hint: majority β<=0 → GA tends to decrease; double-check coding if global OR suggests increase.\n\n")
            +     }
        + }
    GA sign check: fraction of SNPs with beta.outcome > 0 = 50.0%
    Interpretation hint: majority β<=0 → GA tends to decrease; double-check coding if global OR suggests increase.
    
    > 
      > # 5) Export concise report (add post-term where available)
      > current_scenario <- dplyr::case_when(
        +     nrow(pretb_result) && nrow(ga_result) && pretb_result$OR < 1 & ga_result$OR > 1 ~ 1,
        +     nrow(pretb_result) && nrow(ga_result) && pretb_result$OR < 1 & ga_result$OR < 1 ~ 2,
        +     nrow(pretb_result) && nrow(ga_result) && pretb_result$OR > 1 & ga_result$OR < 1 ~ 3,
        +     nrow(pretb_result) && nrow(ga_result) && pretb_result$OR > 1 & ga_result$OR > 1 ~ 4,
        +     TRUE ~ 0
        + )
    > 
      > report <- tibble(
        +     outcome_pretb_OR = if (nrow(pretb_result)) pretb_result$OR else NA_real_,
        +     outcome_ga_OR    = if (nrow(ga_result))    ga_result$OR    else NA_real_,
        +     outcome_posttb_OR= if (nrow(posttb_result)) posttb_result$OR else NA_real_,
        +     logic_coherent   = logic_ok,
        +     n_common_snps_pg = nrow(combined_snps),
        +     pct_opposite_dir_pg = if (nrow(combined_snps)) round(100*mean(combined_snps$opposite_dir),1) else NA_real_,
        +     cor_beta_pretb_beta_ga = if (nrow(combined_snps)) round(cor(combined_snps$beta_pretb, combined_snps$beta_ga, use="complete.obs"),3) else NA_real_,
        +     n_common_snps_pp = if (!is.null(combined_pt)) nrow(combined_pt) else NA_integer_,
        +     pct_opposite_dir_pp = if (!is.null(combined_pt) && nrow(combined_pt)) round(100*mean(combined_pt$opposite_dir),1) else NA_real_,
        +     scenario_code    = current_scenario,
        +     pretb_sign_ok    = isTRUE(pretb_sign_ok),
        +     ga_sign_ok       = isTRUE(ga_sign_ok),
        +     ga_units_meta    = if (length(ga_units)) paste(ga_units, collapse=" | ") else NA_character_
        + ) %>% bind_cols(diag_summary)
    > 
      > write_csv(report, out_csv)
    > cat("Saved report to:", out_csv, "\n")                                                                                     
    Saved report to: results/check_preterm_vs_ga_report.csv 
    > 
      > if (!is.na(logic_ok)) {
        +     if (logic_ok) cat("✅ Preterm ↓ & GA ↑ (or inverse) — logical directions are consistent.\n")
        +     else          cat("🚨 Inconsistency detected — recheck coding & harmonisation.\n")
        + }
    ✅ Preterm ↓ & GA ↑ (or inverse) — logical directions are consistent.