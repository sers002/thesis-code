library(data.table)
library(Matrix.utils)

load(file = "Output/sojo_output/all_betas.RData")
combined <- combined[cis_flag == TRUE]

results_merge <- combined[, c("pid", "snp", "alpha", "maf", "qtl_tag", "b", "b_se", "beta_shrink", "sid_pos", "maf_1000g", "sigma_sq")]


# One sample size test ----------------------------------------------------

# N <- 404
# 
# results_merge[, alpha_hat_var2 := sigma_sq/(2 * N * maf_1000g * (1 - maf_1000g))]
# 
# results_merge[, NCP := (N*2*maf_1000g*(1-maf_1000g)*(alpha^2))/sigma_sq]
# # sig_level = 0.05/nrow(results_merge)
# sig_level <- fdr_p_threshold
# threshold = qchisq(1-sig_level, 1)
# results_merge[, power := 1 - pchisq(threshold, 1, NCP) ]
# mean(results_merge[qtl_tag == TRUE]$power)
# 
# results_merge_an <- results_merge[qtl_tag == TRUE]
# 
# es_bin_cutoffs <- c(-Inf, 0.025, 0.057, 0.13, 0.265, Inf)
# results_merge_an[, es_bin := cut(
#   abs(alpha),
#   es_bin_cutoffs,
#   right = FALSE,
#   labels = c("< 0.025", "0.025-0.057", "0.057-0.13", "0.13-0.265", "0.265+")
# )]
# 
# results_merge_an[, maf_bin := "10-20"]
# results_merge_an[maf > 0.20, maf_bin := "20-30"]
# results_merge_an[maf > 0.30, maf_bin := "30-40"]
# results_merge_an[maf > 0.40, maf_bin := "40-50"]


# Multiple sample sizes ---------------------------------------------------

analytic_results_bon <- results_merge[
  qtl_tag == TRUE, 
  c("snp", "pid", "alpha", "qtl_tag", "maf_1000g", "sigma_sq", "maf")
]

## Calculate significance levels based on number of comparisons
sig_level_filt = 0.05/301548
sig_level_all  = 0.05/11102563

## Thresholds for given adjusted p-values required for significance
threshold_filt = qchisq(1 - sig_level_filt, 1)
threshold_all  = qchisq(1 - sig_level_all, 1)

analytic_results_bon_n <- data.table()

for (N in c(50, 250, 500, 750, 1000)) {
  message("N = ", N)
  
  analytic_results_bon[, NCP := (N * 2 * maf_1000g * (1-maf_1000g) * (alpha^2)) / sigma_sq]
  
  analytic_results_bon[, power_filt := 1 - pchisq(threshold_filt, 1, NCP)]
  analytic_results_bon[, power_all  := 1 - pchisq(threshold_all, 1, NCP)]
  
  analytic_results_bon[, sample_size := N]
  
  analytic_results_bon_n <- rbind(analytic_results_bon_n, analytic_results_bon)
}

es_bin_cutoffs <- c(-Inf, 0.025, 0.057, 0.13, 0.265, Inf)
analytic_results_bon_n[, es_bin := cut(
  abs(alpha),
  es_bin_cutoffs,
  right = FALSE,
  labels = c("< 0.025", "0.025-0.057", "0.057-0.13", "0.13-0.265", "0.265+")
)]

analytic_results_bon_n[, maf_bin := "10-20"]
analytic_results_bon_n[maf > 0.20, maf_bin := "20-30"]
analytic_results_bon_n[maf > 0.30, maf_bin := "30-40"]
analytic_results_bon_n[maf > 0.40, maf_bin := "40-50"]

analytic_power_bon <- analytic_results_bon_n[
  , 
  list(power_filtered = mean(power_filt)), 
  by = c("sample_size", "es_bin", "maf_bin")
]
