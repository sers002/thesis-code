library(data.table)
library(Matrix.utils)

load(file = here::here("Output/sojo_output/all_betas.RData"))
# combined <- combined[cis_flag == TRUE]

results_merge <- combined[, c("pid", "snp", "alpha", "maf", "qtl_tag", "b", "b_se", "beta_shrink", "sid_pos", "cis_flag", "maf_1000g", "sigma_sq")]


# Optimising function -----------------------------------------------------

## Ascertain whether given p-value threshold gives an FDR of q, given an NCP vector;
## if the function returns 0, then desired FDR has been reached
##   pval_threshold = p-value threshold to check against
##   M = number of comparisons
##   q = desired FDR level
##   ncp = vector of NCP values
fdr_p_threshold <- function(pval_threshold, M, q = 0.05, ncp, pi0 = 1) {
  ## Chisquare critical value corresponding to p-value threshold
  chisq_threshold = qchisq(1 - pval_threshold, 1)

  ## Power for each SNP using its NCP
  curr_power <- 1-pchisq(chisq_threshold, 1, ncp = ncp)

  ((M * pi0 * pval_threshold) / sum(curr_power)) - q
  ((M * pi0 * pval_threshold) / (sum(curr_power) + M * pi0 * pval_threshold )) - q
}

analytic_results_fdr <- results_merge[, c("snp", "pid", "alpha", "qtl_tag", "maf_1000g", "sigma_sq", "maf", "cis_flag")]

analytic_results_fdr_n <- data.table()
thresholds.dt <- data.table()

message("Analytic power calculations:")

sample_sizes <- c(50L, 250L, 500L, 750L, 1000L, 1500L, 2000L, 3000L, 5000L)
# sample_sizes <- c(10000, 50000, 100000, 250000, 500000, 750000, 1000000)

n_filt <- 382146
n_all  <- 11102563

timings.dt <- data.table()

for (N in sample_sizes) {
  message("  N = ", N)

  start.time <- proc.time()["elapsed"]

  analytic_results_fdr[, NCP := (2 * N * maf_1000g * (1 - maf_1000g) * (alpha^2))/sigma_sq]

  ## Find p-value threshold that results in a pFDR of q
  threshold_optim <- uniroot(
    fdr_p_threshold,
    interval  = c(.Machine$double.eps, 1),
    tol       = .Machine$double.eps,
    ncp       = analytic_results_fdr$NCP,
    M         = n_filt,
    extendInt = "yes"
  )

  pval_threshold = threshold_optim$root

  critical_val = qchisq(1 - pval_threshold, 1)

  analytic_results_fdr[, power_filt := 1 - pchisq(critical_val, 1, NCP) ]

  end.time <- proc.time()["elapsed"]

  timings.dt <- rbind(
    timings.dt,
    data.table(sample_size = N, time = end.time - start.time)
  )


  threshold_optim <- uniroot(
    fdr_p_threshold,
    interval  = c(.Machine$double.eps, 1),
    tol       = .Machine$double.eps,
    ncp       = analytic_results_fdr$NCP,
    M         = n_all,
    extendInt = "yes"
  )

  pval_threshold = threshold_optim$root

  critical_val = qchisq(1 - pval_threshold, 1)

  analytic_results_fdr[, power_all := 1 - pchisq(critical_val, 1, NCP) ]

  analytic_results_fdr[, sample_size := N]

  analytic_results_fdr_n <- rbind(analytic_results_fdr_n, analytic_results_fdr)

  thresholds.dt <- rbind(
    thresholds.dt,
    data.table(N = N, p_threshold = pval_threshold)
  )
}

es_bin_cutoffs <- c(-Inf, 0.025, 0.057, 0.13, 0.265, Inf)
analytic_results_fdr_n[, es_bin := cut(
  abs(alpha),
  es_bin_cutoffs,
  right = FALSE,
  labels = c("< 0.025", "0.025-0.057", "0.057-0.13", "0.13-0.265", "0.265+")
)]

analytic_results_fdr_n[, maf_bin := "10-20"]
analytic_results_fdr_n[maf > 0.20, maf_bin := "20-30"]
analytic_results_fdr_n[maf > 0.30, maf_bin := "30-40"]
analytic_results_fdr_n[maf > 0.40, maf_bin := "40-50"]

analytic_power_fdr <- analytic_results_fdr_n[
  ,
  list(
    power_filtered_fdr = mean(power_filt),
    power_all_fdr = mean(power_all)
  ),
  by = c("sample_size", "qtl_tag", "es_bin", "maf_bin", "cis_flag")
]

analytic_power_fdr <- analytic_power_fdr[order(sample_size, es_bin, maf_bin, qtl_tag)]
analytic_power_fdr <- analytic_power_fdr[qtl_tag == TRUE]
