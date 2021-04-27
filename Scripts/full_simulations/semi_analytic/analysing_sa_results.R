library(data.table)
library(Matrix.utils)

load(file = here::here("Output/sojo_output/all_betas.RData"))
# combined <- combined[cis_flag == TRUE]

results_merge <- combined[, c("pid", "snp", "alpha", "maf", "qtl_tag", "b", "b_se", "beta_shrink", "sid_pos", "cis_flag", "maf_1000g", "sigma_sq")]


# Analysing results -------------------------------------------------------

## MVN up to 5k
sa_results_big <- readRDS(here::here("Output/semi_analytic/sa_results_mvn.rds"))
sa_type <- "mvn"

## Merge simulated coefficients onto SOJO coefficients/metadata dataset
sa_results_big <- merge(
  results_merge[, c("pid", "snp", "alpha", "qtl_tag", "maf", "maf_1000g", "sigma_sq", "cis_flag")],
  sa_results_big,
  by = c("pid", "snp")
)

sa_results_big[, sigma_sq := 1 - (alpha^2 * 2 * maf_1000g * (1 - maf_1000g))]
sa_results_big[, alpha_hat_var := sigma_sq/(2 * sample_size * maf_1000g * (1 - maf_1000g))]

## Function to calculate power over each realisation of the semi-analytic simulation
sa_power_calc <- function(adjusted_matrix, method_name) {
  ## Combine effect size/MAF bin and QTL metadata with discovery (adjusted p < 0.05),
  ## then calculate the percent detected for each group within the bins for a
  ## given realisation
  power_dt <- data.table(
    sa_curr_n_power[, c("es_bin", "maf_bin", "qtl_tag", "cis_flag")],
    adjusted_matrix < 0.05
  )[,
    lapply(.SD, function(x) sum(x)/.N),
    by = c("qtl_tag", "es_bin", "maf_bin"),
    .SDcols = paste0("V", 1:100)
  ]

  ## Calculate the average power over all realisations. Calculate the Monte Carlo
  ## interval based on the quantiles of power
  power_dt <- melt(power_dt, measure.vars = paste0("V", 1:100))[
    ,
    list(
      power = mean(value),
      power_lwr = quantile(value, probs = 0.025),
      power_upr = quantile(value, probs = 0.975)
    ),
    by = c("qtl_tag", "es_bin", "maf_bin")
  ][order(qtl_tag, es_bin, maf_bin)]

  setnames(
    power_dt,
    c("power", "power_lwr", "power_upr"),
    sprintf("power_%s%s", method_name, c("", "_lwr", "_upr"))
  )

  setkey(power_dt, qtl_tag, es_bin, maf_bin)

  power_dt
}

## Datasets to store results of loop in
sa_power_all   <- data.table()
sa_fdr_cutoffs <- data.table()

sample_sizes <- unique(sa_results_big$sample_size)

message("\nAnalysing semi-analytic results:")

n_filt <- 382146
n_all  <- 11102563

for (curr_n in sample_sizes) {
  message("  N = ", curr_n)

  sa_results_big_n <- sa_results_big[sample_size == curr_n]

  ## For each SNP, divide by its SE - realisations are the V1, V2, ... V100 columns
  z_mat <- as.matrix(sa_results_big_n[, grep("V[0-9]+", names(sa_results_big_n)), with = FALSE])
  z_mat <- sweep(z_mat, MARGIN = 1, sqrt(sa_results_big_n$alpha_hat_var), FUN = "/")

  p_mat <- matrix(
    pnorm(-abs(z_mat), 0, 1, lower = T) + pnorm(abs(z_mat), 0, 1, lower = F),
    ncol = ncol(z_mat)
  )

  rm(z_mat)
  gc()

  ## For each realisation (a column), adjust the p-values
  f_mat_filt  <- apply(p_mat, 2, p.adjust, method = "BH", n = n_filt)
  by_mat_filt <- apply(p_mat, 2, p.adjust, method = "BY", n = n_filt)
  b_mat_filt  <- apply(p_mat, 2, p.adjust, method = "bonferroni", n = n_filt)
  q_mat_filt  <- apply(p_mat, 2, function(p) {
    p <- c(p, runif(n_filt - length(p)))
    qvalue::qvalue(p)$qvalues[1:length(p)]
  })

  f_mat_all  <- apply(p_mat, 2, p.adjust, method = "BH", n = n_all)
  by_mat_all <- apply(p_mat, 2, p.adjust, method = "BY", n = n_all)
  b_mat_all  <- apply(p_mat, 2, p.adjust, method = "bonferroni", n = n_all)

  sa_curr_n_power <- sa_results_big_n[, 1:8]

  es_bin_cutoffs <- c(-Inf, 0.025, 0.057, 0.13, 0.265, Inf)
  sa_curr_n_power[, es_bin := cut(
    abs(alpha),
    es_bin_cutoffs,
    right = FALSE,
    labels = c("< 0.025", "0.025-0.057", "0.057-0.13", "0.13-0.265", "0.265+")
  )]

  sa_curr_n_power[, maf_bin := "10-20"]
  sa_curr_n_power[maf > 0.20, maf_bin := "20-30"]
  sa_curr_n_power[maf > 0.30, maf_bin := "30-40"]
  sa_curr_n_power[maf > 0.40, maf_bin := "40-50"]

  power_dt_filt_bh <- sa_power_calc(f_mat_filt, "filtered_fdr")
  power_dt_all_bh  <- sa_power_calc(f_mat_all, "all_fdr")

  power_dt_filt_by <- sa_power_calc(by_mat_filt, "filtered_by")
  power_dt_all_by  <- sa_power_calc(by_mat_all, "all_by")

  power_dt_filt_bon <- sa_power_calc(b_mat_filt, "filtered_bon")
  power_dt_all_bon  <- sa_power_calc(b_mat_all, "all_bon")

  power_dt_all_q   <- sa_power_calc(q_mat_filt, "filtered_qval")

  ## Join all of the power datasets together
  sa_curr_n_power_merge <- power_dt_filt_bh[
    power_dt_all_bh
  ][
    power_dt_filt_by
  ][
    power_dt_all_by
  ][
    power_dt_filt_bon
  ][
    power_dt_all_bon
  ][
    power_dt_all_q
  ]

  sa_curr_n_power_merge[, sample_size := curr_n]

  sa_power_all <- rbind(sa_power_all, sa_curr_n_power_merge)

  rm(p_mat)
  rm(f_mat_all, f_mat_filt, b_mat_filt, b_mat_all, by_mat_filt, by_mat_all, q_mat_filt)
  rm(sa_curr_n_power, sa_curr_n_power_merge)
  rm(power_dt_filt_bh, power_dt_all_bh, power_dt_filt_by, power_dt_all_by, power_dt_filt_bon, power_dt_all_bon)

  gc()
}

sa_power_all_averaged <- sa_power_all[order(sample_size, es_bin, maf_bin)]

fwrite(sa_power_all_averaged, file = "Output/power/sa_mvn.csv")

