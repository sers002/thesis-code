library(data.table)
library(Matrix.utils)

load(file = "Output/sojo_output/all_betas.RData")

results_merge <- combined[, c("pid", "snp", "alpha", "maf", "qtl_tag", "b", "b_se", "beta_shrink", "sid_pos", "cis_flag", "maf_1000g", "sigma_sq")]

# Analysing results -------------------------------------------------------

## MVN up to 5k
sa_results_big <- readRDS("Output/semi_analytic/sa_results_mvn.rds")
sa_type <- "mvn"

## Merge simulated coefficients onto SOJO coefficients/metadata dataset
sa_results_big <- merge(
  results_merge[, c("pid", "snp", "alpha", "qtl_tag", "maf", "maf_1000g", "sigma_sq", "cis_flag")],
  sa_results_big,
  by = c("pid", "snp")
)

sa_results_big[, alpha_hat_var := sigma_sq/(2 * sample_size * maf_1000g * (1 - maf_1000g))]
sa_results_big <- sa_results_big[, 1:60]
gc()

## Function to calculate power over each realisation of the semi-analytic simulation
sa_power_calc <- function(adjusted_matrix, method_name) {
  ## Combine effect size/MAF bin and QTL metadata with discovery (adjusted p < 0.05),
  ## then calculate the percent detected for each group within the bins for a
  ## given realisation
  power_dt <- data.table(
    sa_curr_n_power[, c("es_bin", "maf_bin", "qtl_tag")],
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

ld.df <- fread(here::here("Output/linkage_disequilibrium/ld_gtex_snps.ld.gz"))

sample_sizes <- unique(sa_results_big$sample_size)

message("\nAnalysing semi-analytic results:")

n_filt <- 382146
n_all  <- 11102563

n_permutations <- 100

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

  rm(z_mat, sa_results_big_n)
  gc()

  rownames(p_mat) <- paste(results_merge$pid, results_merge$snp)

  p.dt <- data.table(reshape2::melt(p_mat, value.name = "p"))
  p.dt[, pid := gsub("(.*) (rs.*)", "\\1", Var1)]
  p.dt[, snp := gsub("(.*) (rs.*)", "\\2", Var1)]

  setnames(p.dt, "Var2", "sa_rep_i")

  min_p.dt <- p.dt[, list(min_p = min(p)), by = c("pid", "sa_rep_i")]

  n_sa_realisations <- ncol(p_mat)

  rm(p_mat, p.dt)
  sa_results_big <- sa_results_big[sample_size != curr_n]
  gc()

  gene_list <- unique(min_p.dt$pid)
  perm_results.dt <- data.table()

  pb <- txtProgressBar(max = length(gene_list), style = 3)

  for (g.i in seq_along(gene_list)) {
    curr_gene <- gene_list[g.i]

    message("  ", g.i, " - ", curr_gene)

    results_merge2 <- results_merge[pid == curr_gene]

    ld_sub_mat <- dMcast(
      ld.df[SNP_A %in% results_merge2$snp & SNP_B %in% results_merge2$snp],
      SNP_A ~ SNP_B,
      value.var = "R"
    )
    colnames(ld_sub_mat) <- gsub("SNP_B", "", colnames(ld_sub_mat))

    ld_sub_mat <- as.matrix(ld_sub_mat[results_merge2$snp, results_merge2$snp])

    perm.dt <- mvtnorm::rmvnorm(
      n = n_permutations * n_sa_realisations,
      mean = rep(0, nrow(results_merge2)),
      sigma = ld_sub_mat
    )

    rm(ld_sub_mat)
    gc()

    perm.dt <- data.table(reshape2::melt(perm.dt, value.name = "z_value"))
    perm.dt[, perm_i := rep(1:n_permutations, times = n_sa_realisations * nrow(results_merge2))]
    perm.dt[, sa_rep_i := rep(rep(1:n_sa_realisations, each = n_permutations), times = nrow(results_merge2))]

    perm.dt <- perm.dt[, list(max_z = max(abs(z_value))), by = c("perm_i", "sa_rep_i")]
    perm.dt[, perm_p :=  pnorm(-abs(max_z), 0, 1, lower = T) + pnorm(abs(max_z), 0, 1, lower = F)]

    perm.dt <- merge(
      min_p.dt[pid == curr_gene],
      perm.dt,
      by = "sa_rep_i"
    )[, list( perm_p = sum(perm_p <= min_p) / .N ), by = c("pid", "sa_rep_i")]

    perm_results.dt <- rbind(perm_results.dt, perm.dt)

    rm(perm.dt)
    gc()
  }

  close(pb)

  perm_results.dt[, sample_size := curr_n]

  sa_power_all <- rbind(sa_power_all, perm_results.dt)

  rm(perm_results.dt)

  saveRDS(sa_power_all, file = sprintf("perm_genes_sa_in_progress_%s.rds", curr_n))

  gc()
}

library(dplyr)
library(ggplot2)

sa_power_all[, list(n = sum(perm_p < 0.05) / .N), by = c("sample_size", "sa_rep_i")][,list( power = mean(n)), by = "sample_size"] %>%
ggplot(aes(sample_size, power)) + geom_point()

sa_power_all[, list(n = sum(perm_p < 0.05) / .N), by = c("sample_size", "sa_rep_i")][,list( power = mean(n)), by = "sample_size"] %>%
 fwrite(file = "Output/power_egenes/sa_perm.csv")
