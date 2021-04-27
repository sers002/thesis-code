library(data.table)
library(Matrix.utils)

load(file = "Output/sojo_output/all_betas_yri.RData")
combined[is.na(qtl_tag), qtl_tag := FALSE]

results_merge <- combined[, c(
  "pid",
  "snp",
  "alpha",
  "maf",
  "qtl_tag",
  "cis_flag",
  "b",
  "b_se",
  "beta_shrink",
  "sid_pos",
  "maf_1000g",
  "sigma_sq",
  "var_y"
)]

# Looping over genes ------------------------------------------------------

## Data table to store results of each simulation (realisations will be columns)
sa_results_big <- data.table()

ld.df <- fread(here::here("Output/linkage_disequilibrium/ld_gtex_snps_yri.ld.gz"))

## Get the genes with large numbers of SNPs over with first to save memory
gene_list <- results_merge[, .(n_snps = .N), by = "pid"][order(-n_snps)]$pid

sample_sizes <- c(50L, 250L, 500L, 750L, 1000L, 1500L, 2000L, 3000L, 4000L, 5000L)
# sample_sizes <- c(10000, 50000, 100000, 250000, 500000, 750000, 1000000)

for (g.i in seq_along(gene_list)) {
  curr_gene <- gene_list[g.i]

  results_merge2 <- results_merge[pid == curr_gene]

  results_merge <- results_merge[pid != curr_gene]

  message(g.i, ": ", curr_gene, " (", sum(results_merge2$qtl_tag), " / ", nrow(results_merge2), ")")

  ld_sub_mat <- dMcast(
    ld.df[SNP_A %in% results_merge2$snp & SNP_B %in% results_merge2$snp],
    SNP_A ~ SNP_B,
    value.var = "R"
  )
  colnames(ld_sub_mat) <- gsub("SNP_B", "", colnames(ld_sub_mat))

  results_merge2 <- results_merge2[snp %in% colnames(ld_sub_mat)]

  ld_sub_mat <- ld_sub_mat[results_merge2$snp, results_merge2$snp]

  gene_sigma <- results_merge2$sigma_sq

  for (N in sample_sizes) {
    begin_time_n <- proc.time()
    message("  N = ", N)

    ## The one SNP situation is slightly different
    if (nrow(results_merge2) == 1) {
      sigma_mat <- as.matrix(gene_sigma/(2 * N * results_merge2$maf_1000g * (1-results_merge2$maf_1000g)))
    } else {
      ## Calculated as sigma_mat[i, j] = sqrt(var_alpha_i) * sqrt(var_alpha_j) * r_ij
      se_alphas <- sqrt(gene_sigma/(2 * N * results_merge2$maf_1000g * (1-results_merge2$maf_1000g)))
      sigma_mat <- sweep(ld_sub_mat, MARGIN = 1, STATS = se_alphas, FUN = "*")
      sigma_mat <- sweep(sigma_mat, MARGIN = 2, STATS = se_alphas, FUN = "*")
    }

    ## Simulate 100 realisations of the semi-analytic betas
    sim_alphas <- mvtnorm::rmvnorm(
      n = 100,
      mean = results_merge2$alpha,
      sigma = sigma_mat,
      checkSymmetry = FALSE
    )

    ## Turn it into rows = SNPs, columns = realisations
    sim_alphas <- t(sim_alphas)

    sim.dt <- data.table(
      pid = results_merge2$pid,
      snp = results_merge2$snp,
      sample_size = N,
      sim_alphas
    )

    sa_results_big <- rbind(
      sa_results_big,
      sim.dt
    )
  }

  rm(ld_sub_mat, sim_alphas, sim.dt)
  gc()

  if (g.i %% 100 == 0) {
    saveRDS(sa_results_big, file = sprintf("sa_afr/bigger_sa_results_mvn_in_progress_%s.rds", g.i))
    sa_results_big <- NULL
    gc()

    sa_results_big <- data.table()
  }
}

saveRDS(sa_results_big, file = "Output/semi_analytic/sa_results_mvn_yri.rds")

