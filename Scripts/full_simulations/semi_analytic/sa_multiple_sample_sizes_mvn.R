library(data.table)
library(Matrix.utils)

load(file = here::here("Output/sojo_output/all_betas.RData"))

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


# Nullifying Genes --------------------------------------------------------

# codes3d_gtex_genes.dt <- fread("Data/lung_eqtl_data/significant_eqtls.txt")
# codes3d_gtex_genes.dt <- codes3d_gtex_genes.dt[gene_chr == "chr19"]
# uniqueN(codes3d_gtex_genes.dt$gencode_id)
#
# codes3d_genes <- unique(results_merge$pid)[unique(results_merge$pid) %in% unique(codes3d_gtex_genes.dt$gencode_id)]
# set.seed(12345)
# null_genes <- sample(unique(results_merge$pid), size = ceiling(0.75*length(unique(results_merge$pid))))
# null_genes <- sample(
#   codes3d_genes,
#   size = ceiling(0.5*length(codes3d_genes)),
#   replace = FALSE
# )

# set.seed(123456)
# null_genes <- null_genes[sample(1:length(null_genes), size = ceiling(0.5 * length(null_genes)))]

# results_merge[pid %in% null_genes, beta_shrink := 0]
# results_merge[pid %in% null_genes, alpha := 0]
# results_merge[pid %in% null_genes, qtl_tag := FALSE]

# results_merge[!(pid %in% codes3d_genes), beta_shrink := 0]
# results_merge[!(pid %in% codes3d_genes), alpha := 0]
# results_merge[!(pid %in% codes3d_genes), qtl_tag := FALSE]

# results_merge[, sigma_sq := var_y - (alpha^2 * 2 * maf_1000g * (1 - maf_1000g))]
# results_merge <- results_merge[pid %in% null_genes | !(pid %in% codes3d_genes)]
# results_merge <- results_merge[pid %in% null_genes]

# rm(codes3d_gtex_genes.dt, combined)

# Looping over genes ------------------------------------------------------

## Data table to store results of each simulation (realisations will be columns)
sa_results_big <- data.table()

ld.df <- fread(here::here("Output/linkage_disequilibrium/ld_gtex_snps.ld.gz"))

## Get the genes with large numbers of SNPs over with first to save memory
gene_list <- results_merge[, .(n_snps = .N), by = "pid"][order(-n_snps)]$pid

sample_sizes <- c(50L, 250L, 500L, 750L, 1000L, 2000L, 3000L, 5000L)
# sample_sizes <- c(10000, 50000, 100000, 250000, 500000, 750000, 1000000)
# sample_sizes <- c(2000L, 5000L, 10000L, 15000L, 20000L, 30000L, 40000L, 50000L)

for (g.i in seq_along(gene_list)[1700:1710]) {
  curr_gene <- gene_list[g.i]

  results_merge2 <- results_merge[pid == curr_gene]
  results_merge  <- results_merge[pid != curr_gene]

  message(g.i, ": ", curr_gene, " (", sum(results_merge2$qtl_tag), " / ", nrow(results_merge2), ")")

  ld_sub_mat <- dMcast(
    ld.df[SNP_A %in% results_merge2$snp & SNP_B %in% results_merge2$snp],
    SNP_A ~ SNP_B,
    value.var = "R"
  )
  colnames(ld_sub_mat) <- gsub("SNP_B", "", colnames(ld_sub_mat))

  ld_sub_mat <- as.matrix(ld_sub_mat[results_merge2$snp, results_merge2$snp])

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
      sigma = sigma_mat
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
}

saveRDS(sa_results_big, file = "Output/semi_analytic/sa_results_mvn.rds")

