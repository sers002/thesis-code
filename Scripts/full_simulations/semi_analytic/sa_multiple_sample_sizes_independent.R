library(data.table)
library(Matrix.utils)

load(file = "Output/sojo_output/all_betas.RData")
# combined <- combined[cis_flag == TRUE]

results_merge <- combined[, c("pid", "snp", "alpha", "maf", "qtl_tag", "b", "b_se", "beta_shrink", "sid_pos", "cis_flag", "maf_1000g", "sigma_sq")]


# Simulating alphas -------------------------------------------------------

sa_results_big <- data.table()

message("Simulating alphas (independent method):")

# sample_sizes <- c(50L, 250L, 500L, 750L, 1000L, 2000L, 3000L, 4000L, 5000L)
sample_sizes <- c(10000, 50000, 100000, 250000, 500000, 750000, 1000000)

for (N in sample_sizes) {
    message("  N = ", N)

  results_merge[, alpha_hat_var := sigma_sq/(2 * N * maf_1000g * (1-maf_1000g))]

  sim_alphas <- rnorm(
    n    = nrow(results_merge) * 100,
    mean = results_merge$alpha,
    sd   = sqrt(results_merge$alpha_hat_var)
  )

  sim.mat <- matrix(sim_alphas, ncol = 100)

  sim.dt <- data.table(
    pid         = results_merge$pid,
    snp         = results_merge$snp,
    sample_size = N,
    sim.mat
  )

  sa_results_big <- rbind(
    sa_results_big,
    sim.dt
  )
}

rm(sample_sizes)
