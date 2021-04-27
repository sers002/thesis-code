library(data.table)
library(Matrix.utils)

ancestry <- "_yri"

load(file = sprintf("Output/sojo_output/all_betas%s.RData", ancestry))
combined <- combined[cis_flag == TRUE]
ids_in_all <- readLines("ancestry_ids_in_all.txt")
combined[, id := paste(pid, snp)]

results_merge <- combined[id %in% ids_in_all, c("pid", "snp", "alpha", "maf", "qtl_tag", "b", "b_se", "beta_shrink", "sid_pos", "cis_flag", "maf_1000g", "sigma_sq")]

# Analysing results -------------------------------------------------------

## MVN up to 5k
sa_results_big <- readRDS(sprintf("Output/semi_analytic/sa_results_mvn%s.rds", ancestry))
sa_type <- "mvn"

sa_results_big <- merge(
  results_merge[, c("pid", "snp", "alpha", "qtl_tag", "maf", "maf_1000g", "sigma_sq", "cis_flag", "beta_shrink")],
  sa_results_big,
  by = c("pid", "snp")
)

sa_results_big[, alpha_hat_var := sigma_sq/(2 * sample_size * maf_1000g * (1 - maf_1000g))]

## Datasets to store results of loop in
sa_power_all <- data.table()
sa_fdr_cutoffs <- data.table()
sa_pi0 <- data.table()

sample_sizes <- unique(sa_results_big$sample_size)

message("\nAnalysing semi-analytic results:")

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
  f_mat_filt <- apply(p_mat, 2, p.adjust, method = "BH", n = 382146)
  by_mat_filt <- apply(p_mat, 2, p.adjust, method = "BY", n = 382146)

  message(nrow(f_mat_filt), ", ", ncol(f_mat_filt))

  b_mat_filt <- apply(p_mat, 2, p.adjust, method = "bonferroni", n = 382146)
  f_mat_all <- apply(p_mat, 2, p.adjust, method = "BH", n = 11102563)

  # sa_curr_n_power <- sa_results_big_n[qtl_tag == TRUE, 1:7]
  sa_curr_n_power_genes <- sa_results_big_n[, 1:10][, list(
    egene = any(beta_shrink != 0)
  ),
  by = "pid"
  ]

  message(" egenes: ", sum(sa_curr_n_power_genes$egene))

  sa_curr_n_power_genes$n_got_fdr_filt <- rowSums(
    apply(f_mat_filt, MARGIN = 2, function(f) {
      fdr.df <- data.frame(
        fdr = f,
        qtl_tag = sa_results_big_n$qtl_tag,
        pid = sa_results_big_n$pid
      )
      by(
        fdr.df,
        factor(fdr.df$pid),
        function(f.df) {
          any(f.df$fdr < 0.05 & f.df$qtl_tag == TRUE)
        })
      })
  )

  sa_curr_n_power_genes$n_got_fdr_all <- rowSums(
    apply(f_mat_all, MARGIN = 2, function(f) {
      fdr.df <- data.frame(
        fdr = f,
        qtl_tag = sa_results_big_n$qtl_tag,
        pid = sa_results_big_n$pid
      )
      by(
        fdr.df,
        factor(fdr.df$pid),
        function(f.df) {
          any(f.df$fdr < 0.05 & f.df$qtl_tag == TRUE)
        })
    })
  )

  sa_curr_n_power_genes$n_got_by_filt <- rowSums(
    apply(by_mat_filt, MARGIN = 2, function(f) {
      fdr.df <- data.frame(
        fdr = f,
        qtl_tag = sa_results_big_n$qtl_tag,
        pid = sa_results_big_n$pid
      )
      by(
        fdr.df,
        factor(fdr.df$pid),
        function(f.df) {
          any(f.df$fdr < 0.05 & f.df$qtl_tag == TRUE)
        })
    })
  )


  sa_curr_n_power_genes$n_got_bon_filt <- rowSums(
    apply(b_mat_filt, MARGIN = 2, function(b) {
      tapply(b, sa_results_big_n$pid, function(b2) any(b2 < 0.05))
    })
  )

  N_realisations <- ncol(f_mat_filt)

  sa_curr_n_power_genes[, pct_correct_fdr_filt := n_got_fdr_filt/N_realisations]
  sa_curr_n_power_genes[, pct_correct_by_filt  := n_got_by_filt/N_realisations]
  sa_curr_n_power_genes[, pct_correct_fdr_all  := n_got_fdr_all/N_realisations]
  sa_curr_n_power_genes[, pct_correct_bon_filt := n_got_bon_filt/N_realisations]

  sa_curr_n_power_genes[, sample_size := curr_n]

  sa_power_all <- rbind(sa_power_all, sa_curr_n_power_genes)

  rm(p_mat, f_mat_all, f_mat_filt, b_mat_filt, b_mat_all, sa_curr_n_power)
  gc()
}

sa_power_all_averaged <- sa_power_all[, .(
  power_filtered_fdr = mean(pct_correct_fdr_filt),
  power_filtered_by = mean(pct_correct_by_filt),
  power_all_fdr      = mean(pct_correct_fdr_all),
  power_filtered_bon = mean(pct_correct_bon_filt)
),
  by = c("sample_size", "egene")
]

sa_power_all_averaged <- sa_power_all_averaged[order(sample_size, egene)]

fwrite(sa_power_all_averaged, file = sprintf("Output/power_egenes/sa_mvn%s.csv", ifelse(ancestry != "", ancestry, "_eur")))

