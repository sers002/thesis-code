library(data.table)
library(Matrix.utils)

load(file = "Output/sojo_output/all_betas.RData")
combined <- combined[cis_flag == TRUE]

results_merge <- combined[, c("pid", "snp", "alpha", "maf", "qtl_tag", "b", "b_se", "beta_shrink", "sid_pos", "cis_flag", "maf_1000g", "sigma_sq")]

justin_sig.dt <- fread("Data/lung_eqtl_data/significant_eqtls.txt")
justin_sig_genes <- unique(justin_sig.dt[gene_chr == "chr19"]$gencode_id)


# Analysing results -------------------------------------------------------

## Read in previously-calculated semi-analytic results

## Alphas treated as independent (runs fast)
# source("Scripts/full_simulations/2_semi_analytic/sa_multiple_sample_sizes_independent.R")
# sa_type <- "independent"

## MVN up to 5k
# sa_results_big <- readRDS("Output/semi_analytic/bigger_sa_results_all_indiv_sigmas.rds")
sa_results_big <- readRDS("Output/semi_analytic/bigger_sa_results_mvn.rds")
sa_type <- "mvn"

## Betas wrapped-up up to 10k?
# sa_results_big <- readRDS("Output/semi_analytic/bigger_sa_results_betas_wrapped.rds")
# sa_type <- "wrapup"

## For bigger_sa_results_betas_wrapped - 204 onwards is non-zeroed out betas
# sa_results_big <- sa_results_big[, c(1:3, 104:203)]

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

n_filt <- 382146
n_all  <- 11102563

message("\nAnalysing semi-analytic results:")

for (curr_n in sample_sizes) {
  message("  N = ", curr_n)
  
  sa_results_big_n <- sa_results_big[sample_size == curr_n]
  sa_results_big <- sa_results_big[sample_size != curr_n]
  
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
  # by_mat_filt <- apply(p_mat, 2, p.adjust, method = "BY", n = 382146)
  
  message(nrow(f_mat_filt), ", ", ncol(f_mat_filt))
  
  # b_mat_filt <- apply(p_mat, 2, p.adjust, method = "bonferroni", n = 382146)
  f_mat_all <- apply(p_mat, 2, p.adjust, method = "BH", n = 11102563)
  
  i <- 0
  
  f_mat_all_alt <- apply(p_mat, 2, function(p) {
    i <<- i + 1
    message(i)
    set.seed(i)
    
    other_values_alt <- sample(
      c(-1, 1), 
      size = n_all - length(p), 
      replace = TRUE, 
      prob = c(0.8, 0.2)
    )

    other_values_alt[other_values_alt == 1] <- sample(
      p,
      size = sum(other_values_alt == 1),
      replace = TRUE
    )
    
    other_values_alt[other_values_alt == -1] <- runif(n = sum(other_values_alt == -1))
    
    p.adjust(c(p, other_values_alt), method = "BH", n = n_all)[1:length(p)]
  })
  
  f_mat_all_alt2 <- apply(p_mat, 2, function(p) {
    i <<- i + 1
    set.seed(i)
    
    other_values_alt <- sample(
      c(-1, 1), 
      size = n_all - length(p), 
      replace = TRUE, 
      prob = c(0.9, 0.1)
    )
    
    other_values_alt[other_values_alt == 1] <- sample(
      p,
      size = sum(other_values_alt == 1),
      replace = TRUE
    )
    
    other_values_alt[other_values_alt == -1] <- runif(n = sum(other_values_alt == -1))
    
    p.adjust(c(p, other_values_alt), method = "BH", n = n_all)[1:length(p)]
  })
  
  f_mat_all_alt3 <- apply(p_mat, 2, function(p) {
    i <<- i + 1
    set.seed(i)
    
    other_values_alt <- sample(
      c(-1, 1), 
      size = n_all - length(p), 
      replace = TRUE, 
      prob = c(0.95, 0.05)
    )
    
    other_values_alt[other_values_alt == 1] <- sample(
      p,
      size = sum(other_values_alt == 1),
      replace = TRUE
    )
    
    other_values_alt[other_values_alt == -1] <- runif(n = sum(other_values_alt == -1))
    
    p.adjust(c(p, other_values_alt), method = "BH", n = n_all)[1:length(p)]
  })
  
  f_mat_all_alt4 <- apply(p_mat, 2, function(p) {
    i <<- i + 1
    set.seed(i)
    
    other_values_alt <- sample(
      c(-1, 1), 
      size = n_all - length(p), 
      replace = TRUE, 
      prob = c(0.99, 0.01)
    )
    
    other_values_alt[other_values_alt == 1] <- sample(
      p,
      size = sum(other_values_alt == 1),
      replace = TRUE
    )
    
    other_values_alt[other_values_alt == -1] <- runif(n = sum(other_values_alt == -1))
    
    p.adjust(c(p, other_values_alt), method = "BH", n = n_all)[1:length(p)]
  })
  
  f_mat_all_same <- apply(p_mat, 2, function(p) {
    i <<- i + 1
    set.seed(i)
    
    other_values_alt <- sample(
      p, 
      size = n_all - length(p), 
      replace = TRUE
    )
    
    p.adjust(c(p, other_values_alt), method = "BH", n = n_all)[1:length(p)]
  })
  
  # b_mat_all <- apply(p_mat, 2, p.adjust, method = "bonferroni", n = 11102563)
  
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
  
  got_mat_fun <- function(f) {
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
  }
  
  sa_curr_n_power_genes$n_got_fdr_all <- rowSums(
    apply(f_mat_all, MARGIN = 2, got_mat_fun)
  )
  
  sa_curr_n_power_genes$n_got_fdr_all_alt <- rowSums(
    apply(f_mat_all_alt, MARGIN = 2, got_mat_fun)
  )
  
  sa_curr_n_power_genes$n_got_fdr_all_alt2 <- rowSums(
    apply(f_mat_all_alt2, MARGIN = 2, got_mat_fun)
  )
  
  sa_curr_n_power_genes$n_got_fdr_all_alt3 <- rowSums(
    apply(f_mat_all_alt3, MARGIN = 2, got_mat_fun)
  )
  
  sa_curr_n_power_genes$n_got_fdr_all_alt4 <- rowSums(
    apply(f_mat_all_alt4, MARGIN = 2, got_mat_fun)
  )
  
  sa_curr_n_power_genes$n_got_fdr_all_same <- rowSums(
    apply(f_mat_all_same, MARGIN = 2, got_mat_fun)
  )
  
  N_realisations <- ncol(f_mat_filt)
  
  sa_curr_n_power_genes[, pct_correct_fdr_filt := n_got_fdr_filt/N_realisations]
  sa_curr_n_power_genes[, pct_correct_fdr_all  := n_got_fdr_all/N_realisations]
  sa_curr_n_power_genes[, pct_correct_fdr_all_alt  := n_got_fdr_all_alt/N_realisations]
  sa_curr_n_power_genes[, pct_correct_fdr_all_alt2  := n_got_fdr_all_alt2/N_realisations]
  sa_curr_n_power_genes[, pct_correct_fdr_all_alt3  := n_got_fdr_all_alt3/N_realisations]
  sa_curr_n_power_genes[, pct_correct_fdr_all_alt4  := n_got_fdr_all_alt4/N_realisations]
  sa_curr_n_power_genes[, pct_correct_fdr_all_same  := n_got_fdr_all_same/N_realisations]
  # sa_curr_n_power[, pct_correct_bon_all  := n_got_bon_all/N_realisations]
  
  sa_curr_n_power_genes[, sample_size := curr_n]
  
  sa_power_all <- rbind(sa_power_all, sa_curr_n_power_genes)
  
  rm(p_mat, f_mat_all, f_mat_filt, b_mat_filt, b_mat_all, sa_curr_n_power, f_mat_all_alt, f_mat_all_alt2, f_mat_all_alt3, f_mat_all_alt4, f_mat_all_same)
  gc()
}
# sa_power_all[, egene_both := egene == TRUE & egene_justin == TRUE]

# sa_power_all[, .(mean(pct_correct_fdr_filt)), by = "sample_size"]
sa_power_all_averaged <- sa_power_all[, .(
  power_filtered_fdr = mean(pct_correct_fdr_filt),
  power_all_fdr      = mean(pct_correct_fdr_all),
  power_all_fdr_alt      = mean(pct_correct_fdr_all_alt),
  power_all_fdr_alt2      = mean(pct_correct_fdr_all_alt2),
  power_all_fdr_alt3      = mean(pct_correct_fdr_all_alt3),
  power_all_fdr_alt4      = mean(pct_correct_fdr_all_alt4),
  power_all_fdr_same      = mean(pct_correct_fdr_all_same)
  # power_all_bon      = mean(pct_correct_bon_all)
), 
by = c("sample_size", "egene")
]

sa_power_all_averaged <- sa_power_all_averaged[order(sample_size, egene)]

# sa_power_all_averaged <- sa_power_all_averaged[egene == TRUE]
# sa_power_all_averaged$qtl_tag <- NULL

# if (sa_type == "independent") {
#   fwrite(sa_power_all_averaged, file = "Output/power_egenes/sa_independent.csv")
# } else if (sa_type == "mvn") {
  fwrite(sa_power_all_averaged, file = "Output/power_egenes/sa_mvn_filtered.csv")
# } else if (sa_type == "wrapup") {
#   fwrite(sa_power_all_averaged, file = "Output/power_egenes/sa_wrapup.csv")
# }
