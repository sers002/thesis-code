library(data.table)
library(Matrix.utils)

load(file = "Output/sojo_output/all_betas.RData")
combined <- combined[cis_flag == TRUE]

results_merge <- combined[, c("pid", "snp", "alpha", "maf", "qtl_tag", "b", "b_se", "beta_shrink", "sid_pos", "cis_flag", "maf_1000g", "sigma_sq")]


# Analysing results -------------------------------------------------------

## MVN up to 5k
sa_results_big <- readRDS("Output/semi_analytic/bigger_sa_results_mvn.rds")
sa_type <- "mvn"

## Merge simulated coefficients onto SOJO coefficients/metadata dataset
sa_results_big <- merge(
  results_merge[, c("pid", "snp", "alpha", "qtl_tag", "maf", "maf_1000g", "sigma_sq", "cis_flag")],
  sa_results_big,
  by = c("pid", "snp")
)

sa_results_big[, alpha_hat_var := sigma_sq/(2 * sample_size * maf_1000g * (1 - maf_1000g))]

## Function to calculate power over each realisation of the semi-analytic simulation
sa_power_calc <- function(
  coef_results,
  adjusted_matrix,
  method_name,
  by = c("qtl_tag", "es_bin", "cis_flag")
) {
  coef_results <- coef_results[, by, with = FALSE]
  ## Combine effect size/MAF bin and QTL metadata with discovery (adjusted p < 0.05),
  ## then calculate the percent detected for each group within the bins for a
  ## given realisation
  power_dt <- data.table(
    coef_results,
    adjusted_matrix < 0.05
  )[,
    lapply(.SD, function(x) sum(x)/.N),
    by = by,
    .SDcols = paste0("V", 1:ncol(adjusted_matrix))
  ]

  ## Calculate the average power over all realisations. Calculate the Monte Carlo
  ## interval based on the quantiles of power
  power_dt <- melt(
    power_dt,
    measure.vars = paste0("V", 1:ncol(adjusted_matrix))
  )[
    ,
    list(
      power = mean(value),
      power_lwr = quantile(value, probs = 0.025),
      power_upr = quantile(value, probs = 0.975)
    ),
    by = by
  ]

  power_dt[, method := method_name]

  setkeyv(power_dt, cols = by)

  power_dt
}

## Datasets to store results of loop in
sa_power_all   <- data.table()
sa_fdr_cutoffs <- data.table()

sample_sizes <- unique(sa_results_big$sample_size)

message("\nAnalysing semi-analytic results:")

n_filt <- 382146
n_all  <- 11102563 ### cis
# n_all <- 301177637 ### trans

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

  message(nrow(f_mat_filt), ", ", ncol(f_mat_filt))

  f_mat_all <- apply(p_mat, 2, p.adjust, method = "BH", n = n_all)

  i <- 0

  message("alt1")
  f_mat_all_alt <- apply(p_mat, 2, function(p) {
    i <<- i + 1
    set.seed(i)
    message(i)

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

  message("alt2")
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

  message("alt3")
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

  message("alt4")
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

  message("same")
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

  power_dt_filt_bh <- sa_power_calc(sa_curr_n_power, f_mat_filt, "filtered_fdr")
  power_dt_all_bh  <- sa_power_calc(sa_curr_n_power,f_mat_all, "all_fdr")

  power_dt_all_bh_alt1 <- sa_power_calc(sa_curr_n_power,f_mat_all_alt, "all_fdr_alt")
  power_dt_all_bh_alt2 <- sa_power_calc(sa_curr_n_power,f_mat_all_alt2, "all_fdr_alt2")
  power_dt_all_bh_alt3 <- sa_power_calc(sa_curr_n_power,f_mat_all_alt3, "all_fdr_alt3")
  power_dt_all_bh_alt4 <- sa_power_calc(sa_curr_n_power,f_mat_all_alt4, "all_fdr_alt4")
  power_dt_all_bh_same <- sa_power_calc(sa_curr_n_power,f_mat_all_same, "all_fdr_same")

  ## Join all of the power datasets together
  sa_curr_n_power_merge <- rbind(
    power_dt_filt_bh,
    power_dt_all_bh,
    power_dt_all_bh_alt1,
    power_dt_all_bh_alt2,
    power_dt_all_bh_alt3,
    power_dt_all_bh_alt4,
    power_dt_all_bh_same
  )

  sa_curr_n_power_merge[, sample_size := curr_n]


  sa_power_all <- rbind(sa_power_all, sa_curr_n_power_merge)

  rm(p_mat)
  rm(f_mat_all, f_mat_filt, b_mat_filt, b_mat_all, by_mat_filt, by_mat_all, q_mat_filt)
  rm(sa_curr_n_power, sa_curr_n_power_merge)
  rm(power_dt_filt_bh, power_dt_all_bh, power_dt_filt_by, power_dt_all_by, power_dt_filt_bon, power_dt_all_bon)

  gc()
}

sa_power_all_averaged <- sa_power_all[order(sample_size, cis_flag)]

# fwrite(sa_power_all_averaged, file = "Output/power/sa_mvn_filtering_cis.csv")
# sa_power_all_averaged <- fread("Output/power/sa_mvn_filtering_cis.csv")

sa_power_all_averaged[, es_bin_labelled := sprintf("beta*': %s'", es_bin)]

sa_power_all_averaged[, method2 := factor(method, levels = c("filtered_fdr", "all_fdr", "all_fdr_same", "all_fdr_alt", "all_fdr_alt2", "all_fdr_alt3", "all_fdr_alt4"))]

sa_power_all_averaged[qtl_tag == TRUE & cis_flag == TRUE] %>%
  ggplot(aes(x = sample_size, y = power, colour = method2)) +
  geom_point() +
  geom_line() +
  facet_wrap(~es_bin_labelled, labeller = label_parsed) +
  scale_colour_discrete(
    name = "Filtering",
    labels = c(
      "filtered_fdr"        = "Hi-C Filtered cis Variants",
      "all_fdr"  = "All cis Variants (Null Distribution)",
      "all_fdr_alt"  = "All cis Variants (20% Distribution)",
      "all_fdr_alt2"  = "All cis Variants (10% Distribution)",
      "all_fdr_alt3"  = "All cis Variants (5% Distribution)",
      "all_fdr_alt4"  = "All cis Variants (1% Distribution)",
      "all_fdr_same"  = "All cis Variants (Same Distribution)"
    ),
    guide = guide_legend(override.aes = list(size = 2, linetype = NA), nrow = 3)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power to Detect cis-eQTLs by Hi-C Filtering",
    subtitle = "BH FDR < 0.05",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1))

ggsave("Plots/power_comparisons/power_filtering.pdf", width = 7.8, height = 5)
