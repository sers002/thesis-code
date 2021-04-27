library(data.table)
library(Matrix.utils)

load(file = "Output/sojo_output/all_betas.RData")
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

  if (M != length(ncp)) {
    # message("Omitted are null")
    ((M * pi0 * pval_threshold) / (sum(curr_power) + M * pi0 * pval_threshold )) - q
  } else {
    # message("Omitted have dist")
    ((M * pi0 * pval_threshold) / sum(curr_power)) - q
  }
}

analytic_results_fdr <- results_merge[, c("snp", "pid", "alpha", "qtl_tag", "maf_1000g", "sigma_sq", "maf", "cis_flag")]

analytic_results_fdr_n <- data.table()
thresholds.dt <- data.table()

message("Analytic power calculations:")

sample_sizes <- c(50L, 250L, 500L, 750L, 1000L, 1500L, 2000L, 3500L, 5000L)
# sample_sizes <- c(10000, 50000, 100000, 250000, 500000, 750000, 1000000)

n_filt <- 382146
n_all  <- 11102563

analytic_results_fdr[, q2 := (2 * maf_1000g * (1 - maf_1000g) * (alpha^2))/sigma_sq]

set.seed(12345)

other_values_alt <- rbinom(n = n_all - nrow(analytic_results_fdr), size = 1, prob = 0.2)
other_values_alt[other_values_alt == 1] <- sample(
  analytic_results_fdr$q2,
  size = sum(other_values_alt == 1),
  replace = TRUE
)

other_values_alt2 <- rbinom(n = n_all - nrow(analytic_results_fdr), size = 1, prob = 0.1)
other_values_alt2[other_values_alt2 == 1] <- sample(
  analytic_results_fdr$q2,
  size = sum(other_values_alt2 == 1),
  replace = TRUE
)

other_values_alt3 <- rbinom(n = n_all - nrow(analytic_results_fdr), size = 1, prob = 0.05)
other_values_alt3[other_values_alt3 == 1] <- sample(
  analytic_results_fdr$q2,
  size = sum(other_values_alt3 == 1),
  replace = TRUE
)

other_values_alt4 <- rbinom(n = n_all - nrow(analytic_results_fdr), size = 1, prob = 0.01)
other_values_alt4[other_values_alt4 == 1] <- sample(
  analytic_results_fdr$q2,
  size = sum(other_values_alt4 == 1),
  replace = TRUE
)

other_values_same <- sample(
  analytic_results_fdr$q2,
  size = n_all - nrow(analytic_results_fdr),
  replace = TRUE
)

rbnm <- function(n, shape1, shape2, pi) {
  values <- rbinom(n = n, size = 1, prob = pi)
  values[values != 0] <- rbeta(
    n = sum(values != 0),
    shape1 = shape1,
    shape2 = shape2
  )

  values
}

for (N in sample_sizes) {
  message("  N = ", N)

  analytic_results_fdr[, NCP := (2 * N * maf_1000g * (1 - maf_1000g) * (alpha^2))/sigma_sq]

  message("    Filtered")
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


# All null ----------------------------------------------------------------

  message("    Null")
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


# Alternate distribution --------------------------------------------------

  # other_values <- rbnm(n = n_all - nrow(analytic_results_fdr), 0.1, 40, 0.4)
  # other_values <- other_values * N

  message("    Alt 1")
  threshold_optim <- uniroot(
    fdr_p_threshold,
    interval  = c(.Machine$double.eps, 1),
    tol       = .Machine$double.eps,
    ncp       = c(analytic_results_fdr$NCP, N * other_values_alt),
    M         = n_all,
    extendInt = "yes"
  )

  pval_threshold = threshold_optim$root

  critical_val = qchisq(1 - pval_threshold, 1)

  analytic_results_fdr[, power_all_alt_dist := 1 - pchisq(critical_val, 1, NCP) ]


  # Alternate distribution 2 --------------------------------------------------

  # other_values <- rbnm(n = n_all - nrow(analytic_results_fdr), 0.1, 40, 0.4)
  # other_values <- other_values * N

  message("    Alt 2")
  threshold_optim <- uniroot(
    fdr_p_threshold,
    interval  = c(.Machine$double.eps, 1),
    tol       = .Machine$double.eps,
    ncp       = c(analytic_results_fdr$NCP, N * other_values_alt2),
    M         = n_all,
    extendInt = "yes"
  )

  pval_threshold = threshold_optim$root

  critical_val = qchisq(1 - pval_threshold, 1)

  analytic_results_fdr[, power_all_alt2_dist := 1 - pchisq(critical_val, 1, NCP) ]


  # Alternate distribution 3 --------------------------------------------------

  # other_values <- rbnm(n = n_all - nrow(analytic_results_fdr), 0.1, 40, 0.4)
  # other_values <- other_values * N

  message("    Alt 3")
  threshold_optim <- uniroot(
    fdr_p_threshold,
    interval  = c(.Machine$double.eps, 1),
    tol       = .Machine$double.eps,
    ncp       = c(analytic_results_fdr$NCP, N * other_values_alt3),
    M         = n_all,
    extendInt = "yes"
  )

  pval_threshold = threshold_optim$root

  critical_val = qchisq(1 - pval_threshold, 1)

  analytic_results_fdr[, power_all_alt3_dist := 1 - pchisq(critical_val, 1, NCP) ]


  # Alternate distribution 4 --------------------------------------------------

  # other_values <- rbnm(n = n_all - nrow(analytic_results_fdr), 0.1, 40, 0.4)
  # other_values <- other_values * N

  message("    Alt 4")
  threshold_optim <- uniroot(
    fdr_p_threshold,
    interval  = c(.Machine$double.eps, 1),
    tol       = .Machine$double.eps,
    ncp       = c(analytic_results_fdr$NCP, N * other_values_alt4),
    M         = n_all,
    extendInt = "yes"
  )

  pval_threshold = threshold_optim$root

  critical_val = qchisq(1 - pval_threshold, 1)

  analytic_results_fdr[, power_all_alt4_dist := 1 - pchisq(critical_val, 1, NCP) ]

# Same distribution -------------------------------------------------------

  message("    Same")
  threshold_optim <- uniroot(
    fdr_p_threshold,
    interval  = c(.Machine$double.eps, 1),
    tol       = .Machine$double.eps,
    ncp       = c(analytic_results_fdr$NCP, N * other_values_same),
    M         = n_all,
    extendInt = "yes"
  )

  pval_threshold = threshold_optim$root

  critical_val = qchisq(1 - pval_threshold, 1)

  analytic_results_fdr[, power_all_same_dist := 1 - pchisq(critical_val, 1, NCP) ]



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

# analytic_results_fdr_n <- analytic_results_fdr_n[cis_flag == TRUE]

analytic_power_fdr <- analytic_results_fdr_n[
  ,
  list(
    power_filtered_fdr = mean(power_filt),
    power_all_fdr = mean(power_all),
    power_all_alt_fdr = mean(power_all_alt_dist),
    power_all_alt2_fdr = mean(power_all_alt2_dist),
    power_all_alt3_fdr = mean(power_all_alt3_dist),
    power_all_alt4_fdr = mean(power_all_alt4_dist),
    power_all_same_fdr = mean(power_all_same_dist)
  ),
  by = c("sample_size", "qtl_tag", "es_bin", "maf_bin", "cis_flag")
]

analytic_power_fdr <- analytic_power_fdr[order(sample_size, es_bin, maf_bin, qtl_tag)]
analytic_power_fdr <- analytic_power_fdr[qtl_tag == TRUE]

library(dplyr)
library(ggplot2)

analytic_power_fdr %>%
  melt(measure.vars = c(
    "power_filtered_fdr",
    "power_all_fdr",
    "power_all_alt_fdr",
    "power_all_alt2_fdr",
    "power_all_alt3_fdr",
    "power_all_alt4_fdr",
    "power_all_same_fdr"
  )) %>%
  subset(qtl_tag == TRUE & cis_flag == TRUE & sample_size <= 5000) %>%
  ggplot(aes(sample_size, value, colour = variable, linetype = variable)) +
  geom_point() +
  geom_line() +
  facet_grid(es_bin~maf_bin)

maf_labeller <- function(maf_bin) {
  sprintf("MAF: %s%%", maf_bin)
}

analytic_power_fdr[, maf_bin_labelled := sprintf("'%s'", maf_labeller(maf_bin))]
analytic_power_fdr[, es_bin_labelled := sprintf("alpha*': %s'", es_bin)]

analytic_power_fdr[qtl_tag == TRUE & cis_flag == TRUE] %>%
subset(es_bin != "< 0.025") %>%
  melt(
    measure.vars = c(
      "power_filtered_fdr",
      "power_all_fdr",
      "power_all_alt_fdr",
      "power_all_alt2_fdr",
      "power_all_alt3_fdr",
      "power_all_alt4_fdr",
      "power_all_same_fdr"
    ),
    variable.name = "Adjustment",
    value.name = "power"
  ) %>%
  ggplot(aes(x = sample_size, y = power, colour = Adjustment, lty = Adjustment)) +
  geom_point() +
  geom_line() +
  facet_grid(es_bin_labelled~maf_bin_labelled, labeller = label_parsed) +
  scale_colour_discrete(
    name = "Filtering",
    labels = c(
      "power_filtered_fdr"        = "Hi-C Filtered cis Variants",
      "power_all_fdr"  = "All cis Variants (Null Distribution)",
      "power_all_alt_fdr"  = "All cis Variants (20% Distribution)",
      "power_all_alt2_fdr"  = "All cis Variants (10% Distribution)",
      "power_all_alt3_fdr"  = "All cis Variants (5% Distribution)",
      "power_all_alt4_fdr"  = "All cis Variants (1% Distribution)",
      "power_all_same_fdr"  = "All cis Variants (Same Distribution)"
    ),
    guide = guide_legend(override.aes = list(size = 2, linetype = NA), nrow = 3)
  ) +
  scale_linetype_discrete(
    name = "Filtering",
    labels = c(
      "power_filtered_fdr"        = "Hi-C Filtered cis Variants",
      "power_all_fdr"  = "All cis Variants (Null Distribution)",
      "power_all_alt_fdr"  = "All cis Variants (20% Distribution)",
      "power_all_alt2_fdr"  = "All cis Variants (10% Distribution)",
      "power_all_alt3_fdr"  = "All cis Variants (5% Distribution)",
      "power_all_alt4_fdr"  = "All cis Variants (1% Distribution)",
      "power_all_same_fdr"  = "All cis Variants (Same Distribution)"
    )
  ) +
  scale_shape_discrete(
    name = "Filtering",
    labels = c(
      "power_filtered_fdr"        = "Hi-C Filtered cis Variants",
      "power_all_fdr"  = "All cis Variants (Null Distribution)",
      "power_all_alt_fdr"  = "All cis Variants (20% Distribution)",
      "power_all_alt2_fdr"  = "All cis Variants (10% Distribution)",
      "power_all_alt3_fdr"  = "All cis Variants (5% Distribution)",
      "power_all_alt4_fdr"  = "All cis Variants (1% Distribution)",
      "power_all_same_fdr"  = "All cis Variants (Same Distribution)"
    )
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power to Detect eSNPs by Hi-C Filtering",
    subtitle = "B-H FDR < 0.05",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1))

ggsave("Plots/power_comparisons/power_filtered_alt_dist.pdf", width = 7, height = 6)


