library(dplyr)
library(ggplot2)
library(data.table)

maf_labeller <- function(maf_bin) {
  sprintf("MAF: %s%%", maf_bin)
}

me_power <- fread(here::here("Output/power/matrixEQTL.csv"))
me_power <- me_power[qtl_tag == FALSE]

mvn_power    <- fread(here::here("Output/power/sa_mvn.csv"))

## Calculate the analytic power on-the-fly
sample_sizes <- c(50L, 250L, 500L, 750L, 1000L, 1500L, 2000L)
sample_sizes <- unique(mvn_power$sample_size)
source(here::here("Scripts/full_simulations/analytic/analytic_power_fdr.R"))

all_power <- rbindlist(list(
  data.table(me_power, method = "me"),
  data.table(mvn_power[cis_flag == TRUE], method = "sa_mvn"),
  data.table(analytic_power_fdr[cis_flag == TRUE], method = "an")
), fill = TRUE)

all_power[, maf_bin_labelled := sprintf("'%s'", maf_labeller(maf_bin))]
all_power[, es_bin_labelled := sprintf("beta*': %s'", es_bin)]

all_power %>%
  # subset((method %in% c("an"))) %>%
  # subset((es_bin %in% c("< 0.025", "0.025-0.057"))) %>%
  ggplot(aes(x = sample_size, y = power_filtered_fdr, colour = method, lty = qtl_tag)) +
  geom_point() +
  geom_line() +
  # geom_errorbar(
  #   aes(ymin = power_filtered_fdr_lwr, ymax = power_filtered_fdr_upr)
  #   # position = position_dodge(width = 0.5)
  # ) +
  facet_grid(es_bin_labelled~maf_bin_labelled, labeller = label_parsed) +
  scale_colour_manual(
    name = "Estimation Method",
    labels = c(
      "sa_mvn"        = "Semi-analytic",
      "me"  = "Empirical",
      "an"        = "Analytic"
    ),
    values = c(
      "sa_mvn" = "#5FA1F7",
      "an" = "#83A552",
      "me" = "#9B1F1A"
    ),
    guide = guide_legend(override.aes = list(size = 2), nrow = 1)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Comparison of Power to Detect eQTLs by Estimation Method ",
    subtitle = "BH FDR < 0.05",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1))

ggsave("Plots/power_comparisons_newer/all_causal_snp/power_all_methods.pdf", width = 7, height = 6.5)
# ggsave("Plots/power_comparisons_newer/all_causal_snp/power_all_intervals.pdf", width = 7, height = 6.5)
# ggsave("Plots/power_comparisons_newer/all_causal_snp/power_extended.pdf", width = 7, height = 5)


all_power_trans <- rbindlist(list(
  data.table(mvn_power[qtl_tag == TRUE],        method = "sa_mvn"),
  data.table(analytic_power_fdr[qtl_tag == TRUE], method = "an")
), fill = TRUE)

all_power_trans[, maf_bin_labelled := sprintf("'%s'", maf_labeller(maf_bin))]
all_power_trans[, es_bin_labelled := sprintf("beta*': %s'", es_bin)]

all_power_trans %>%
  subset(method %in% c("an")) %>%
  ggplot(aes(x = sample_size, y = power_filtered_fdr, colour = cis_flag)) +
  geom_point() +
  geom_line() +
  facet_grid(es_bin_labelled~maf_bin_labelled, labeller = label_parsed) +
  scale_colour_manual(
    name = "Method",
    labels = c(
      "FALSE"        = "Trans",
      "TRUE"  = "Cis",
      "an"        = "Analytic",
      "me"        = "Empirical",
      "me_zeroed" = "Empirical (beta=0 removed)",
      "an_pm"     = "Analytic PowerMap",
      "sa_wrapup" = "Semi-analytic (Wrap-up)"
    ),
    values = c(
      "FALSE" = "#00ba38",
      "TRUE" = "#619cff",
      "me" = "#f8766d",
      "an_pm" = "magenta",
      "sa_kerby" = "green",
      "sa_wrapup" = "violetred",
      "me_zeroed" = "chocolate4",
      "sa_kerby_non1" = "blue"
    ),
    guide = guide_legend(override.aes = list(size = 2), nrow = 2)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power to Detect trans-eQTLS",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1)) +
  xlim(c(0, 2000))




## Comparing filtering/adjustment method
all_power %>%
  filter(method %in% c("an")) %>%
  subset(es_bin != "< 0.025") %>%
  melt(
    measure.vars = c("power_filtered_fdr", "power_all_fdr"),
    variable.name = "Adjustment",
    value.name = "power"
  ) %>%
  ggplot(aes(x = sample_size, y = power, colour = Adjustment, lty = Adjustment, shape = Adjustment)) +
  geom_point() +
  geom_line() +
  facet_grid(es_bin_labelled~maf_bin_labelled, labeller = label_parsed) +
  scale_colour_manual(
    name = "Filtering",
    labels = c(
      "power_filtered_fdr"        = "cis Variants in Hi-C Contact",
      "power_all_fdr"  = "All cis Variants"
    ),
    values = c(
      "power_filtered_fdr" = "#A7473A",
      "power_all_fdr" = "#B09B37"
    ),
    guide = guide_legend(override.aes = list(size = 2, linetype = NA), nrow = 1)
  ) +
  scale_linetype_discrete(
    name = "Filtering",
    labels = c(
      "power_filtered_fdr"        = "cis Variants in Hi-C Contact",
      "power_all_fdr"  = "All cis Variants"
    )
  ) +
  scale_shape_discrete(
    name = "Filtering",
    labels = c(
      "power_filtered_fdr"        = "cis Variants in Hi-C Contact",
      "power_all_fdr"  = "All cis Variants"
    )
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power to Detect eQTLs by Hi-C Filtering",
    subtitle = "BH FDR < 0.05",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1))

ggsave("Plots/power_comparisons_newer/all_causal_snp/power_filtered_5k.pdf", width = 7, height = 6)
# ggsave("Plots/power_comparisons_newer/all_causal_snp/power_filtered_extended.pdf", width = 7, height = 6)

## Comparing filtering/adjustment method
all_power %>%
  filter(method %in% c("sa_mvn")) %>%
  subset(es_bin != "< 0.025") %>%
  melt(
    measure.vars = c("power_filtered_fdr", "power_filtered_by", "power_filtered_bon", "power_filtered_qval"),
    variable.name = "Adjustment",
    value.name = "power"
  ) %>%
  ggplot(aes(x = sample_size, y = power, colour = Adjustment)) +
  geom_point() +
  geom_line() +
  facet_grid(es_bin_labelled~maf_bin_labelled, labeller = label_parsed) +
  scale_linetype(
    labels = c(
      "power_filtered_fdr" = "FDR (BH)",
      "power_filtered_by"      = "FDR (BY)",
      "power_filtered_bon"      = "Bonferroni",
      "power_filtered_qval"        = "FDR (Storey)"
    ),
    guide = guide_legend(nrow =1)
  ) +
  scale_shape(
    labels = c(
      "power_filtered_fdr" = "FDR (BH)",
      "power_filtered_by"      = "FDR (BY)",
      "power_filtered_bon"      = "Bonferroni",
      "power_filtered_qval"        = "FDR (Storey)"
    ),
    guide = guide_legend(nrow = 1)
  ) +
  scale_colour_manual(
    name = "Multiple Testing\nAdjustment",
    labels = c(
      "power_filtered_fdr" = "FDR (BH)",
      "power_filtered_by"  = "FDR (BY)",
      "power_filtered_bon"        = "Bonferroni",
      "power_filtered_qval"        = "FDR (Storey)"
    ),
    values = c(
      "power_filtered_fdr" = "#6C803A",
      "power_filtered_by" = "#AB7C47",
      "power_filtered_bon" = "#CCAE42",
      "power_filtered_qval" = "#D73202"
    ),
    guide = guide_legend(override.aes = list(size = 2), nrow = 2)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Comparison of Power to Detect eQTLs by Multiple Testing Adjustment Method",
    # subtitle = "All causal variants",
    x = "Sample Size",
    y = "Power",
    shape = "Multiple Testing\nAdjustment",
    linetype = "Multiple Testing\nAdjustment"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1))

ggsave("Plots/power_comparisons/power_adjustment_comparison.pdf", width = 7, height = 5.5)



# Comparing semi-analytic variance ----------------------------------------------------

sa_power   <- fread("Output/power/sa_mvn.csv")
sa_power2  <- fread("Output/power/sa_mvn_sigma1.csv")
sa_power3  <- fread("Output/power/sa_mvn_sigmasq1varexplained.csv")

sa_comp.df <- rbindlist(list(
  data.table(sa_power[qtl_tag == TRUE & cis_flag == TRUE],  method = "sa_mvn"),
  data.table(sa_power2[qtl_tag == TRUE],  method = "sa_mvn_sigma1"),
  data.table(sa_power3[qtl_tag == TRUE],  method = "sa_mvn_sigmavarexp")
), fill = TRUE) %>%
  subset(es_bin %in% c("0.13-0.265", "0.265+")) %>%
  subset(sample_size %in% c(50, 250, 500, 750, 1000, 1500, 2000))

sa_comp.df[, maf_bin_labelled := sprintf("'%s'", maf_labeller(maf_bin))]
sa_comp.df[, es_bin_labelled := sprintf("beta*': %s'", es_bin)]

sa_comp.df[, method := factor(method, levels = c("sa_mvn_sigma1", "sa_mvn_sigmavarexp", "sa_mvn"))]

ggplot(sa_comp.df, aes(x = sample_size, y = power_filtered_fdr, colour = method)) +
  geom_line() +
  geom_point(alpha = 0.5) +
  facet_grid(es_bin_labelled~maf_bin_labelled, labeller = label_parsed) +
  scale_colour_manual(
    name = "Error Variance",
    labels = c(
      "sa_mvn_sigma1"  = bquote(sigma[e]^2 == 1),
      "sa_mvn_sigmavarexp"  = bquote(sigma[e]^2 == "1 - " * beta^2 * 2 * p * (1-p)),
      "sa_mvn" = bquote(sigma[e]^2 == hat("Var(Y)")*" - " * beta^2 * 2 * p * (1-p))
    ),
    values = c(
      "sa_mvn_sigma1" = "green",
      "sa_mvn_sigmavarexp" = "purple",
      "sa_mvn" = "blue"
    ),
    guide = guide_legend(override.aes = list(size = 1), ncol = 1)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Comparison of Error Variance Terms",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1)) +
  xlim(c(0, 2000))

ggsave(
  "Plots/power_comparisons/residual_error_power_comparison.pdf",
  width = 7,
  height = 4
)


