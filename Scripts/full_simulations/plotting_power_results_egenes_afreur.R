library(dplyr)
library(ggplot2)
library(data.table)

maf_labeller <- function(maf_bin) {
  sprintf("MAF: %s%%", maf_bin)
}

## Without alpha=0 removed
me_power_eur     <- fread("Output/power_egenes/sa_mvn_eur.csv")
# me_power_eur$sample_size <- me_power_eur$n

me_power_afr     <- fread("Output/power_egenes/sa_mvn_yri.csv")
# me_power_afr$sample_size <- me_power_afr$n

# me_power_eas     <- fread("Output/power_egenes/sa_mvn_eas.csv")
# me_power_eas$sample_size <- me_power_eas$n


# source("Scripts/full_simulations/3_analytic/analytic_power_fdr.R")

all_power <- rbindlist(list(
  data.table(me_power_eur[egene == TRUE],     method = "EUR"),
  data.table(me_power_afr[egene == TRUE & !(sample_size %in% c(1500, 4000))],     method = "YRI")
  # data.table(me_power_eas[egene == TRUE],     method = "EAS")
), fill = TRUE)
  
all_power %>% 
  # filter(method %in% c("sa_mvn", "sa_mvn2", "sa_mvn3")) %>% 
  ggplot(aes(x = sample_size, y = power_filtered_fdr, colour = method)) +
  geom_point() +
  geom_line() +
  # facet_grid(es_bin~maf_bin, labeller = labeller(maf_bin = maf_labeller)) +
  scale_colour_manual(
    name = "Ancestry",
    labels = c(
      "YRI" = "YRI",
      "EAS" = "East Asian",
      "EUR" = "EUR"
    ),
    values = c(
      "YRI" = "#C582B2",
      "EAS" = "#4d5f8e",
      "EUR" = "#7d9fc2"
    ),
    guide = guide_legend(override.aes = list(size = 2), nrow = 1)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power to Detect cis-eGenes by Ancestry",
    subtitle = "BH FDR < 0.05",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1))

ggsave("Plots/power_comparisons_newer/power_egenes_ancestry.pdf", width = 7, height = 4)




all_power %>% 
  filter(method %in% c("sa_mvn", "me")) %>% 
  ggplot(aes(x = sample_size, y = power_filtered_fdr, colour = method)) +
  geom_point() +
  geom_line() +
  # facet_grid(es_bin~maf_bin, labeller = labeller(maf_bin = maf_labeller)) +
  scale_colour_manual(
    name = "Method",
    labels = c(
      "sa_mvn3"        = "Only GTex-significant eGenes",
      "sa_kerby"  = "Semi-analytic (Kerby)",
      "me"        = "Empirical",
      "sa_mvn"        = "Semi-Analytic",
      "me_zeroed" = "Empirical (alpha=0 removed)",
      "an_pm"     = "Analytic PowerMap",
      "sa_perm" = "Permutation",
      "sa_mvn2" = "50% of GTEx-significant eGenes"
    ),
    values = c(
      "sa_mvn" = "#00ba38",
      "sa_mvn2" = "#619cff",
      "me" = "#f8766d",
      "sa_mvn3" = "magenta",
      "sa_kerby" = "green",
      "sa_perm" = "violetred",
      "me_zeroed" = "chocolate4",
      "sa_kerby_non1" = "blue"
    ),
    guide = guide_legend(override.aes = list(size = 2), nrow = 1)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Comparison of Power to Detect eGenes by Estimation Method",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1)) + 
  xlim(c(0, 2000))

 ggsave("Plots/power_comparisons_newer/power_egenes_empirical_sa.pdf", width = 7, height = 4)


all_power %>% 
  filter(method %in% c("sa_mvn", "sa_perm")) %>% 
  ggplot(aes(x = sample_size, y = power_filtered_fdr, colour = method)) +
  geom_point() +
  geom_line() +
  # facet_grid(es_bin~maf_bin, labeller = labeller(maf_bin = maf_labeller)) +
  scale_colour_manual(
    name = "Method",
    labels = c(
      "sa_mvn3"        = "Only GTex-significant eGenes",
      "sa_kerby"  = "Semi-analytic (Kerby)",
      "me"        = "Empirical",
      "sa_mvn"        = "Semi-Analytic",
      "me_zeroed" = "Empirical (alpha=0 removed)",
      "an_pm"     = "Analytic PowerMap",
      "sa_perm" = "Permutation (eGene-MVN)",
      "sa_mvn2" = "50% of GTEx-significant eGenes"
    ),
    values = c(
      "sa_mvn" = "#00ba38",
      "sa_mvn2" = "#619cff",
      "me" = "#f8766d",
      "sa_mvn3" = "magenta",
      "sa_kerby" = "green",
      "sa_perm" = "violetred",
      "me_zeroed" = "chocolate4",
      "sa_kerby_non1" = "blue"
    ),
    guide = guide_legend(override.aes = list(size = 2), nrow = 1)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power to Detect eGenes Using eGene-MVN",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1)) + 
  xlim(c(0, 2000))

ggsave("Plots/power_comparisons_newer/power_egenes_permutation.pdf", width = 7, height = 4)


all_power %>% 
  filter(method %in% c("sa_mvn")) %>% 
  melt(measure.vars = c("power_filtered_fdr", "power_all_fdr")) %>% 
  ggplot(aes(x = sample_size, y = value, colour = variable)) +
  geom_point() +
  geom_line() +
  # facet_grid(es_bin~maf_bin, labeller = labeller(maf_bin = maf_labeller)) +
  scale_colour_manual(
    name = "Method",
    labels = c(
      "power_filtered_fdr"        = "Filtered",
      "power_all_fdr"  = "Unfiltered",
      "an"        = "Analytic",
      "sa_mvn"        = "SOJO eGenes",
      "me_zeroed" = "Empirical (alpha=0 removed)",
      "an_pm"     = "Analytic PowerMap",
      "sa_perm" = "Permutation",
      "sa_mvn2" = "50% of GTEx-significant eGenes"
    ),
    values = c(
      "power_filtered_fdr" = "#00ba38",
      "power_all_fdr" = "#619cff",
      "me" = "#f8766d",
      "sa_mvn3" = "magenta",
      "sa_kerby" = "green",
      "sa_perm" = "violetred",
      "me_zeroed" = "chocolate4",
      "sa_kerby_non1" = "blue"
    ),
    guide = guide_legend(override.aes = list(size = 2), nrow = 1)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power to Detect eGenes by Hi-C Filtering",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1))

ggsave("Plots/power_comparisons_newer/power_egenes_filtering.pdf", width = 7, height = 4)



all_power %>% 
  filter(method %in% c("sa_mvn")) %>% 
  melt(measure.vars = c("power_filtered_fdr", "power_filtered_by", "power_filtered_bon")) %>% 
  ggplot(aes(x = sample_size, y = value, colour = variable)) +
  geom_point() +
  geom_line() +
  # facet_grid(es_bin~maf_bin, labeller = labeller(maf_bin = maf_labeller)) +
  scale_colour_manual(
    name = "Method",
    labels = c(
      "power_filtered_fdr"        = "FDR (BH)",
      "power_filtered_by"  = "FDR (BY)",
      "power_filtered_bon"        = "Bonferroni",
      "sa_mvn"        = "SOJO eGenes",
      "me_zeroed" = "Empirical (alpha=0 removed)",
      "an_pm"     = "Analytic PowerMap",
      "sa_perm" = "Permutation",
      "sa_mvn2" = "50% of GTEx-significant eGenes"
    ),
    values = c(
      "power_filtered_fdr" = "#00ba38",
      "power_filtered_by" = "#619cff",
      "power_filtered_bon" = "#f8766d",
      "sa_mvn3" = "magenta",
      "sa_kerby" = "green",
      "sa_perm" = "violetred",
      "me_zeroed" = "chocolate4",
      "sa_kerby_non1" = "blue"
    ),
    guide = guide_legend(override.aes = list(size = 2), nrow = 1)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power to Detect eGenes by Multiple Testing Adjustment",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1))

ggsave("Plots/power_comparisons_newer/power_egenes_multipletesting.pdf", width = 7, height = 4)

