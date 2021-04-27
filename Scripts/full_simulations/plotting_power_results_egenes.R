library(dplyr)
library(ggplot2)
library(data.table)

maf_labeller <- function(maf_bin) {
  sprintf("MAF: %s%%", maf_bin)
}

## Without alpha=0 removed
me_power     <- fread("Output/power_egenes/matrixEQTL.csv")
# me_power <- me_power[qtl_tag == TRUE]
## With alpha=0 removed
# me_power2     <- fread("Output/power/matrixEQTL2.csv")
me_power$sample_size <- me_power$n
# me_power2$sample_size <- me_power2$n

mvn_power3    <- fread("Output/power_egenes/sa_mvn_nullgenes_anything_not_justin.csv")
mvn_power2    <- fread("Output/power_egenes/sa_mvn_nullgenes_50pct_of_justin.csv")
mvn_power4    <- fread("Output/power_egenes/sa_mvn_nullgenes_random_50.csv")
mvn_power5    <- fread("Output/power_egenes/sa_mvn_nullgenes_random2_25.csv")
mvn_power6    <- fread("Output/power_egenes/sa_mvn_nullgenes_random3_75.csv")

# mvn_power <- fread("Output/power_egenes/sa_mvn_extended2.csv")
mvn_power <- rbind(
  fread("Output/power_egenes/sa_mvn.csv"),
  fread("Output/power_egenes/sa_mvn_extended2.csv")[sample_size != 2000]
)
perm_power <- fread("Output/power_egenes/sa_permall.csv")
perm_power[, power_filtered_fdr := power]

# source("Scripts/full_simulations/3_analytic/analytic_power_fdr.R")

all_power <- rbindlist(list(
  data.table(me_power[!(n %in% c(404, 100))],     method = "me"),
  data.table(mvn_power[egene == TRUE],    method = "sa_mvn"),
  data.table(mvn_power2[egene == TRUE],    method = "sa_mvn2"),
  data.table(mvn_power3[egene == TRUE],    method = "sa_mvn3"),
  data.table(mvn_power4[egene == TRUE],    method = "sa_mvn4"),
  data.table(mvn_power5[egene == TRUE],    method = "sa_mvn5"),
  data.table(mvn_power6[egene == TRUE],    method = "sa_mvn6"),
  data.table(perm_power, method = "sa_perm")
), fill = TRUE)
  
all_power %>% 
  filter(method %in% c("sa_mvn", "sa_mvn2", "sa_mvn3")) %>% 
  ggplot(aes(x = sample_size, y = power_filtered_fdr, colour = method)) +
  geom_point() +
  geom_line() +
  # facet_grid(es_bin~maf_bin, labeller = labeller(maf_bin = maf_labeller)) +
  scale_colour_manual(
    name = "True eGenes",
    labels = c(
      "sa_mvn3"        = "Sig. in GTEx After Hi-C Filtering",
      "sa_mvn"        = "SOJO eGenes",
      "sa_mvn2" = "Sig. in GTEx After Hi-C Filtering (50% Sample)"
    ),
    values = c(
      "sa_mvn" = "#00ba38",
      "sa_mvn2" = "#619cff",
      "sa_mvn3" = "magenta"
    ),
    guide = guide_legend(override.aes = list(size = 2), nrow = 2)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power to Detect eGenes by Proportion of True eGenes",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1))

 ggsave("Plots/power_comparisons_newer/power_egenes_proportions.pdf", width = 7, height = 4)

 all_power %>% 
   filter(method %in% c("sa_mvn", "sa_mvn4", "sa_mvn5", "sa_mvn6")) %>% 
   ggplot(aes(x = sample_size, y = power_filtered_fdr, colour = method)) +
   geom_point() +
   geom_line() +
   # facet_grid(es_bin~maf_bin, labeller = labeller(maf_bin = maf_labeller)) +
   scale_colour_manual(
     name = "True eGenes",
     labels = c(
       "sa_mvn5"        = "SOJO eGenes (75% Sample)",
       "sa_mvn6"        = "SOJO eGenes (25% Sample)",
       "sa_mvn"        = "SOJO eGenes",
       "sa_mvn4" = "SOJO eGenes (50% Sample)"
     ),
     values = c(
       "sa_mvn" = "#00ba38",
       "sa_mvn5" = "#619cff",
       "sa_mvn4" = "magenta",
       "sa_mvn6" = "red"
     ),
     guide = guide_legend(override.aes = list(size = 2), nrow = 2)
   ) +
   theme(legend.position = "bottom") +
   labs(
     title = "Power to Detect eGenes by Proportion of True eGenes",
     x = "Sample Size",
     y = "Power"
   ) +
   guides(x =  guide_axis(angle = -90)) +
   ylim(c(0, 1))
 
 ggsave("Plots/power_comparisons_newer/power_egenes_proportions_random.pdf", width = 7, height = 4)
 
 
 
 
 
 
 


all_power %>% 
  filter(method %in% c("sa_mvn", "me")) %>% 
  ggplot(aes(x = sample_size, y = power_filtered_fdr, colour = method)) +
  geom_point() +
  geom_line() +
  # facet_grid(es_bin~maf_bin, labeller = labeller(maf_bin = maf_labeller)) +
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
  ggplot(aes(x = sample_size, y = value, colour = variable, linetype = variable, shape = variable)) +
  geom_point() +
  geom_line() +
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
  rbind(perm_power[, .(variable = "perm", value = power), by = "sample_size"], fill = TRUE) %>% 
  ggplot(aes(x = sample_size, y = value, colour = variable)) +
  geom_point() +
  geom_line() +
  # facet_grid(es_bin~maf_bin, labeller = labeller(maf_bin = maf_labeller)) +
  scale_colour_manual(
    name = "Multiple Testing\nAdjustment",
    labels = c(
      "power_filtered_fdr" = "FDR (BH)",
      "power_filtered_by"  = "FDR (BY)",
      "power_filtered_bon"        = "Bonferroni",
      "perm"        = "Permutation (eGene-MVN)"
    ),
    values = c(
      "power_filtered_fdr" = "#6C803A",
      "power_filtered_by" = "#AB7C47",
      "power_filtered_bon" = "#CCAE42",
      "perm" = "#D73202"
    ),
    guide = guide_legend(override.aes = list(size = 2), nrow = 2)
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



library(dplyr)
library(ggplot2)
sa_power_filtered <- fread("Output/power_egenes/sa_mvn_filtered.csv")

sa_power_filtered[egene == TRUE] %>% 
  melt(
    measure.vars = c("power_filtered_fdr", "power_all_fdr", "power_all_fdr_same", "power_all_fdr_alt", "power_all_fdr_alt2", "power_all_fdr_alt3", "power_all_fdr_alt4"),
    variable.name = "Adjustment",
    value.name = "power"
  ) %>% 
  ggplot(aes(sample_size, power, colour = Adjustment)) +
  geom_point() + 
  geom_line() + 
  scale_colour_discrete(
    name = "Filtering",
    labels = c(
      "power_filtered_fdr"        = "Hi-C Filtered cis Variants",
      "power_all_fdr"  = "All cis Variants (Null Distribution)", 
      "power_all_fdr_alt"  = "All cis Variants (20% Distribution)", 
      "power_all_fdr_alt2"  = "All cis Variants (10% Distribution)", 
      "power_all_fdr_alt3"  = "All cis Variants (5% Distribution)", 
      "power_all_fdr_alt4"  = "All cis Variants (1% Distribution)", 
      "power_all_fdr_same"  = "All cis Variants (Same Distribution)"
    ),
    # values = c(
    #   "power_filtered_fdr" = "#A7473A",
    #   "power_all_fdr" = "#B09B37",
    #   "power_all_alt_fdr" = "#4b5f6c", 
    #   "power_all_same_fdr" = "#955f47"
    # ),
    guide = guide_legend(override.aes = list(size = 2, linetype = NA), nrow = 3)
  ) +
  scale_linetype_discrete(
    name = "Filtering",
    labels = c(
      "power_filtered_fdr"        = "Hi-C Filtered cis Variants",
      "power_all_fdr"  = "All cis Variants (Null Distribution)", 
      "power_all_fdr_alt"  = "All cis Variants (20% Distribution)", 
      "power_all_fdr_alt2"  = "All cis Variants (10% Distribution)", 
      "power_all_fdr_alt3"  = "All cis Variants (5% Distribution)", 
      "power_all_fdr_alt4"  = "All cis Variants (1% Distribution)", 
      "power_all_fdr_same"  = "All cis Variants (Same Distribution)"
    )
  ) + 
  scale_shape_discrete(
    name = "Filtering",
    labels = c(
      "power_filtered_fdr"        = "Hi-C Filtered cis Variants",
      "power_all_fdr"  = "All cis Variants (Null Distribution)", 
      "power_all_fdr_alt"  = "All cis Variants (20% Distribution)", 
      "power_all_fdr_alt2"  = "All cis Variants (10% Distribution)", 
      "power_all_fdr_alt3"  = "All cis Variants (5% Distribution)", 
      "power_all_fdr_alt4"  = "All cis Variants (1% Distribution)", 
      "power_all_fdr_same"  = "All cis Variants (Same Distribution)"
    )
  ) + 
  theme(legend.position = "bottom") +
  labs(
    title = "Power to Detect cis-eGenes by Hi-C Filtering",
    subtitle = "BH FDR < 0.05",
    x = "Sample Size",
    y = "Power"
  ) +
  ylim(c(0, 1))

ggsave("Plots/power_comparisons_newer/power_egenes_filtering.pdf", width = 7.8, height = 4)


all_power %>% 
  filter(method %in% c("sa_mvn", "me")) %>% 
  filter(!(method == "me" & sample_size == 2000)) %>% 
  ggplot(aes(x = sample_size, y = power_filtered_fdr)) +
  geom_point(colour = "#5FA1F7") +
  geom_line(colour = "#5FA1F7") +
  theme(legend.position = "bottom") +
  labs(
    title = "Comparison of Power to Detect eGenes",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1)) 

ggsave("Plots/power_comparisons_newer/power_egenes_extended.pdf", width = 7, height = 4)

