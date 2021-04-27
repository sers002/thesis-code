library(dplyr)
library(ggplot2)
library(data.table)

maf_labeller <- function(maf_bin) {
  sprintf("MAF: %s%%", maf_bin)
}

## Without alpha=0 removed
sa_power     <- fread("Output/power/sa_mvn_eur.csv")
sa_power <- sa_power[qtl_tag == TRUE & cis_flag == TRUE]
# sa_power$sample_size <- sa_power$n

sa_power_afr     <- fread("Output/power/sa_mvn_afr.csv")
sa_power_afr <- sa_power_afr[qtl_tag == TRUE & cis_flag == TRUE]
# sa_power_afr$sample_size <- sa_power_afr$n

sa_power_eas     <- fread("Output/power/sa_mvn_eas.csv")
sa_power_eas <- sa_power_eas[qtl_tag == TRUE & cis_flag == TRUE]
# sa_power_eas$sample_size <- sa_power_eas$n

all_power <- rbindlist(list(
  data.table(sa_power[qtl_tag == TRUE & sample_size != 404],           method = "EUR"),
  data.table(sa_power_afr[qtl_tag == TRUE & sample_size != 404],           method = "AFR"),
  data.table(sa_power_eas[qtl_tag == TRUE & sample_size != 404],           method = "EAS")
), fill = TRUE)

all_power[, maf_bin_labelled := sprintf("'%s'", maf_labeller(maf_bin))]
all_power[, es_bin_labelled := sprintf("alpha*': %s'", es_bin)]


all_power %>% 
  subset(method %in% c("EUR", "AFR", "EAS")) %>%
  ggplot(aes(x = sample_size, y = power_filtered_fdr, colour = method)) +
  geom_point() +
  geom_line() +
  # geom_errorbar(aes(ymin = power_filtered_fdr_lwr, ymax = power_filtered_fdr_upr), position = position_dodge(width = 0.5)) + 
  facet_wrap(~es_bin_labelled, labeller = label_parsed) +
  scale_colour_manual(
    name = "Ancestry",
    labels = c(
      "AFR"        = "African",
      "EUR"  = "European",
      "EAS"        = "East Asian",
      "me"        = "Empirical",
      "sa_zeroed" = "Empirical (alpha=0 removed)",
      "an_pm"     = "Analytic PowerMap",
      "sa_wrapup" = "Semi-analytic (Wrap-up)"
    ),
    values = c(
      "EUR" = "#00ba38",
      "AFR" = "#619cff",
      "EAS" = "#f8766d",
      "an_pm" = "magenta",
      "sa_kerby" = "green",
      "sa_wrapup" = "violetred",
      "sa_zeroed" = "chocolate4",
      "sa_kerby_non1" = "blue"
    ),
    guide = guide_legend(override.aes = list(size = 2), nrow = 1)
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power Curve Comparison Between Ancestries",
    # subtitle = "All causal variants; BH < 0.05",
    x = "Sample Size",
    y = "Power"
  ) +
  guides(x =  guide_axis(angle = -90)) +
  ylim(c(0, 1)) 

# ggsave(filename = "Plots/power_comparisons_newer/power_ancestry", width = 7, height = 5)





