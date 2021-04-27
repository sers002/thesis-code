## This file takes PLINK's Hardy-Weinberg results, adjusts the p-values and 
## exports a list of SNPs where the adjusted p is < 0.05

library(data.table)

hardy.dt = fread("Output/quality_control/plink.hwe.gz")
hardy.dt[, p.adjust := p.adjust(P, method = "bonferroni")]

cat(
  hardy.dt[p.adjust < 0.05]$SNP,
  sep = "\n",
  file = "Output/quality_control/bad_hardy_snps.txt"
)







library(data.table)
library(ggplot2)

#### 1000G MAF

og.geno = fread(
  "Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz", 
  skip = 100
)

og.geno = og.geno[nchar(REF) == 1]
og.geno = og.geno[nchar(ALT) == 1]

og.geno = og.geno[!duplicated(ID)]

setnames(
  og.geno,
  c("ID",  "ALT", "REF"),
  c("snp", "a1", "a2")
)

## VCFs have genotypes stored in "x|y" format for each SNP, so these need to be 
## added to get 0, 1, 2 
og.geno[1:5, 1:12]

og.geno[, 10:ncol(og.geno)] <- lapply(
  og.geno[, 10:ncol(og.geno)], 
  function(x) {
    as.numeric(substr(x,1,1)) + as.numeric(substr(x,3,3))
  }
)

og.geno[1:5, 1:12]

og.maf <- rowMeans(og.geno[, 10:ncol(og.geno)]) / 2 
og.var <- apply(og.geno[, 10:ncol(og.geno)], MARGIN = 1, FUN = var)
og.maf.dt <- data.table(
  snp = og.geno$snp, 
  maf_1000g = round(og.maf, 4), 
  var_1000g = og.var,
  a2_1000g = og.geno$a2, 
  a1_1000g = og.geno$a1
)
og.maf.dt <- og.maf.dt[!duplicated(og.maf.dt$snp)]

gtex_snps <- readLines("Output/quality_control/gtex_all_snps.txt")

og.maf.dt <- og.maf.dt[snp %in% gtex_snps]

og.maf.dt[, var_expected := 2 * maf_1000g * (1 - maf_1000g)]

og.maf.dt <- merge(og.maf.dt, hardy.dt, by.x = c("snp"), by.y = "SNP")

ggplot(og.maf.dt, aes(var_1000g, var_expected, colour = p.adjust < 0.05)) + 
  geom_point()

hardy_limit_labels <- c(
  sprintf("Adj. p < 0.05 (N = %s)", sum(og.maf.dt$p.adjust < 0.05)),
  sprintf("Adj. p ≥ 0.05 (N = %s)", prettyNum(sum(og.maf.dt$p.adjust >= 0.05), big.mark = ",", scientific = FALSE))
)

og.maf.dt2 <- rbind(
  og.maf.dt[p.adjust < 0.05],
  dplyr::slice_sample(og.maf.dt[p.adjust >= 0.05], prop = 0.1)
)

ggplot(
  og.maf.dt, 
  aes(
    var_1000g, 
    var_expected, 
    colour = p.adjust >= 0.05, 
    shape = p.adjust >= 0.05,
    alpha = p.adjust >= 0.05
  )
) + 
  geom_point() + 
  scale_color_manual(
    values = c("red", "black"),
    labels = hardy_limit_labels,
    guide = guide_legend(override.aes = list(alpha = 1))
  ) + 
  scale_shape_manual(
    values = c("triangle", "circle"),
    labels = hardy_limit_labels,
    guide = guide_legend(override.aes = list(alpha = 1))
  ) + 
  scale_alpha_manual(
    values = c(
      "FALSE" = 0.25,
      "TRUE" = 0.01
    ),
    guide = NULL
  ) +
  labs(
    title = "Comparison of Observed vs Expected Variance",
    x = "Observed Variance",
    y = "Expected Variance",
    colour = "Hardy-Weinberg Equilibrium",
    shape = "Hardy-Weinberg Equilibrium"
  ) + 
  xlim(0, NA) +
  ylim(0, NA) +
  # coord_fixed() + 
  geom_abline(slope = 1, colour = "blue", lty = "dashed")

ggsave(
  filename = "Plots/quality_control/hardy_weinberg_comparison.png",
  width = 6, 
  height = 4
)






