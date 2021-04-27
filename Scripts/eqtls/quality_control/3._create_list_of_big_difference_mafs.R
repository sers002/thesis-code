## This file identifies SNPs where the AF is significantly different between the
## 1000G estimate and GTEx estimate (>10%).

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
og.maf.dt <- data.table(
  snp = og.geno$snp,
  maf_1000g = round(og.maf, 4),
  a2_1000g = og.geno$a2,
  a1_1000g = og.geno$a1
)
og.maf.dt <- og.maf.dt[!duplicated(og.maf.dt$snp)]


#### GTEX MAFS

gtex_data <- fread("Data/lung_eqtl_data/eqtls.txt")
gtex_data = gtex_data[!duplicated(gtex_data[, c("pid", "sid")])]
genes = fread("Data/lung_eqtl_data/genes.txt")
gtex_data = merge(gtex_data, genes[gene_chr == "chr19", ], by.x = c("sid", "pid"), by.y = c("variant_id","gencode_id"))
rm(genes)

gtex.mafs <- unique(
  gtex_data[, list(snp = snp, maf_gtex = round(maf, 4))]
)

anyDuplicated(og.maf.dt$snp)
anyDuplicated(gtex.mafs$snp)

#### Merging GTEx and 1000G

mafs.dt <- merge(og.maf.dt, gtex.mafs, by = "snp")
mafs.dt[, maf_diff := maf_1000g - maf_gtex]
mafs.dt[, maf_add := maf_1000g + maf_gtex]

# mafs.dt[abs(maf_add - 1) < 0.05 & !abs(round(maf_diff, 3)) < 0.1 ]

## SNPs with MAFs within 10% between GTEx and 1000G
mafs.dt[, within_range := abs(maf_1000g - maf_gtex) < 0.1 | abs((1 - maf_1000g) - maf_gtex) < 0.1]

## Export list of SNPs where MAFs differ by 10% or more
cat(
  mafs.dt[within_range == FALSE]$snp,
  sep = "\n",
  file = "Output/quality_control/out_of_range_snps.txt"
)


# Frequency comparison plots ----------------------------------------------

af_limit_labels <- c(
  sprintf("F (N = %s)", sum(mafs.dt$within_range == FALSE)),
  sprintf("T (N = %s)", prettyNum(sum(mafs.dt$within_range == TRUE), big.mark = ",", scientific = FALSE))
)

ggplot(mafs.dt, aes(
  maf_1000g,
  maf_gtex,
  colour = within_range,
  shape = within_range,
  alpha = within_range
  )
) +
  geom_point() +
  scale_color_manual(
    values = c("red", "black"),
    labels = af_limit_labels,
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  scale_shape_manual(
    values = c("triangle", "circle"),
    labels = af_limit_labels,
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
    title = "Comparison of Allele Frequencies",
    x = "AF (1000 Genomes)",
    y = "AF (GTEx)",
    colour = "Within 10% AF",
    shape = "Within 10% AF"
  )

ggsave(filename = "Plots/quality_control/allele_freq_gtex_vs_1000g.png", width = 6, height = 4)
