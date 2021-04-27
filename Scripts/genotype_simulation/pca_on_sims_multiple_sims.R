library(data.table)

# HAPGEN ------------------------------------------------------------------

source(here::here("Scripts/full_simulations/investigation/helper_functions.R"))

hapgen.dip1 <- read_genotype_hapgen("Output/genotype_simulation/rep_sims/n_404/1000g_sim_n404_1.controls.haps.gz")
hapgen.dip2 <- read_genotype_hapgen("Output/genotype_simulation/rep_sims/n_404/1000g_sim_n404_2.controls.haps.gz")

# OG ----------------------------------------------------------------------

original.vcf <- fread("./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz", skip = 100)

## VCFs have "0|0" for the genotyep for each participant, so turn these into 0,1,2
original.vcf[, 10:ncol(original.vcf)] <- lapply(
  original.vcf[, 10:ncol(original.vcf)],
  function(col) {
    as.numeric(substr(col, 1, 1)) + as.numeric(substr(col, 3, 3))
  }
)

# Combining ---------------------------------------------------------------

all.equal(hapgen.dip1$snp, hapgen.dip2$snp)

hapgen.dip <- cbind(hapgen.dip1, hapgen.dip2[, -(1:5)])

rm(hapgen.dip1, hapgen.dip2)

in.both <- intersect(original.vcf$ID, hapgen.dip$snp)

original.vcf <- original.vcf[ID %in% in.both]
hapgen.dip   <- hapgen.dip[snp %in% in.both]

original.vcf <- original.vcf[!duplicated(ID)]
hapgen.dip   <- hapgen.dip[!duplicated(snp)]

setorder(original.vcf, ID)
setorder(hapgen.dip, snp)

## Check SNPs are in the same order
all.equal(hapgen.dip$V2, original.vcf$snp)

## Do the SNPs need flipping?
joined.df <- original.vcf[, 1:9][hapgen.dip[, 1:5], on = c("ID" = "snp")]
joined.df[, match1 := REF == a1 & ALT == a2]
joined.df[, match2 := REF == a2 & ALT == a1]
table(joined.df$match2)

library(Matrix)
library(mltools)

original.mat <- t(sparsify(original.vcf[, 10:ncol(original.vcf)]))
simulated.mat <- t(sparsify(hapgen.dip[, 6:ncol(hapgen.dip)]))

dim(original.mat)
dim(simulated.mat)

rm(
  original.vcf,
  hapgen.dip,
  joined.df
)

gc()

n.original <- nrow(original.mat)
n.simulated <- nrow(simulated.mat)

all.mat <- rbind(original.mat, simulated.mat)
all.maf <- colMeans(all.mat) / 2

all.mat <- all.mat[, all.maf %between% c(0.1, 0.9)]
all.maf <- all.maf[all.maf %between% c(0.1, 0.9)]

cols <- sample(1:ncol(all.mat), size = 50000)
# cols <- 1:ncol(all.mat)

all.mat2 <- all.mat[, cols]
all.maf2 <- all.maf[cols]

## Standardise genotype matrix
all_stand.mat <- sweep(all.mat2, 2, all.maf2)
all_stand.mat <- sweep(all_stand.mat, 2, sqrt(all.maf2 * (1 - all.maf2)), FUN = "/")

xxmat <- tcrossprod(all_stand.mat)
evv <- eigen(xxmat, symmetric=TRUE)

all_stand.pc <- prcomp(all_stand.mat)
all_stand.pc2 <- princomp(all_stand.mat)

# plot(
#   evv$vectors[, 1],
#   evv$vectors[, 2],
#   col = rep(1:2, c(n.original, n.simulated))
# )

library(ggplot2)
library(patchwork)
library(dplyr)

eigen.df <- data.frame(
  population = rep(c("Raw", "Simulation 1", "Simulation 2"), c(n.original, n.simulated / 2, n.simulated / 2)),
  all_stand.pc$x[, 1:4]
)
# colnames(eigen.df) = c("population", "PC1", "PC2", "PC3", "PC4")
# eigen.df$PC1 <- eigen.df$PC1 * sqrt(evv$values[1])
# eigen.df$PC2 <- eigen.df$PC2 * sqrt(evv$values[2])

eigen.df <- eigen.df %>%
  group_by(population) %>%
  sample_n(sum(eigen.df$population == "Raw"))

eigen.plot1 <- ggplot(eigen.df, aes(x = PC1, y = PC2, colour = population)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(
    colour = "Dataset",
    title  = "Comparison of Genotype Datasets"
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  coord_fixed()

eigen.plot1

# ggsave("Plots/genotype_simulation/pca_comparison.pdf", width = 200, height = 150, units = "mm", dpi = "print")
ggsave("Plots/genotype_simulation/pca_comparison.pdf", width = 132, height = 99, units = "mm", dpi = "print")

eigen.plot2 <- ggplot(eigen.df, aes(x = PC3, y = PC4, colour = population)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(
    colour = "Dataset"
  )

eigen.plot2

eigen.plot1 + eigen.plot2

