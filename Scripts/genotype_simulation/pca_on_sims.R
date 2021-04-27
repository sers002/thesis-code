library(data.table)

# HAPGEN ------------------------------------------------------------------

## Import the first 5 columns of the gen file (metadata about variants)
hapgen.gen <- fread(
  "./Output/genotype_simulation/1000g_chr19_eur_simulation.controls.gen", 
  select = 1:5
)

hapgen.haps <- data.table::fread(
  "./Output/genotype_simulation/1000g_chr19_eur_simulation.controls.haps"
)

## Add consecutive pairs of columns (i.e. col1+col2, col3+col4, col5+col6, ...)
hapgen.dip <- hapgen.haps[, .SD, .SDcols=seq(1, ncol(hapgen.haps) - 1, 2)] + hapgen.haps[, .SD, .SDcols=seq(2, ncol(hapgen.haps), 2)]

hapgen.dip <- cbind(hapgen.gen, hapgen.dip)

rm(
  hapgen.haps,
  hapgen.gen
)

gc()

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

in.both <- intersect(original.vcf$ID, hapgen.dip$V2)

original.vcf <- original.vcf[ID %in% in.both]
hapgen.dip   <- hapgen.dip[V2 %in% in.both]

original.vcf <- original.vcf[!duplicated(ID)]
hapgen.dip   <- hapgen.dip[!duplicated(V2)]

setorder(original.vcf, ID)
setorder(hapgen.dip, V2)

## Check SNPs are in the same order
all.equal(hapgen.dip$V2, original.vcf$ID)

## Do the SNPs need flipping?
joined.df <- original.vcf[, 1:9][hapgen.dip[, 1:5], on = c("ID" = "V2")]
joined.df[, match1 := REF == V4 & ALT == V5]
joined.df[, match2 := REF == V5 & ALT == V4]
table(joined.df$match2)

library(Matrix)
library(mltools)

original.mat <- t(sparsify(original.vcf[, 10:ncol(original.vcf)]))
simulated.mat <- t(sparsify(2 - hapgen.dip[, 6:ncol(hapgen.dip)]))

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


# plot(
#   evv$vectors[, 1], 
#   evv$vectors[, 2], 
#   col = rep(1:2, c(n.original, n.simulated))
# )

library(ggplot2)
library(patchwork)
library(dplyr)

eigen.df <- data.frame(
  population = rep(c("Raw", "Simulated"), c(n.original, n.simulated)),
  evv$vectors[, 1:4]
)
colnames(eigen.df) = c("population", "PC1", "PC2", "PC3", "PC4")

eigen.df <- eigen.df %>% 
  group_by(population) %>% 
  sample_n(sum(eigen.df$population == "Raw"))

eigen.plot1 <- ggplot(eigen.df, aes(x = PC1, y = PC2, colour = population)) + 
  geom_point(alpha = 0.5) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  labs(
    colour = "Dataset"
  )

eigen.plot1

eigen.plot2 <- ggplot(eigen.df, aes(x = PC3, y = PC4, colour = population)) + 
  geom_point(alpha = 0.5) + 
  theme_bw() + 
  labs(
    colour = "Dataset"
  )

eigen.plot2

eigen.plot1 + eigen.plot2

