library(data.table)
library(Matrix.utils)

# Read in declumped eqtls
file_match = fread(here::here("Output/eqtls_for_sojo/declumped_eqtls.csv"))


# 1000 Genomes LD matrix --------------------------------------------------

## PLINK command used in newest-project/Scripts/linkage_disequilibrium
ld.df <- fread(here::here("Output/linkage_disequilibrium/ld_gtex_snps.ld.gz"))

## Check LD with Swedish data from lassosum
# load(file = "Data/LD_chr19.rda")
#
# in.both <- ld.df[SNP_A %in% rownames(LD_mat) & SNP_B %in% colnames(LD_mat)]
# in.both$R_SWEDE <- NA
# for (i in 1:NROW(in.both))
#   in.both$R_SWEDE[i] <- with(in.both, LD_mat[SNP_A[i], SNP_B[i]])
#
# ggplot(in.both, aes(x = abs(R), y = abs(R_SWEDE))) + geom_point()

## Extract all SNP names from LD file that are also in the eQTL data
snp.names <- unique(c(ld.df$SNP_A, ld.df$SNP_B))
snp.names <- snp.names[snp.names %in% file_match$snp]

## Take only the SNPs in the eQTL data
ld.df <- subset(ld.df, SNP_A %in% snp.names & SNP_B %in% snp.names)

## Add in SNPs that are in GTEx that are not in LD (and make zero)
snp.names2 <- setdiff(file_match$snp, snp.names)
snp.names3 = c(snp.names, snp.names2)

## Make sure SNPs with no LD values are still in the matrix later by changing
## to a factor (some will have no observations)
ld.df$SNP_A <- factor(ld.df$SNP_A, levels = snp.names3)
ld.df$SNP_B <- factor(ld.df$SNP_B, levels = snp.names3)

# length(levels(ld.df$SNP_A))
# length(levels(ld.df$SNP_B))

## Pivot to a sparse matrix
ld_mat <- dMcast(
  ld.df,
  SNP_A ~ SNP_B,
  value.var = "R",
  fun.aggregate = "sum",
  drop.unused.levels = FALSE
)

## dMcast puts the variable name in the column names, so remove it to give names
## like rs12345 (vs SNP_Brs12345)
colnames(ld_mat) <- gsub("^SNP_B", "", colnames(ld_mat))

## Some SNPs are not contained within the SNP_A column of the LD dataset, so
## these have to be appended onto the matrix
missing.rows <- snp.names3[which(!(snp.names3 %in% rownames(ld_mat)))]

missing.mat <- Matrix(
  data = 0,
  nrow = length(missing.rows),
  ncol = ncol(ld_mat),
  dimnames = list(missing.rows, colnames(ld_mat)),
  sparse = TRUE
)

ld_mat <- rbind(ld_mat, missing.mat)
# dim(ld_mat)

## Subset to only SNPs in the LD matrix
snp_ref.df <- unique(file_match[snp %in% snp.names3, c("snp", "sid")])

## Set reference allele from GTEX SNP name
snp_ref.df$ref <- gsub("[A-z0-9]+_[0-9]+_([A-Z]+)_.*", "\\1", snp_ref.df$sid)

## SOJO needs a named vector in the form of "rs12345" = "G" where "G" is the
## reference allele
snp_ref2 <- setNames(snp_ref.df$ref, snp_ref.df$snp)

## Checks
# table(names(snp_ref2) %in% colnames(ld_mat))

## Reorder the reference vector
snp_ref2 <- snp_ref2[colnames(ld_mat)]

## Reorder the LD matrix to be the same as the reference vector
ld_mat <- ld_mat[names(snp_ref2), ]
ld_mat <- ld_mat[, names(snp_ref2)]

## Check they have the same ordering
# snp_ref2[1:5]
# ld_mat[1:5, 1:5]

rm(ld.df, missing.mat, missing.rows, snp_ref.df, snp.names, snp.names2, snp.names3)

## Save LD matrix and reference vector for use in SOJO
saveRDS(ld_mat, here::here("Output/linkage_disequilibrium/ld_mat.rds"))
saveRDS(snp_ref2, here::here("Output/linkage_disequilibrium/snp_ref2.rds"))
