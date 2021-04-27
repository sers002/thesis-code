library(data.table)
library(Matrix.utils)

source(here::here("Scripts/eqtls/dataprep_and_exploratory/1._process_eqtl_data.R"))

## Read in LD matrix 
##   (created by Scripts/linkage_disequilibrium/plink_ld_process)
ld.df <- fread(here::here("Output/linkage_disequilibrium/ld_gtex_snps.ld.gz"))

## Read in set of declumped SNPs
##   (created by Scripts/variant_pruning/plink_pruning_gene)
declump = fread(here::here("Output/variant_pruning/plink_gene_prune.in"), header = FALSE)

setnames(declump, "V1", "pid")
setnames(declump, "V2", "snp")

uniqueN(declump$pid)
uniqueN(declump$snp)

## Remove pairwise LD, so only one of each pair remains
ld.df[, ID := paste(pmin(SNP_A, SNP_B), pmax(SNP_A, SNP_B)) ]
ld.df = ld.df[!duplicated(ID)]

## Subset to only SNPs within the 1000G/GTEx merged dataset
ld.df = ld.df[SNP_A %in% file_match$snp & SNP_B %in% file_match$snp]

length(unique(c(ld.df$SNP_A, ld.df$SNP_B)))
length(unique(file_match$snp))

file_match = file_match[snp %in% c(ld.df$SNP_A, ld.df$SNP_B)]

## Subset to only associations in declumped set
file_match = merge(file_match, declump, by = c("snp", "pid"))

## Remove any duplicate betas (randomly picks first) because of LD
file_match = file_match[!duplicated(file_match[, list(pid = pid, b = abs(b))])]

ld.df = ld.df[SNP_A %in% file_match$snp & SNP_B %in% file_match$snp]
length(unique(c(ld.df$SNP_A, ld.df$SNP_B)))
length(unique(file_match$snp))

file_match[, h2_pruned := sum(q2_raw), by = pid]

fwrite(file_match, here::here("Output/eqtls_for_sojo/declumped_eqtls.csv"))

