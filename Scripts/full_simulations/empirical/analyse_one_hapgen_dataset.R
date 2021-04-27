###
### THIS SCRIPT SHOULD ONLY BEEN RUN FROM callr WITH A FILENAME ARGUMENT
###

library(data.table)
library(Matrix.utils)
# devtools::load_all("../MatrixEQTL-master")

source(here::here("Scripts/helper_functions.R"))

filename <- commandArgs(trailingOnly = TRUE)

load(file = here::here("Output/sojo_output/all_betas.RData"))

combined <- combined[, c("snp", "alpha", "qtl_tag", "beta_shrink", "pid", "sid", "sid_pos", "b", "b_se", "maf", "cis_flag", "gene_start", "gene_end", "a1", "a2")]

# ## Each repetition will use the same variant/gene locations so only save these once
location_filenames <- write_matrixeqtl_files(combined = combined)

message(" Reading in genotype data")
## Read in genotype data
genotype_data <- read_genotype_hapgen(filename)

## Subset to only SNPs in both genotype and betas
in_both <- intersect(combined$snp, genotype_data$snp)
genotype_data <- genotype_data[snp %in% in_both]

message(" Generate expression")
## Create gene expression matrix
##   Start at column 6 for HAPGEN outputs
expression_data <- generate_expression(
  genotype_data,
  combined,
  genotype_start_col = 6
)

## Write out genotype/gene expression matrices and SNP/gene locations for use
## in MatrixEQTL
matrixeqtl_files <- write_matrixeqtl_files(
  genotype_mat    = genotype_data,
  gene_expression = expression_data,
  directory       = tempdir()
)

## Extract out the N and replicate number for the simulation dataset (so we know
## what to name MatrixEQTL results files)
n <- gsub(".*n([0-9]+).*", "\\1", filename)
i <- gsub(".*_([0-9]+)\\.controls.*", "\\1", filename)

me <- run_matrixeqtl(
  SNP_file_name           = matrixeqtl_files$snp,
  snps_location_file_name = location_filenames$snp_loc,

  expression_file_name    = matrixeqtl_files$expression,
  gene_location_file_name = location_filenames$gene_loc,

  output_file_name_cis    = sprintf("Output/matrixEQTL/rep_sims/n_%s/output_cis_%s.txt.gz", n, i),
  output_file_name_tra    = sprintf("Output/matrixEQTL/rep_sims/n_%s/output_trans_%s.txt.gz", n, i)
)

## Delete the files used for MatrixEQTL
unlink(
  c(
    matrixeqtl_files$snp,
    matrixeqtl_files$expression
  )
)

rm(genotype_data, expression_data, me)
gc()

