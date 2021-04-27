#load libraries
library(data.table)

source("Scripts/eqtls/dataprep_and_exploratory/1._process_eqtl_data.R")

snp_list = unique(file_match$snp)

snp.df <- file_match[, c("pid", "snp")]
setorder(snp.df, pid, snp)

fwrite(
  snp.df, 
  file = "./Output/variant_pruning/snp_list_from_gtex2.txt",
  sep = " ",
  col.names = FALSE
)
