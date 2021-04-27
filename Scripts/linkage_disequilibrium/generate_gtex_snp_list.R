#load libraries
library(data.table)

source("Scripts/eqtls/dataprep_and_exploratory/1._process_eqtl_data.R")

snp_list = unique(file_match$snp)

cat(
  snp_list, 
  file = "./Output/linkage_disequilibrium/snp_list_from_gtex.txt",
  sep = "\n"
)