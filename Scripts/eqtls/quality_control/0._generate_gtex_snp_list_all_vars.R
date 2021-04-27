## Load libraries
library(data.table)
library(ggplot2)

## Read in Justin's data - all associations of snps on Chr 19
eqtls = fread(here::here("Data/lung_eqtl_data/eqtls.txt"))

## Duplicate checking
# nrow(unique(eqtls[, c("sid", "pid")])) == nrow(eqtls)

# test = eqtls[duplicated(eqtls[, c("sid", "pid")])]
# test[, max_b := max(b), by = c("sid", "pid")]
# test[, min_b := min(b), by = c("sid", "pid")]
# test[, range_b := max_b - min_b]
# test2 = unique(test[, c("pid", "sid", "sid_pos", "min_b", "max_b", "range_b")])
# setorder(test2, pid, sid_pos)
# summary(test2$range_b)
# rm(test, test2)

## Difference in b's is negligible - pick a record at random for duplicates
eqtls = eqtls[!duplicated(eqtls[, c("pid", "sid")])]

## Read in Justin's data - gene information
genes = fread("Data/lung_eqtl_data/genes.txt")

## Duplicate checking
# nrow(unique(genes[, c("variant_id", "gencode_id")])) == nrow(genes)

## Merge to get all cis and intrachromosomal trans associations on Chr19
ss = merge(eqtls, genes[gene_chr == "chr19", ], by.x = c("sid", "pid"), by.y = c("variant_id","gencode_id"))

rm(eqtls, genes)

# add cis/trans flag
ss[, cis_flag := FALSE]
ss[sid_pos >= (gene_start - 1000000 - 1) & sid_pos <= (gene_end + 1000000), cis_flag := TRUE]
table(ss$cis_flag, exclude = NULL)

ss$chr <- 19

## Extract out alt allele
ss$a1 <- gsub("[A-z0-9]+_[0-9]+_([A-Z]+)_([A-Z]+).*", "\\2", ss$sid)

## Extract out ref allele
ss$a2 <- gsub("[A-z0-9]+_[0-9]+_([A-Z]+)_.*", "\\1", ss$sid)

## Subset to SNPs only
ss = ss[nchar(a1) == 1]
ss = ss[nchar(a2) == 1]

## Read in metadata about variants from 1000G
ref_data = fread(
  here::here("Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz"), 
  skip = 100, 
  select = 2:5
)

## Subset to SNPs only
ref_data = ref_data[nchar(REF) == 1]
ref_data = ref_data[nchar(ALT) == 1]

## Some IDs are duplicated - take the first ID in each case where this happens
ref_data = ref_data[!duplicated(ID)]

file_match = merge(ss, ref_data, by.x = "snp", by.y = "ID")

snp_list = unique(file_match$snp)

cat(
  snp_list, 
  file = here::here("Output/quality_control/gtex_all_snps.txt"),
  sep = "\n"
)
