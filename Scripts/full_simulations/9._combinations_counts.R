library(data.table)
library(ggplot2)
library(Matrix.utils)

load(file = here::here("Output/sojo_output/all_betas.RData"))

test_genotype = fread(here::here("Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz"), skip = 100)

setnames(
  test_genotype,
  c("ID",  "ALT", "REF"),
  c("snp", "a1", "a2")
)

#check for duplicates and remove
to_remove = test_genotype[duplicated(snp)]$snp
test_genotype = test_genotype[!snp %in% to_remove]
combined = combined[!snp %in% to_remove]

rm(to_remove)

#subset to biallelic snps only
test_genotype = test_genotype[nchar(a1) == 1]
test_genotype = test_genotype[nchar(a2) == 1]

test_genotype = test_genotype[, 2:5]

genes = unique(combined[, c("pid", "gene_start", "gene_end")])


genes[, lower_end := gene_start - 1000000 - 1]
genes[, upper_end := gene_end + 1000000]
combined = test_genotype[genes, on = list(POS >= lower_end, POS <= upper_end), nomatch = 0]

#cis
nrow(combined)
#11102563

#trans
nrow(genes)*nrow(test_genotype) - nrow(combined)
#301177637

