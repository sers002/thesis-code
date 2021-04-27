## Load libraries
library(data.table)
library(ggplot2)


# GTEx Data ---------------------------------------------------------------

## Read in CoDeS3D data - all associations of snps on Chr 19
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
genes = fread(here::here("Data/lung_eqtl_data/genes.txt"))

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


# 1000G Data --------------------------------------------------------------

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


# Merging 1000G and GTEx --------------------------------------------------

file_match = merge(ss, ref_data, by.x = "snp", by.y = "ID")

## Checks
# sum(is.na(file_match$REF))
# length(unique(ref_data$ID))
# length(unique(ss$snp))
# length(unique(file_match$snp))

# Check direction of alleles between 1000G and GTEx
file_match[, orig1 := paste0(a1, a2)]
file_match[, new1  := paste0(REF, ALT)]
file_match[, new2  := paste0(ALT, REF)]

## 4 rows have no match between 1000G and GTEx, so these are removed
# nrow(file_match[(orig1 != new1) & (orig1 != new2)])
file_match = file_match[!((orig1 != new1) & (orig1 != new2))]

## All match direction of a1 = ALT, a2 = REF
# nrow(file_match[(orig1 == new1)])
# nrow(file_match[(orig1 == new2)])

# length(unique(file_match$sid))
# length(unique(file_match$pid))
# length(unique(ss$sid))
# length(unique(ss$pid))


# Quality Control ---------------------------------------------------------

## Remove SNPs where Bonferroni-adjusted HWE p-value is < 0.05
bad_hardy <- readLines(here::here("Output/quality_control/bad_hardy_snps.txt"))
file_match <- file_match[!(snp %in% bad_hardy)]

## Remove SNPs where allele frequency differs by more than 10%
big_diff <- readLines(here::here("Output/quality_control/out_of_range_snps.txt"))
file_match <- file_match[!(snp %in% big_diff)]

rm(bad_hardy, big_diff)

## Remove columns not being used
file_match[, orig1 := NULL]
file_match[, new1 := NULL]
file_match[, new2 := NULL]
file_match[, sid_chr := NULL]
file_match[, tissue := NULL]
file_match[, chrom := NULL]
file_match[, locus := NULL]
file_match[, gene_chr := NULL]
file_match[, interactions := NULL]
file_match[, replicates := NULL]
file_match[, cell_line := NULL]
file_match[, cell_line_replicates := NULL]
file_match[, sum_interactions := NULL]
file_match[, sum_replicates := NULL]
file_match[, sum_cell_lines := NULL]

rm(ref_data)


# Gene Investigations -----------------------------------------------------

file_match[, snps_per_gene := length(sid), by = pid]
# summary(file_match$snps_per_gene)

file_match[, gene_length    := gene_end - gene_start]
file_match[, q2_raw         := (b^2)*2*maf*(1-maf)]
file_match[, h2_raw         := sum(q2_raw), by = pid]
file_match[, n_sig_per_gene := sum(adj_pval < 0.05), by = pid]

gene_heritability = unique(file_match[, c("pid", "n_sig_per_gene",  "snps_per_gene", "gene_length", "h2_raw")])

gene_heritability[, snp_density := gene_length/snps_per_gene]
# summary(gene_heritability$snp_density)

# summary(gene_heritability$h2_raw)

gene_heritability[, avg_q2 := h2_raw/snps_per_gene]
# summary(gene_heritability$avg_q2)


# Heritability Quartiles --------------------------------------------------

sumq = summary(gene_heritability$avg_q2)
sumq = unname(sumq)
gene_heritability[, avg_q2_Q := "Q1"]
gene_heritability[avg_q2 > sumq[2], avg_q2_Q := "Q2"]
gene_heritability[avg_q2 > sumq[3], avg_q2_Q := "Q3"]
gene_heritability[avg_q2 > sumq[5], avg_q2_Q := "Q4"]

file_match = merge(
  file_match,
  gene_heritability[, c("pid", "avg_q2", "avg_q2_Q")],
  by = "pid",
  all.x = TRUE
)

# summary(gene_heritability[causalpropQ == "Q4"]$n_sig_per_gene)
# summary(gene_heritability[avg_q2_Q == "Q4"]$n_sig_per_gene)
# summary(gene_heritability[avg_q2_Q == "Q4"]$causalprop)
# summary(gene_heritability[avg_q2_Q == "Q4"]$avg_q2)

rm(sumq, gene_heritability)
