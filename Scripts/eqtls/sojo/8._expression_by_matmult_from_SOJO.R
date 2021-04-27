library(data.table)
library(ggplot2)
library(Matrix.utils)

load(file = here::here("Output/sojo_output/all_betas.RData"))

combined[is.na(beta_shrink), beta_shrink := 0]


# Genotype ----------------------------------------------------------------

test_genotype = fread(here::here("Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz"), skip = 100)

setnames(
  test_genotype,
  c("ID",  "ALT", "REF"),
  c("snp", "a1", "a2")
)

## Subset to only SNPs common to eQTL dataset and genotype
in_both <- intersect(combined$snp, test_genotype$snp)
test_genotype <- test_genotype[snp %in% in_both]
combined <- combined[snp %in% in_both]

## See if any SNPs need to be flipped
genotype_flip.df <- merge(
  test_genotype[, c("snp", "a1", "a2")],
  unique(combined[, c("snp", "a1", "a2")]),
  by = "snp"
)

nrow(genotype_flip.df[a1.x != a1.y])
to_remove = genotype_flip.df[a1.x != a1.y]$snp

test_genotype = test_genotype[!snp %in% to_remove]
combined = combined[!snp %in% to_remove]

# Check for duplicates and remove
to_remove = test_genotype[duplicated(snp)]$snp
test_genotype = test_genotype[!snp %in% to_remove]
combined = combined[!snp %in% to_remove]

rm(to_remove, in_both, genotype_flip.df)

## Change combined data frame to be a sparse matrix with rows=genes and cols=snps
## with values being the betas
all_assoc_genes <- dMcast(
  combined[, c("snp", "pid", "beta_shrink")],
  pid ~ snp,
  value.var = "beta_shrink",
  fun.aggregate = "sum"
)

## dMcast appends the variable name to column names, so remove the "snp"
colnames(all_assoc_genes) <- sub("^snp", "", colnames(all_assoc_genes))
all_assoc_genes <- all_assoc_genes[, order(colnames(all_assoc_genes))]

## Reorder SNPs so they match up with gene-snp beta matrix
test_genotype <- test_genotype[order(snp)]

## All SNPs are in same order between the two datasets
all.equal(test_genotype$snp, colnames(all_assoc_genes))

dim(test_genotype)

## Convert the "0 | 0" from the VCF to a number
test_genotype[, 10:ncol(test_genotype)] <- lapply(
  test_genotype[, 10:ncol(test_genotype)],
  function(x) {
    as.numeric(substr(x, 1, 1)) + as.numeric(substr(x, 3, 3))
  }
)

## Rearrange datasets to get required matrices
genotype_mat <- Matrix(as.matrix(test_genotype[, 10:ncol(test_genotype)]), sparse = TRUE)
genotype_mat <- t(genotype_mat)

colnames(genotype_mat) <- test_genotype$snp
all_assoc_genes <- t(all_assoc_genes)

## Check dimensions are as required for matrix multiplication
# dim(genotype_mat)
# dim(all_assoc_genes)

gene_expression <- genotype_mat %*% all_assoc_genes

dim(gene_expression)

set.seed(12345)
gene_expression = gene_expression + rnorm(n = ncol(gene_expression)*nrow(gene_expression), 0, 1)

## Save expression and genotype for MatrixEQTL to use
g_out = t(genotype_mat)
g_out = data.table(id = row.names(g_out), as.matrix(g_out))
g_out_trans = g_out[id %in% combined[cis_flag == FALSE]$snp]
fwrite(g_out, "Output/matrixEQTL/SNP.txt", sep = "\t")
fwrite(g_out_trans, "Output/matrixEQTL/SNP_trans.txt", sep = "\t")

g_out = t(gene_expression)
g_out = data.table(id = row.names(g_out), as.matrix(g_out))
g_out_trans = g_out[id %in% combined[cis_flag == FALSE]$pid]
fwrite(g_out, "Output/matrixEQTL/GE.txt", sep = "\t")
fwrite(g_out_trans, "Output/matrixEQTL/GE_trans.txt", sep = "\t")

g_out = unique(data.table(snp = combined$snp, chr = "chr19", pos = combined$sid_pos))
g_out_trans = g_out[snp %in% combined[cis_flag == FALSE]$snp]
fwrite(g_out, "Output/matrixEQTL/snpsloc.txt", sep = "\t")
fwrite(g_out_trans, "Output/matrixEQTL/snpsloc_trans.txt", sep = "\t")

g_out = unique(data.table(geneid = combined$pid, chr = "chr19", s1 = combined$gene_start, s2 = combined$gene_end))
g_out_trans = g_out[geneid %in% combined[cis_flag == FALSE]$pid]
fwrite(g_out, "Output/matrixEQTL/geneloc.txt", sep = "\t")
fwrite(g_out_trans, "Output/matrixEQTL/geneloc_trans.txt", sep = "\t")

rm(all_assoc_genes, combined, test_genotype)
gc()
