library(data.table)

load(file = here::here("Output/sojo_output/all_betas.RData"))


# Calculating 1000G MAFs --------------------------------------------------

## Import 1000G dataset (to calculate the observed MAFs)
og.geno = fread(
  here::here("Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz"),
  skip = 100
)

setnames(
  og.geno,
  c("ID",  "ALT", "REF"),
  c("snp", "a1", "a2")
)

og.geno[, 10:ncol(og.geno)] <- lapply(
  og.geno[, 10:ncol(og.geno)],
  function(x) {
    as.numeric(substr(x,1,1)) + as.numeric(substr(x,3,3))
  }
)

## Calculate MAFs for each SNP in the 1000G
og.maf <- rowMeans(og.geno[, 10:ncol(og.geno)]) / 2
og.varx <- apply(og.geno[, 10:ncol(og.geno)], MARGIN = 1, FUN = var)
og.maf.dt <- data.table(
  snp = og.geno$snp,
  maf_1000g = og.maf,
  a2_1000g = og.geno$a2,
  a1_1000g = og.geno$a1,
  varx_1000g = og.varx
)
og.maf.dt <- og.maf.dt[!duplicated(og.maf.dt$snp)]

og.maf.dt <- og.maf.dt[snp %in% combined$snp]

combined <- merge(combined, og.maf.dt[, c("snp", "maf_1000g", "varx_1000g")], by = "snp")


# Calculating Var(Y) ------------------------------------------------------

## Calculate (observed) Var(Y) - this is based on N=404 (1000G data)
gene_expression <- fread(here::here("Output/matrixEQTL/GE.txt"))
expression_vars <- data.table(
  pid = gene_expression$id,
  var_y = apply(gene_expression[, -1], MARGIN = 1, var)
)

combined <- merge(combined, expression_vars, by = "pid")
combined[, sigma_sq := var_y - (alpha^2 * 2 * maf_1000g * (1 - maf_1000g))]

save(combined, file = "Output/sojo_output/all_betas.RData")
