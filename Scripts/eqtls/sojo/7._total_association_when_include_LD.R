library(data.table)

## SOJO causal effects
load(file = "Output/sojo_output/set_betas.RData")

## Set one causal variant per gene
# set.seed(12345)
# 
# for(g in unique(combined$pid)) {
#   non_zero <- combined[pid == g & beta_shrink != 0, ]$snp
#   
#   if (length(non_zero) > 1) {
#     chosen_one <- sample(non_zero, size = 1)
#     combined[pid == g & !(snp %in% chosen_one), beta_shrink := 0]
#   }
# }

source("Scripts/eqtls/dataprep_and_exploratory/1._process_eqtl_data.R")

ld.df <- fread("Output/linkage_disequilibrium/ld_gtex_snps.ld.gz")
ld.df[, ID := paste(pmin(SNP_A, SNP_B), pmax(SNP_A, SNP_B))]
ld.df = ld.df[!duplicated(ID)]
ld.df[, ID := NULL]

gc()

## Combine the coefficients and the original dataset - anything that was removed 
## via declumping will have causal effect of zero
filematch_sub = file_match[!paste0(snp, pid) %in% paste0(combined$snp, combined$pid)]
combined = rbindlist(list(combined, filematch_sub), fill = TRUE)
combined[is.na(beta_shrink), beta_shrink := 0]

rm(file_match, filematch_sub)

genes_list = unique(combined$pid)

full_out <- NULL

pb <- txtProgressBar(max = length(genes_list), style = 3)

## Loop over genes, calculating the total association for each variant
for(i.ind in seq_along(genes_list)) {
  setTxtProgressBar(pb, i.ind)
  
  i <- genes_list[i.ind]
  test = combined[pid == i]
  
  ld_sub = ld.df[SNP_A %in% test$snp & SNP_B %in% test$snp]

  ld_sub2 = ld_sub[, c("SNP_A","SNP_B", "R")]
  ld_sub3 = ld_sub[, c("SNP_B","SNP_A", "R")]
  
  setnames(ld_sub2, "SNP_A", "SNP")
  setnames(ld_sub3, "SNP_B", "SNP")
  setnames(ld_sub2, "SNP_B", "SNP_other")
  setnames(ld_sub3, "SNP_A", "SNP_other")
  
  ld_sub = rbind(ld_sub2, ld_sub3)
  ld_sub = unique(ld_sub)
  
  rm(ld_sub2, ld_sub3)
  
  ## Merge the dataset as a cartesian join
  ld_sub = merge(
    ld_sub, 
    test[, c("pid", "b", "maf", "sid_pos", "beta_shrink", "snp")], 
    by.x = "SNP_other", 
    by.y = "snp", 
    all.x = TRUE
  )
  
  ## Total association based on LD
  ld_sub[, br := beta_shrink * R]
  ld_sub[, alpha := sum(br), by = SNP]
  ld_sub[, r_count := length(R), by = SNP]
  ld_sub[, r_total := sum(abs(R)), by = SNP]
  ld_sub[, R2 := R^2]
  
  ## Anything with R^2 over 0.8 with a causal SNP is considered an eQTL
  ld_sub[, qtl_tag := any(R2 >= 0.8 & SNP != SNP_other & beta_shrink != 0), by = SNP]
  
  assoc = unique(ld_sub[, c("SNP", "alpha", "r_count", "r_total", "qtl_tag")])

  assoc = merge(assoc, test, by.x = "SNP", by.y = "snp", all.y = TRUE)
  assoc[is.na(alpha), alpha := 0]
  assoc[is.na(qtl_tag), qtl_tag := FALSE]
  
  message("test: ", nrow(test), " ; assoc: ", nrow(assoc))

  full_out <- rbind(full_out, assoc)
  
  gc()
}

close(pb)

combined = copy(full_out)

setnames(combined, "SNP", "snp")

combined[beta_shrink != 0, qtl_tag := TRUE]

save(combined, file = "Output/sojo_output/all_betas.RData")
