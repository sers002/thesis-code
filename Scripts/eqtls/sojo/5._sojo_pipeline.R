library(sojo)
library(data.table)
library(Matrix.utils)
library(Matrix)

## Generating the LD matrix takes too much memory if it's calculated at the time,
## so the results are saved as RDS files
# source("Scripts/eqtls/dataprep_and_exploratory/4._LD_matrix.R")

## These are the results from 4._LD_matrix.R (as of Nov 14)
file_match = fread("Output/eqtls_for_sojo/declumped_eqtls.csv")
ld_mat     = readRDS("Output/linkage_disequilibrium/ld_mat.rds")
snp_ref2   = readRDS("Output/linkage_disequilibrium/snp_ref2.rds")

num_to_try = c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

file_list = vector("list", length = length(num_to_try))
names(file_list) = paste0("output_", num_to_try)

gene_list = unique(file_match$pid)

for(j in num_to_try) {
  message(j)
  message("  ", length(gene_list))
  
  full_out = NULL

  pb <- txtProgressBar(max = length(gene_list), style = 3)

  for(i in gene_list) {
    setTxtProgressBar(pb, which(gene_list == i))
    
    ## Subset to associations to the one gene
    ss <- file_match[pid == i]
    sum.stat.discovery  = ss[, c("snp", "a1", "a2", "maf", "b", "b_se")]
    setnames(sum.stat.discovery, "snp", "SNP")
    setnames(sum.stat.discovery, "a1", "A1")
    setnames(sum.stat.discovery, "a2", "A2")
    setnames(sum.stat.discovery, "maf", "Freq1")
    setnames(sum.stat.discovery, "b_se", "se")
    
    ## N = 515 as that's the sample size of GTEx
    sum.stat.discovery[, N := 515]
  
    ## Gene-specific LD matrix and reference vector
    ld_mat2 <- as.matrix(ld_mat[ss$snp, ss$snp])
    snp_ref3 <- snp_ref2[names(snp_ref2) %in% colnames(ld_mat2)]
    
    ## Run SOJO with gene's data and current proportion allowed
    sumstats_out = tryCatch({
      res2 <- sojo(
        sum.stat.discovery, 
        LD_ref = ld_mat2, 
        snp_ref = snp_ref3, 
        nvar = max(round(nrow(sum.stat.discovery) * j, 0), 1),   
        standardize = TRUE
      )
      
      b <- res2$beta.mat
      
      ## Calculate percent change from original coefficient from GTEx
      size_check1 = apply(b, MARGIN = 2, function(x) {
        ifelse(round(x, 5) == 0 | round(ss$b, 5) == 0, 0, x/ss$b - 1)
      })
      
      ## Check the coefficient hasn't changed direction or increased in size
      size_check = apply(size_check1, MARGIN = 2, function(x) {
        all(x>=-1 & x <= 0)
      })
      
      ## Calculate the total heritability for each set of coefficients
      q2_gamma_vec = apply(b, MARGIN = 2, function(x) {
        sum((x^2)*2*ss$maf*(1-ss$maf))
      })
      
      ## Choose the column closest to 80% heritability
      h2_max = min(ss$h2_pruned[1], 0.8)
      chosen.one = which.min(abs(h2_max - (q2_gamma_vec + (!size_check * Inf))))
      
      if (q2_gamma_vec[chosen.one] > h2_max) {
        chosen.one  = chosen.one - 1
      }
      
      ## Store chosen betas 
      betas = b[, chosen.one]
      out = data.table(snp = names(betas), beta_shrink = betas)
      rm(betas)
      out
      }, 
      error = function(e) {
        ## If SOJO resulted in an error, recover from it
        data.table(snp = sum.stat.discovery$SNP, beta_shrink = NA)
      })
    
    sumstats_out = merge(
      sumstats_out, 
      ss, 
      by.x = "snp", 
      by.y = "snp",
      all.x = TRUE
    )
    
    full_out = rbind(full_out, sumstats_out)
  }
  
  ## Remove genes where SOJO failed
  full_out = full_out[!is.na(beta_shrink)]
  file_list[[paste0("output_", j)]] = copy(full_out)
  gene_list = unique(full_out$pid)
  
  gc()
}


output_0.05 = file_list$output_0.05
output_0    = file_match[!pid %in% output_0.05$pid]
output_0[, beta_shrink := NA]
output_0.1  = file_list$output_0.1
output_0.15 = file_list$output_0.15
output_0.2  = file_list$output_0.2
output_0.3  = file_list$output_0.3
output_0.4  = file_list$output_0.4
output_0.5  = file_list$output_0.5
output_0.6  = file_list$output_0.6
output_0.7  = file_list$output_0.7
output_0.8  = file_list$output_0.8
output_0.9  = file_list$output_0.9
output_1    = file_list$output_1

rm(file_list)

save(
  list = c(
    "output_0",
    "output_0.05",
    "output_0.1", 
    "output_0.15",
    "output_0.2", 
    "output_0.3",
    "output_0.4",
    "output_0.5",
    "output_0.6", 
    "output_0.7", 
    "output_0.8", 
    "output_0.9",
    "output_1"
  ),
  file = "Output/sojo_output/sojo_all_nvar_passes.RData"
)
