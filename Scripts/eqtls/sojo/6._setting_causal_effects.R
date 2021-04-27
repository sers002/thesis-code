library(data.table)

## Read in the results from SOJO
load(file = "Output/sojo_output/sojo_all_nvar_passes.RData")

seed = 3475387

## Take the largest proportion that worked for each gene
combined = copy(output_1)
combined = rbind(combined, output_0.9[!pid %in% combined$pid])
combined = rbind(combined, output_0.8[!pid %in% combined$pid])
combined = rbind(combined, output_0.7[!pid %in% combined$pid])
combined = rbind(combined, output_0.6[!pid %in% combined$pid])
combined = rbind(combined, output_0.5[!pid %in% combined$pid])
combined = rbind(combined, output_0.4[!pid %in% combined$pid])
combined = rbind(combined, output_0.3[!pid %in% combined$pid])
combined = rbind(combined, output_0.2[!pid %in% combined$pid])
combined = rbind(combined, output_0.1[!pid %in% combined$pid])
combined = rbind(combined, output_0.15[!pid %in% combined$pid])
combined = rbind(combined, output_0.05[!pid %in% combined$pid])
combined = rbind(combined, output_0[!pid %in% combined$pid])

## Calculate total heritability based on SOJO causal effects
combined[, q2 := (beta_shrink^2) * 2 * maf * (1 - maf)]
combined[, h2 := sum(q2, na.rm = TRUE), by = pid]

combined[is.na(beta_shrink), beta_shrink := 0]

## Checks of original and SOJO causal effects
# summary(abs(combined$b))
# plot(density(abs(combined$b)))
# summary(abs(combined[beta_shrink != 0]$beta_shrink))
# plot(density(abs(combined[beta_shrink != 0]$beta_shrink)))
# min(abs(combined[beta_shrink != 0]$beta_shrink))

## Round beta shrink to 7 dp
combined[, beta_shrink_raw := beta_shrink]
combined[, beta_shrink := round(beta_shrink, 7)]
combined[, picked_snps_per_gene := sum(beta_shrink != 0), by = pid]

# summary(combined$picked_snps_per_gene)

## Checks for cis causal effects
# nrow(combined[cis_flag == FALSE & beta_shrink != 0])
# nrow(combined[ beta_shrink != 0])
# length(unique(combined[cis_flag == FALSE & beta_shrink != 0]$pid))
# summary(abs(combined[cis_flag == FALSE ]$beta_shrink))
# summary(abs(combined[cis_flag == TRUE ]$beta_shrink))
# summary(abs(combined[cis_flag == FALSE ]$b))
# summary(abs(combined[cis_flag == TRUE ]$b))
# plot(density(abs(combined$b)))
# plot(density(abs(combined$beta_shrink)))
# plot(density(abs(combined[beta_shrink != 0]$beta_shrink)))

## Remove everything from the environment except the coefficient dataset
rm(list = setdiff(ls(), "combined"))

save(combined, file = "Output/sojo_output/set_betas.RData")
