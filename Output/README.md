# Output

This directory contains output from the various steps in the framework. This documentation describes each file and which script file that generates the file.

## Output Traceback

### Quality Control (`quality_control/`)

`bad_hardy_snps.txt`
: Source: *`eqtls/quality_control/2._create_list_of_hwe_fails.R`*
: A list of SNP IDs where the HWE exact test $p$-value was lower than the Bonferroni threshold.

`gtex_all_snps.txt`
: Source: *`eqtls/quality_control/0._generate_gtex_snp_list_all_vars.R`*
: A list of SNP IDs contained in the GTEx eQTL dataset.

`out_of_range_snps.txt`
: Source: *`eqtls/quality_control/3._create_list_of_big_difference_mafs.R`*
: A list of SNP IDs where the GTEx and 1000G allele frequencies differ by more than 10%.

`plink.hwe.gz`
: Source: *`eqtls/quality_control/1._hardy_weinberg.sh`*
: A dataset of HWE exact tests with p-values.

### Linkage Disequilibrium (`linkage_disequilibrium/`)

`ld_gtex_snps.ld.gz`
: Source: *`eqtls/quality_control/1._hardy_weinberg.sh`*
: A dataset containing pairwise LD values (correlations, $R$) for all SNPs within 250kb.

`ld_gtex_snps_decay.ld.gz`
: Source: *`eqtls/quality_control/1._hardy_weinberg.sh`*
: A dataset containing pairwise LD values (correlations, $R$) for all SNPs within 750kb. For use in plotting LD.

`ld_mat.rds`, `snp_ref2.rds`
: Source: *`eqtls/quality_control/1._hardy_weinberg.sh`*
: R datasets required for SOJO algorithm. `ld_mat` is a sparse matrix created from `ld_gtex_snps.ld.gz` and subset to only pruned variants. `snp_ref2` is a vector indicating the reference allele for each SNP.

`ld_ref_duplicate_vars.dupvar`
: Source: *`eqtls/quality_control/1._hardy_weinberg.sh`*
: A list of SNP IDs that are duplicated in the 1000G dataset.

`snp_list_from_gtex.txt`
: Source: *`eqtls/quality_control/1._hardy_weinberg.sh`*
: A list of SNP IDs that are contained in the CoDeS3D dataset.

### Variant Pruning (`variant_pruning/`)

`plink_gene_prune.in`
: Source: *`variant_pruning/plink_pruning_gene.sh`*
: A list of gene-variant associations to retain after variant pruning by gene.

`1000g_chr19_eur_subset_gtex.vcf.gz`
: Source: *`variant_pruning/plink_pruning_gene.sh`*
: A subset of the 1000 Genomes European subset with only variants contained in the CoDeS3D/GTEx dataset.

`snp_list_from_gtex2.txt`
: Source: *`variant_pruning/generate_gtex_snp_list2.R`*
: A list of gene-variant associations that are contained in the CoDeS3D dataset.

### Genotype Simulation (`genotype_simulation/`)

`1000g_chr19_eur_subset.vcf.gz`
: Source: *`genotype_simulation/hapgen2_full_process_repeat.sh`*
: A subset of the 1000 Genomes chromosome 19 dataset with only non-Finnish European participants.

`1000g_chr19_eur_subset_for_hapgen.haps`
: Source: *`genotype_simulation/hapgen2_full_process_repeat.sh`*
: The same dataset as `1000g_chr19_eur_subset.vcf.gz` but converted to the Oxford genomic format required by *HAPGEN2*.

`rep_sims/n_{n}/1000g_sim_n{n}_{i}.controls.haps.gz`
: Source: *`genotype_simulation/hapgen2_full_process_repeat.sh`*
: A set of simulated genotype datasets for sample size `{n}` and index `{i}`.

### MatrixEQTL (`matrixEQTL/`)

`geneloc.txt`, `snpsloc.txt`
: Source: *`eqtls/sojo/8._expression_by_matmult_from_SOJO.R`*
: Gene and SNP positions on the genome. 

`GE.txt`
: Source: *`eqtls/sojo/8._expression_by_matmult_from_SOJO.R`*
: Gene expression values for each individual and gene.

`SNP.txt`
: Source: *`eqtls/sojo/8._expression_by_matmult_from_SOJO.R`*
: Genotype values for each individual and SNP.

`rep_sims/n_#/output_cis_#.txt.gz`, `rep_sims/n_#/output_trans_#.txt.gz`
: Source: *`full_simulations/1_matrixEQTL/10._loop_over_hapgen_sims.R`*
: Summary statistics from marginal regressions between each variant and gene used in the analysis. Each file is for a given simulation realisation. 

### Power Results (`power/`)

`matrixEQTL.csv`
: Source: *``*
: The estimated power to detect eQTLs using the empirical method.

`sa_mvn.csv`
: Source: *``*
: The estimated power to detect eQTLs using the semi-analytic method.

`sa_mvn_extended.csv`
: Source: *``*
: The estimated power to detect eQTLs for large sample sizes using the semi-analytic method.

### Semi-Analytic Results (`semi_analytic/`)

`sa_results_mvn.rds`
: Source: *`full_simulations/3_semi_analytic/sa_multiple_sample_sizes_mvn.R`*
: A set of realisations from the semi-analytic method. Each row is a gene-variant association and contains 100 realisations of marginal effect estimates.

`sa_results_mvn_extended.rds`
: Source: *`full_simulations/3_semi_analytic/sa_multiple_sample_sizes_mvn.R`*
: A set of realisations from the semi-analytic method. Each row is a gene-variant association and contains 100 realisations of marginal effect estimates for larger sample sizes. 

### SOJO Results (`sojo_output/`)

`all_betas.RData`
: Source: *`eqtls/sojo/7._total_association_when_include_LD.R`*
: A dataset of all variants containing the causal effects and total association values as given by SOJO. 
