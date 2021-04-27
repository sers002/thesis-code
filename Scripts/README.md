# Scripts

This directory contains all of the scripts (*R* and *bash*) required to run the framework. The order in which the script files should be run is detailed below (some scripts rely on data generated in previous scripts). Note the documentation in the `Output/` directory details which script file each output file is generated within for additional information. 

## Script Descriptions

### eQTL Data Cleaning and Causal Variant Generation (`eqtls/`)

#### Quality Control (`quality_control/`)

##### **`0._generate_gtex_snp_list_all_vars.R`**

This file creates a list of all SNP IDs that are contained within the eQTL file
(`Data/lung_eqtl_data/eqtls.txt`), and outputs to a text file in 
`Output/quality_control/gtex_all_snps.txt`. This list of SNPs is then used in 
subsequent steps, such as HWE, AF checking, etc.
  
##### **`1._hardy_weinberg.sh`**

This script runs an exact test for each SNP using the 1000G data and the list 
of SNPs created in the previous step. It outputs a file with the statistics 
from each test in `Output/quality_control/plink.hwe.gz`. 

##### **`2._create_list_of_hwe_fails.R`**

This script reads in the results of `1._hardy_weinberg.sh`, adjusts the 
$p$-values using a Bonferroni adjustment, and then outputs a list of SNP IDs 
where the adjusted $p$-value is below 0.05 to the 
`Output/quality_control/bad_hardy_snps.txt` file. 

##### **`3._create_list_of_big_difference_mafs.R`**

This script reads in both the 1000 Genomes and CoDeS3D dataset to compare their 
allele frequencies. The allele frequencies are calculated in the case of the 
1000 Genomes dataset as it has been subset to only European participants. 
Variant IDs where the difference in allele frequency is above 10\% are output 
to `Output/quality_control/out_of_range_snps.txt`. This script additionally 
creates a plot which compares the allele frequencies 
(`Plots/quality_control/allele_freq_gtex_vs_1000g.png`).

#### Data Preparation (`dataprep_and_exploratory/`)

##### **`1._process_eqtl_data.R`**

This script reads in the eQTL dataset (CoDeS3D data: 
`Data/lung_eqtl_data/eqtls.txt`, `genes.txt`) and genotype dataset (1000 
Genomes European subset: 
`Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz`). It subsets the 
eQTL dataset to only those that are also contained in the genotype dataset and 
removes variants that were identified as low quality in the quality control 
scripts. Finally, it calculates some summary information about each 
variant/gene (e.g. narrow-sense heritability, SNP density, etc). The results 
are **not** saved to any file; the resulting dataset is stored in the *R* 
environment.

##### **`2._checking_plink_files_and_flipping.R`**

This script runs the preceding script file (`1._process_eqtl_data.R`) and 
checks whether both the 1000 Genomes European subset and *HAPGEN2*-simulated 
datasets have the same reference allele as the GTEx/CoDeS3D data.

##### **`3._LD_declumping.R`**

This script runs the eQTL processing script (`1._process_eqtl_data.R`) and 
reads in the linkage disequilibrium and pruned variant list. It removes SNPs 
that were not retained by the pruning routine and saves the dataset in 
`Output/eqtls_for_sojo/declumped_eqtls.csv`. 

##### **`4._LD_matrix.R`**

This script imports the linkage disequilibrium dataset, subsets it to only 
those SNPs in the pruned dataset, and saves a sparse matrix representation 
of it (in binary RDS format) for use in the *SOJO* routine. This is required 
due to memory constraints when runing the *SOJO* loop (otherwise it resulted 
in an error, even if the raw linkage disequilibrium dataset is removed from 
the *R* environment after the sparse matrix is constructed). 

#### Determining Plausible Set of Causal Effects (`sojo/`)

##### **`5._sojo_pipeline.R`**

This script iterates over a series of proportions, running *SOJO* for each 
gene and storing the results of the highest set of coefficients where the set 
of conditions apply ($h^2 < 0.8$, no sign change, not larger than input 
coefficients). The results for each proportion are saved as 
`Output/sojo_output/sojo_all_nvar_passes.RData`.

##### **`6._setting_causal_effects.R`**

This script reads in the datasets created in the preceding script 
(`5._sojo_pipeline.R`) and sets the causal effects for each gene to the set of 
coefficients from the highest proportion possible for that gene. The resulting 
set of coefficients for all genes are saved in a dataset, 
`Output/sojo_output/set_betas.RData`.

##### **`7._total_association_when_include_LD.R`**

This script takes the dataset created in the preceding script 
(`6._setting_causal_effects.R`) and determines the marginal effects/total associations 
based on the causal variants ($\beta_i = \sum_j \gamma_j r_{ij}$). The results 
are saved in `Output/sojo_output/all_betas.RData`.

##### **`8._expression_by_matmult_from_SOJO.R`**

This script takes the causal effects and generates gene expression values for 
each gene for individuals in the 1000 Genomes European subset. This code was 
generalised to run over multiple datasets for the empirical simulations, and 
remains useful for the subsequent script file.

##### **`9._observed_var.R`**

This script takes the gene expression calculated in the preceding script and 
calculates the observed variance for each gene's expression. These observed 
variances are used later to verify the assumptions of the semi-analytic method 
and to generate estimates. 


### eQTL Study Power Estimation Methods (`full_simulations/`)

##### **`plotting_power_results.R`**

This script imports the power results contained in the `Output/power` directory 
and creates plots of the comparisons of interest within the thesis. 

##### **`9._combinations_counts.R`**

This script determines the total number of associations (cis and trans) based 
on the number of associations present in the CoDeS3D dataset and the number of 
variants in the 1000 Genomes dataset.


#### Empirical Method (`empirical/`)

##### **`10._loop_over_hapgen_sims.R`**

This script iterates over the simulated genotype datasets and runs 
`analyse_one_hapgen_dataset.R` using *callr* on each, storing the results for 
each simulated genotype dataset in an individual file.

##### **`11._loop_over_matrixeqtl_results.R`**

This script iterates over the *MatrixEQTL* results generated by the preceding 
script (`10._loop_over_hapgen_sims.R`) and estimates the empirical average 
power of detection of eQTLs based on them.

##### **`11.5._loop_over_matrixeqtl_results_egenes.R`**

This script iterates over the *MatrixEQTL* results generated by the preceding 
script (`10._loop_over_hapgen_sims.R`) and estimates the empirical average 
power of detection of eGenes based on them.

##### **`analyse_one_hapgen_dataset.R`**

This script is run within the `10._loop_over_hapgen_sims.R` script using the 
*callr* pacakge. This script takes a filename as a command-line argument, 
imports that file and runs the required data manipulation/analysis. 

#### Analytic Method (`2_analytic/`)

##### **`analytic_power_bonferroni.R`**

This script takes a set of gene-variant associations and calculates the 
estimated statistical power at a set of sample sizes using a Bonferroni 
adjustment with a user-defined number of associations. The estimates are 
calculated using statistical theory, with calculation of the non-centrality 
parameter based on the sample size, effect size, and minor allele frequency. 

##### **`analytic_power_fdr.R`**

This script takes a set of gene-variant associations and calculates the 
estimated statistical power at a set of sample sizes controlling the FDR using 
the Benjamini-Hochberg method. The estimates are calculated using statistical 
theory where the $p$-value threshold corresponding to the desired FDR is 
determined by optimisation with the `uniroot` function. Determining the 
$p$-value threshold is necessary due to its dependence on the $p$-value 
distribution.

##### **`analytic_power_alternate.R`**

This script uses the same methodology as `analytic_power_fdr.R` to calculate 
the power of detecting eQTLs under a set of assumptions about the omitted 
variants (i.e. the ones not in physical contact).

#### Semi-Analytic Method (`3_semi_analytic/`)

##### **`sa_multiple_sample_sizes_mvn.R`**

This script takes associations and linkage disequilibrium datasets and 
simulates a number of realisations of the associations based on the correlation 
between variants. It achieves this by iterating over the genes, and generating 
coefficient estimates for each gene separately. The resulting dataset is saved 
to a file for subsequent analysis. 

##### **`analysing_sa_results.R`**

This script takes the results of the preceding script and calculates the 
estimated statistical power to detect eQTLs based on the realisations. 

##### **`analysing_sa_results_egenes.R`**

This script takes the results of the preceding script and calculates the 
estimated statistical power to detect eGenes based on the realisations. 


### Genotype Simulation (`genotype_simulation/`)

##### **`generate_id_file.R`**

This script outputs a list of sample IDs from the 1000 Genomes dataset based 
on populations (e.g. a list of the sample IDs of all participants from the 
YRI population).

##### **`hapgen2_full_process_repeat.sh`**

This script performs some preliminary data cleaning and manipulation on the 
1000 Genomes dataset (subsetting the dataset to only EUR participants, removing 
variants with MAF \< 1% or missing data, biallelic SNPs). It also converts the 
genotype dataset from VCF format to the required Oxford format using *PLINK2*. 
It then iterates over sample sizes, generating 100 datasets of each sample size 
and storing them in `Output/genotype_simulation/rep_sims/n_#/` where `#` is the 
sample size.

##### **`pca_on_sims_multiple_sims.R`**

This script reads in the original 1000 Genomes dataset and two simulated 
genotype datasets of size $N = 404$. After the required data manipulation, the 
script outputs the PCA plot (`Plots/genotype_simulation/pca_comparison.pdf`) 
and ANOVA results. 


### Linkage Disequilibrium (`linkage_disequilibrium/`)

##### **`generate_gtex_snp_list.R`**

This script generates a list of variant IDs from the CoDeS3D dataset and 
outputs them to a file that is used for calculating linkage disequilibrium 
and variant pruning. 

##### **`plink_ld_process.sh`**

This script uses *PLINK* to estimate the linkage disequilibrium between each 
pair of variants contained in the CoDeS3D dataset based on the 1000 Genomes 
European subset. This script also contains similar calculations for the Yoruba 
ancestry group. Finally, this script contains the *PLINK* command used to 
create the extended linkage disequilibrium dataset used to demonstrate the 
decay of linkage disequilibrium according to distance.  

### Variant Pruning (`variant_pruning/`)

##### **`generate_gtex_snp_list2.R`**

This script generates a list of gene-variant pairs from the CoDeS3D dataset 
and outputs them to `Output/variant_pruning/snp_list_from_gtex2.txt` which is 
used in the subsequent variant pruning script. 

##### **`plink_pruning_gene.sh`**

This script iterates over the gene-variant pairs generated in the preceding 
script file. It iteratres over the unique gene IDs contained in the file, and 
uses the `grep` program to extract out the variants that have a pairing with 
each gene. The variant pruning routine from *PLINK* is run using only those 
variants. As *PLINK* only saves the list of SNP IDs, `awk` is used to append 
the gene ID to each line and the results are saved into a file, 
`Output/variant_pruning/plink_gene_prune.in`.

## *R* Packages Required

* `data.table`
* `Matrix` 
* `Matrix.utils` 
* `sojo`
* `ggplot2`
* `callr`
* `magrittr`
* `here`
* `MatrixEQTL`
* `mvtnorm`
* `qvalue`

To install all required *R* packages, the following code can be used.

```r
required_packages = c(
  "data.table",
  "Matrix.utils",
  "sojo",
  "ggplot2",
  "callr",
  "magrittr",
  "here",
  "MatrixEQTL",
  "mvtnorm"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

## The qvalue package is on Bioconductor, not CRAN.
if (!requireNamespace("qvalue", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  BiocManager::install("qvalue")
}
```

## Script Order

Not all scripts are included in this list as other scripts that depend on them 
may automatically run them (e.g. `hapgen2_full_process_repeat.sh` runs the 
`generate_id_file.R` script automatically). 

1. `genotype_simulation/hapgen2_full_process_repeat.sh`
2. `quality_control/0._generate_gtex_snp_list_all_vars.R`
3. `quality_control/1._hardy_weinberg.sh`
4. `quality_control/2._create_list_of_hwe_fails.R`
5. `quality_control/3._create_list_of_big_difference_mafs.R`
6. `linkage_disequilibrium/plink_ld_process.sh`
7. `variant_pruning/plink_pruning_gene.sh`
8. `eqtls/dataprep_and_exploratory/2._checking_plink_files_and_flipping.R`
9. `eqtls/dataprep_and_exploratory/3._LD_declumping.R`
10. `eqtls/dataprep_and_exploratory/4._LD_matrix.R`
11. `eqtls/sojo/5._sojo_pipeline.R`
12. `eqtls/sojo/6._setting_causal_effects.R`
13. `eqtls/sojo/7._total_association_when_include_LD.R`
14. `eqtls/sojo/8._expression_by_matmult_from_SOJO.R`
15. `eqtls/sojo/9._observed_var.R`
16. `eqtls/full_simulations/9._combinations_counts.R`
17. `eqtls/full_simulations/empirical/10._loop_over_hapgen_sims.R`
18. `eqtls/full_simulations/empirical/11._loop_over_matrixeqtl_results.R`
19. `eqtls/full_simulations/empirical/11.5._loop_over_matrixeqtl_results_egenes.R`
20. `eqtls/full_simulations/semi_analytic/sa_multiple_sample_sizes_mvn.R`
21. `eqtls/full_simulations/semi_analytic/analysing_sa_results.R`
22. `eqtls/full_simulations/semi_analytic/analysing_sa_results_egenes.R`
23. `eqtls/full_simulations/analytic/analytic_power_fdr.R`
24. `eqtls/full_simulations/plotting_power_results.R`
25. `eqtls/full_simulations/plotting_power_results_egenes.R`

