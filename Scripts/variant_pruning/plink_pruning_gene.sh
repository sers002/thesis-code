
## Generate list of SNP IDs that are present within the GTEx file
Rscript ./Scripts/linkage_disequilibrium/generate_gtex_snp_list.R

## Generate list of SNP IDs for each gene
Rscript ./Scripts/variant_pruning/generate_gtex_snp_list2.R

## Get SNP IDs of duplicate variants
##  from: https://stackoverflow.com/questions/9863208/how-do-i-remove-duplicated-snps-using-plink
gunzip -c  ./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz |
  grep "^[^##]" |
  cut -f3 |
  sort |
  uniq -d > ./Output/linkage_disequilibrium/ld_ref_duplicate_vars.dupvar

./Programs/bcftools view -i'ID=@Output/linkage_disequilibrium/snp_list_from_gtex.txt' ./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz |
  ./Programs/bgzip > ./Output/variant_pruning/1000g_chr19_eur_subset_gtex.vcf.gz

## Remove file if it already exists
rm ./Output/variant_pruning/plink_gene_prune.in

touch ./Output/variant_pruning/plink_gene_prune.in

## Extract a list of unique genes using the following method:
##   Take the 1st column (-f1) of snp list file that's delimited with spaces (-d' ')
##   Remove duplicates (uniq)
## Then loop over this list of unique genes
for gene in $(cut -d' ' -f1 ./Output/variant_pruning/snp_list_from_gtex2.txt | uniq)
do
  echo $gene

  ## Extract gene-snp associations from snp list file for current gene (fgrep),
  ## Take the 2nd column (which is SNP rsids) and output them to the curr_gene_snps.txt file
  fgrep "$gene" ./Output/variant_pruning/snp_list_from_gtex2.txt | cut -d' ' -f2- > ./Output/variant_pruning/curr_gene_snps.txt

  ## Run plink variant pruning using only SNPs in curr_gene_snps.txt (made in previous step)
  ./Programs/plink \
    --vcf ./Output/variant_pruning/1000g_chr19_eur_subset_gtex.vcf.gz \
    --keep-allele-order \
    --extract ./Output/variant_pruning/curr_gene_snps.txt \
    --indep-pairwise 50 5 0.81

  ## Append the gene name to each line of plink.prune.in (which is a list of rsids)
  ## and append onto the plink_gene_prune.in file
  awk -v gene="$gene" '{print gene, $0}' plink.prune.in >> ./Output/variant_pruning/plink_gene_prune.in
done
