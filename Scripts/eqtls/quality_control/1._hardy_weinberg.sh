## Get SNP IDs of duplicate variants
##  from: https://stackoverflow.com/questions/9863208/how-do-i-remove-duplicated-snps-using-plink

gunzip -c  ./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz |
  grep "^[^##]" |
  cut -f3 |
  sort |
  uniq -d > ./Output/linkage_disequilibrium/ld_ref_duplicate_vars.dupvar

./Programs/plink \
  --vcf ./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz \
  --exclude ./Output/linkage_disequilibrium/ld_ref_duplicate_vars.dupvar \
  --extract ./Output/quality_control/gtex_all_snps.txt \
  --hardy gz \
  --keep-allele-order

mv plink.hwe.gz Output/quality_control/plink.hwe.gz
rm plink.log
rm plink.nosex
