
## Generate list of SNP IDs that are present within the GTEx file (from snps.txt)

Rscript ./Scripts/linkage_disequilibrium/generate_gtex_snp_list.R

## Get SNP IDs of duplicate variants
##  from: https://stackoverflow.com/questions/9863208/how-do-i-remove-duplicated-snps-using-plink

gunzip -c  ./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz |
  grep "^[^##]" |
  cut -f3 |
  sort |
  uniq -d > ./Output/linkage_disequilibrium/ld_ref_duplicate_vars.dupvar

## Change R^2 and kb cutoff - this is LD now in use (Nov 13)

./Programs/plink \
  --vcf ./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz \
  --exclude ./Output/linkage_disequilibrium/ld_ref_duplicate_vars.dupvar \
  --extract ./Output/linkage_disequilibrium/snp_list_from_gtex.txt \
  --keep-allele-order \
  --r gz \
  --ld-snp-list ./Output/linkage_disequilibrium/snp_list_from_gtex.txt \
  --ld-window-kb 250 \
  --ld-window 50000 \
  --ld-window-r2 0.0 \
  --out ./Output/linkage_disequilibrium/ld_gtex_snps

## AFR LD
./Programs/plink \
  --vcf ./Output/genotype_simulation/1000g_chr19_yri_subset.vcf.gz \
  --exclude ./Output/linkage_disequilibrium/ld_ref_duplicate_vars.dupvar \
  --extract ./Output/linkage_disequilibrium/snp_list_from_gtex.txt \
  --keep-allele-order \
  --r gz \
  --ld-snp-list ./Output/linkage_disequilibrium/snp_list_from_gtex.txt \
  --ld-window-kb 250 \
  --ld-window 50000 \
  --ld-window-r2 0.0 \
  --out ./Output/linkage_disequilibrium/ld_gtex_snps_yri


## LD for LD decay plot

./Programs/plink \
  --vcf ./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz \
  --r gz \
  --keep-allele-order \
  --exclude ./Output/linkage_disequilibrium/ld_ref_duplicate_vars.dupvar \
  --extract ./Output/linkage_disequilibrium/snp_list_from_gtex.txt \
  --ld-snp-list ./Output/linkage_disequilibrium/snp_list_from_gtex.txt \
  --ld-window-kb 700 \
  --ld-window 50000 \
  --ld-window-r2 0 \
  --out ./Output/linkage_disequilibrium/ld_gtex_snps_decay
