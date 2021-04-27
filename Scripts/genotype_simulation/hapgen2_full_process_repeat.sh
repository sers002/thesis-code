################################################################################
### Preliminaries                                                            ###
################################################################################

## Generate a text file of EUR participant IDs from the 1000 Genomes metadata

Rscript ./Scripts/genotype_simulation/generate_eur_id_file.R

## Download bcftools from:
##   http://www.htslib.org/download/
## and compile

## Download htslib from:
##   http://www.htslib.org/download/
## and perform the same steps as for bcftools. This contains a few general
## functions that are used in the bcftools universe.

## Note the "|" operator is the pipe operator - it feeds the output of the left
## side into the command on the right side (so we don't have to store interim
## files)

## Note the \ at the end of some lines is so the command is broken into separate
## lines to improve readability. There should be no character after the \ (even
## a space)

################################################################################
### bcftools and subsetting full dataset                                     ###
################################################################################

## These bcftools commands first subset to only the sample IDs contained within
## 1000g_eur_ids.txt, then filter to only SNPs where allele frequency > 0, then it
## omits missing genotypes (i.e. NA's), then compresses the output using the
## special (used for VCFs) bgzip compression and stores it in an output file of
## type ".vcf.gz".

./Programs/bcftools view -S ./Output/1000g_eur_ids.txt --force-samples ./Data/1000g_grch38/ALL.chr19_GRCh38.genotypes.20170504.vcf.gz |
  ./Programs/bcftools filter -i 'AF>=0.01 && EUR_AF>=0.01' |
  ./Programs/bcftools view -g  ^miss |
  ./Programs/bgzip > ./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz

## PLINK and other programs sometimes require an index, so this generates one
## based on the previously outputted VCF file.

./Programs/bcftools index ./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz

################################################################################
### Using PLINK to convert VCF to HAPGEN2 formats                            ###
################################################################################

## Download and decompress PLINK 2 from:
##   https://www.cog-genomics.org/plink/2.0/
## If you download the binary version at the top, no compilation is required.

## More info about the file formats that PLINK 2 can handle and what will be
## generated for HAPGEN2:
##   https://www.cog-genomics.org/plink/2.0/formats

## Use PLINK 2 to remove non-biallelic SNPs and then convert the VCF into .haps
## and .legend files.

## PLINK2 is used because PLINK 1.9 cannot create haps/legend files.

./Programs/plink2 \
  --vcf ./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz \
  --max-alleles 2 \
  --export hapslegend \
  --out ./Output/genotype_simulation/1000g_chr19_eur_subset_for_hapgen

## Remove unrequired sample file that PLINK makes

rm ./Output/genotype_simulation/1000g_chr19_eur_subset_for_hapgen.sample

## To make bim file for lassosum (not needed anymore)

# ./Programs/plink2 \
#   --vcf ./Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz \
#   --snps-only \
#   --max-alleles 2 \
#   --make-just-bim \
#   --out ./Output/genotype_simulation/1000g_chr19_eur_subset_bim


################################################################################
### HAPGEN2 simulations - single                                             ###
################################################################################

## Recombination rate .map file from:
##   https://github.com/naumanjaved/fingerprint_maps/blob/master/references/genetic_map_hg38.tar.gz
## original from:
##   http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/

## Download HAPGEN2 from:
##   https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html
## It contains the program already compiled so no compilation is required.

## For HAPGEN2 options (e.g. -m -l -h, etc):
##   https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html#Running_HAPGEN_top

## Run HAPGEN2 on previously created files (.map, .legend, .haps) and generate
## 1000 people (-n 1000). One SNP has to be chosen as a "disease SNP" (-dl),
## but this can be set to have a relative risk of 1 for both the reference and
## alternate allele (1 1).

./Programs/hapgen2 \
  -m ./Data/1000g_grch38_map/plink.chr19.GRCh38.map \
  -l ./Output/genotype_simulation/1000g_chr19_eur_subset_for_hapgen.legend \
  -h ./Output/genotype_simulation/1000g_chr19_eur_subset_for_hapgen.haps \
  -o ./Output/genotype_simulation/1000g_chr19_eur_simulated \
  -n 100 \
  -dl 69984 0 1 1 \
  -Ne 11418


################################################################################
### HAPGEN2 simulations - multiple sims for sample sizes                     ###
################################################################################

## Sample sizes to simulate
n_to_sim=(50 100 250 500 750 1000 1500 2000)

## Loop over each sample size (n)
for n in "${n_to_sim[@]}"
do
  echo N=$n

  ## Run simulation 100 times for each sample size
  for i in {1..100}
  do
    echo -n "  $i"
    start=`date +%s`

    ## Same HAPGEN code as previous single run
    ./Programs/hapgen2 \
      -m ./Data/1000g_grch38_map/plink.chr19.GRCh38.map \
      -l ./Output/genotype_simulation/1000g_chr19_eur_subset_for_hapgen.legend \
      -h ./Output/genotype_simulation/1000g_chr19_eur_subset_for_hapgen.haps \
      -o ./Output/genotype_simulation/rep_sims/n_${n}/1000g_sim_n${n}_${i} \
      -n ${n} \
      -dl 69984 0 1 1 \
      -Ne 11418 > /dev/null

    ## Extract the first 5 columns/fields from the gen file (contains metadata about variants, e.g.
    ## rsid, alt/ref, position)
    cut -d' ' -f1-5 Output/genotype_simulation/rep_sims/n_${n}/1000g_sim_n${n}_${i}.controls.gen | \
      gzip > Output/genotype_simulation/rep_sims/n_${n}/1000g_sim_n${n}_${i}.controls.info.gz

    ## Remove files that are not required for subseqeunt analysis
    rm Output/genotype_simulation/rep_sims/n_${n}/1000g_sim_n${n}_${i}.cases.*
    rm Output/genotype_simulation/rep_sims/n_${n}/1000g_sim_n${n}_${i}.controls.gen
    rm Output/genotype_simulation/rep_sims/n_${n}/1000g_sim_n${n}_${i}.controls.sample
    rm Output/genotype_simulation/rep_sims/n_${n}/1000g_sim_n${n}_${i}.legend

    ## Compress simulated haplotypes file
    gzip Output/genotype_simulation/rep_sims/n_${n}/1000g_sim_n${n}_${i}.controls.haps

    end=`date +%s`
    runtime=$((end-start))

    echo -e $(date +%F\ %T) '\t' ${n} '\t' ${i} '\t' ${runtime} >> Output/genotype_simulation/rep_sims/time_taken.txt
  done
  echo ""
done

