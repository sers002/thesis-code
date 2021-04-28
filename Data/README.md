# Data

This directory contains the *unmodified* input datasets that are used through the framework. The files in these directories are never modified to maintain reproducibility. 

### `1000g_grch38`

`1000g_grch38` contains the genotype files from the 1000 Genomes project with GRCh38 positions. These files are in compressed VCF format (.vcf.gz). Other chromosomes can be used instead of chromosome 19, but script files in the `Scripts/genotype_simulation` directory would need to be updated.

Source: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/

#### Files: 

* `ALL.chr19_GRCh38.genotypes.20170504.vcf.gz`

### `1000g_grch38_map`

`1000g_grch38_map` contains the recombination maps in PLINK format for GRCh38. These files are sourced from the Beagle program resources.

Source: http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/

#### Files:

* `plink.chr19.GRCh38.map`

### `1000g_metadata`

`1000g_metadata` contains the metadata (population, superpopulation, etc.) for the individuals sampled in the 1000 Genomes project. 

Source: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/

#### Files:

* `20130606_g1k.ped`

### `lung_eqtl_data`

`lung_eqtl_data` contains the eQTL results from GTEx after processing by the CoDeS3D pipeline. These files are output by the CoDeS3D pipeline using the GTEx GRCh38 dataset and Schmitt lung Hi-C data.

Source: O'Sullivan Lab (Liggins Institute)

#### Files:

* `eqtls.txt`
* `genes.txt`
* `snps.txt`
