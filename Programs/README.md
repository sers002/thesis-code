# Programs Required

This directory contains the programs required to manipulate and simulate the genomic datasets. Each heading in the following list is what the program should be named (which are the default after compilation).

## `bcftools`

`bcftools` is used to manipulate genomic files in the VCF or BCF formats. 

Source: http://www.htslib.org/download/

## `bgzip`

`bgzip` is use to compress the genomic file that results from `bcftools`. It is contained in the `htslib` library (which contains many other utilities that are not used in this repository).

Source: http://www.htslib.org/download/

## `hapgen2`

`hapgen2` is the simulator used to generate genomic datasets based on a reference population when using the empirical method. 

Source: https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html

## `plink`

`plink` is used for manipulation and analysis of genomic datasets. 

Source: https://www.cog-genomics.org/plink/

## `plink2`

`plink2` is an updated version of `plink`. 

Source: https://www.cog-genomics.org/plink/2.0/

