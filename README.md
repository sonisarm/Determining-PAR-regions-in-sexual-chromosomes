# Determining-PAR-regions-in-sexual-chromosomes

Author: Sonia Sarmiento

### This repository aims at defining PAR and non-PAR regions of sexual chromosomes. I used mainly a HPC cluster with bash scripts and R code at the end for plotting.

## Introduction
Welcome to my repository for identifying PAR (Pseudoautosomal Regions) and non-PAR regions of sexual chromosomes!
Sexual chromosomes determine biological sex of an organisms (e.g. in mammals X and Y chromosomes). In many species, sexual chromosomes are differente in size and gene content. However, there are regions of homology between the X and Y chromosomes known as PAR, where, accordingly, genes are present on both X and Y chromosomes. This repository aims to provide tools and algorithms for identifying PAR regions in the sexual chromosomes. I developed this code in order to phase and obtain an accurate reference panel that I used for imputation. However, identifying PAR regions can also be used for instance to help understand sexual chromosome divergence and to distinguish between X-linked and Y-linked genes in medial genetics. 

The workflow includes SNP calling (both with and without ploidy), filtering steps and finally plotting sequencing depth distributions. 

## Workflow

### Step 1: SNP Calling
* Input: 
1) ```${ref}```  --> Reference genome (e.g. from NCBI)
2) ```${BAMFile}```  --> BAM file from individual
3) ```Z_ss.list``` --> list of sexual chromosomes/scaffolds (one per line) (e.g. ChrX for humans)

* Script: 
```bash
gatk HaplotypeCaller \
   -R {ref} \
   -I ${BAMFile} \
   -L Z_ss.list \
   --pcr-indel-model NONE \
   -ERC GVCF \
   --sample-ploidy 1 \ #INCLUDE THIS LINE CODE TO PRECISE THE CHROMOSOME IS HAPLOID. 
   -O ${output}.g.vcf.gz 
```
Run this script for the sexual chromosome one time specifying for haploidy vs without this part of the code. 

* Output: GVCF containing called SNPs.

### Step 2: Combining GVCFs and genotype
* Input: 
1) ```${ref}```  --> Reference genome (e.g. from NCBI)
2) ```GVCF_files.list```  --> List of GVCFs to combine (one per sample).
3) ```Z_ss.list``` --> list of sexual chromosomes/scaffolds (one per line) (e.g. ChrX for humans)

* Script:
```bash
# Combine GVCFs
gatk CombineGVCFs \
         -R ${ref} \
         -V GVCF_files.list \
         -L Z_ss.list \
         -O ${intermediate_file}
         
gatk --java-options "-Djava.io.tmpdir=${myDIR}/tmp/" GenotypeGVCFs -R ${ref} -V ${intermediate_file} -O ${output} --tmp-dir ${myDIR}/tmp/

````
* Output: One VCF with called SNPs for all samples


### Step 3: Filtering VCFs

* Input: VCF file ($vcf) and reference genome ($ref)
* Script: ```1_FilteringVCFs.sh```
* Output: filtered VCF

**NOTE**: when filtering for sequencing depth, the haploid individual's lower cutoff has to be at least halved in the sexual chromosome only (e.g. males XY in human, females ZW in birds) 

