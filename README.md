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
         -O ${intermediate_GVCF}
         
gatk --java-options "-Djava.io.tmpdir=${myDIR}/tmp/" GenotypeGVCFs -R ${ref} -V ${intermediate_file} -O ${output} --tmp-dir ${myDIR}/tmp/

````
* Output: One VCF with called SNPs for all samples


### Step 3: Filtering VCFs
Filter VCF according to GATK best practices.
* Input:
1) ```${ref}```  --> Reference genome (e.g. from NCBI)
2) ```${input_VCF}```  --> VCF (output of previous step)
3) ```${maked_regions}``` --> List of maked regions per chromosome 
4) ```${high_cutoff}``` and ```${low_cutoff}``` --> 

* Script:
```bash
#Techfilters SNPs selection
gatk VariantFiltration \
     	-R ${ref} \
     	-V ${input_VCF} \
     	-O ${intermediate_VCF} \
     	--filter-expression "QD < 2.0" --filter-name "QD2" \
     	--filter-expression "FS > 60.0" --filter-name "FS60" \
     	--filter-expression "SOR > 3.0" --filter-name "SOR3" \
     	--filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS-8" \
     	--filter-expression "MQRankSum < -12.5" --filter-name "MQRS-12.5" \
     	--filter-expression "MQ < 40.0" --filter-name "MQ40"

#Subsetting
gatk SelectVariants \
     	-R ${ref} \
     	-V ${intermediate_VCF} \
     	-O ${intermediate_VCF2} \
     	--exclude-filtered true \
     	--exclude-non-variants

# Remove masked regions
bcftools view -O z -o ${intermediate_VCF3} -R ${maked_regions} \
	${intermediate_VCF2}
	
# Filter for read depth
gatk VariantFiltration \
        -V ${intermediate_VCF3} \
        -G-filter "DP > ${high_cutoff}" \
        -G-filter-name 'GDPhigh' \
        -G-filter "DP < ${low_cutoff}" \
        -G-filter-name 'GDPlow' \
        -O  ${output} \
        --set-filtered-genotype-to-no-call true
   
`````
* Output: Filtered VCF (GATK and read-depth) with masked regions removed.

**IMPORTANT**: when filtering for sequencing depth, the haploid individual's lower cutoff (```${low_cutoff}```) has to be at least halved in the sexual chromosome only (e.g. males XY in human, females ZW in birds). For instance, 2.0 for female birds (ZW) and 5.0 for male birds (ZZ). 

