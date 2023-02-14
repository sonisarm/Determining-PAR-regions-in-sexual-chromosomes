# Identifying-PAR-regions-in-sexual-chromosomes

Author: Sonia Sarmiento

### Description: This repository aims at defining PAR and non-PAR regions of sexual chromosomes for diploid species. 

#### Coding: I used mainly a HPC cluster with bash scripts and R code at the end for analysis and plotting.

## Introduction
Welcome to my repository for identifying PAR (Pseudoautosomal Regions) and non-PAR regions of sexual chromosomes!
In diploid species, sexual chromosomes are typically different in size and gene content. However, there are regions of homology between the heterozygote sexual chromosomes known as PAR, in which, accordingly, genes are present on both X and Y chromosomes. This repository aims to provide tools and algorithms for identifying PAR regions in the sexual chromosomes. I developed this code in order to phase and obtain an accurate reference panel that I used for imputation. However, identifying PAR regions can also be used for instance to help understand sexual chromosome divergence and to distinguish between X-linked and Y-linked genes in medial genetics. 

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

### Step 2: Combining GVCFs and joint genotyping.
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
  
# Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
gatk --java-options "-Djava.io.tmpdir=${myDIR}/tmp/" GenotypeGVCFs -R ${ref} -V ${intermediate_file} -O ${output} --tmp-dir ${myDIR}/tmp/

````
* Output: A final VCF in which all samples have been jointly genotyped


### Step 3: Filtering VCFs
Filter VCF according to GATK best practices.
* Input:
1) ```${ref}```  --> Reference genome (e.g. from NCBI)
2) ```${input_VCF}```  --> VCF (output of previous step)
3) ```${maked_regions}``` --> List of maked regions per chromosome
The mask file should contain the start and end coordinates of each masked region. 

4) ```${high_cutoff}``` and ```${low_cutoff}``` --> Cutt off for sequencing depth filtering. The high cutoff depends on the average sequencing depth (e.g. high cutoff should not be double than the average depth). The lower cutoff depends on the ploidy of the chromosome. 

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
bcftools view -O z -o ${intermediate_VCF3} -R ${maked_regions} ${intermediate_VCF2}
	
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
**IMPORTANT**: when filtering for sequencing depth, the haploid individual's lower cutoff (```${low_cutoff}```) has to be at least halved in the sexual chromosome only (e.g. males XY in human, females ZW in birds). For instance, 2.0 for female birds (ZW) and 5.0 for male birds (ZZ). 

* Output: Filtered VCF (GATK and read-depth) with masked regions removed. As the lower cutoff is different for females/males due to heterozygosity in sexual chromosomes, the last step is done separately for each sex and we obtain one VCF per sex. 


### Step 4: Extracting read depth
The PAR and non-PAR regions can be detected by differences in sequencing depth between regions in the sexual chromosome. To agilize further analysis, we extract only this information from the VCF. 
* Input:
1) ```${input_VCF}```  --> Filtered VCF (output of previous step)

* Script:
```bash
bcftools annotate -i 'TYPE="snp" && FORMAT/DP!="."' -x^FORMAT/DP -O z -o ${output} ${input_VCF} 
````
* Output: Zipped VCF with read-depth per individual and SNP.

### Step 5: Analyse read-depth information
After obtaining a VCF with read-depth information, we analyse this data in R. The R code provided was written to analyse two VCF (one per sex) separately but can also be modified to get information from one VCF with both female+male samples. 
* Input:
1) VCF (one per sex) with read-depth information per individual (column) and SNP (rows).
2) Coverage of BAM file per individual (TXT file - col1=ind; col2=coverage, col3;sex (or one file per sex))
3) List of individuals to analyse

* Script: ```1_Coverage_Analysis.R```
* 
This script divides each SNP position for each individual by the specific-coverage of the individual (for the whole genome). 

* Output: Distribution plots of sequencing depth (sexual chr) / coverage (whole genome) per SNP and individual.

On the one hand, I expect that the DIPLOID homogametic individual for the sexual chromosome present a distribution of around 1 (e.g.```SD_distribution_homozygous_sex.jpg``` shows the distribution of 280 individuals) . On the other hand, I expect that the distribution for the heterogametic individual follows a bimodal distribution between 0.5 and 1, corresponding to the non-PAR and PAR regions, accordingly. The non-PAR (haploid) region will be 0.5 when divided by the coverage of the whole genome because its haploid. 


