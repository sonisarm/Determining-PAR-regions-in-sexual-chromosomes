##########################################################
### Author: Sarmiento Cabello, Sonia                   ###
### Version: 1.0.                                      ###
### Objective: Plot coverage distribution males vs     ###
###           females.                                 ###
##########################################################

rm(list=ls())
setwd('/Users/soniasarmiento/Downloads/BarnOwl')

# Load library
library(vcfR)
library(tidyverse)

########################
##  Section for males ##
########################
# Load data
vcf <- read.vcfR('MalesOnly_DP.vcf.gz')
coverage <- read.table('coverage_refpanel.txt', header=T)
coverage <- na.omit(coverage)   #Remove indv with no coverage data (N=7)
#female <- read.table('RefPanel_femaleonly.txt', header=T)
male <- read.table('RefPanel_male.only.txt', header=T)  #get male indvs in refpanel
#Extract coverages of our samples
cov_final <- coverage %>% filter(Library_name %in% male$Library_name) #274 males

# Extract read depth
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)

#Visualize the data
dp[1:50, 1:5]
#Extract samples of interest
dp_samples <- dp[,cov_final$Library_name]

# A peak of the read depth of some samples:
#subsample <- dp[1:50, 1:274]
#par(mar=c(12,4,4,2))
#boxplot(dp_samples, col=2:8, las=3)
#title(ylab = "Depth (DP)")


# Divide the individual SD y individual coverage
result <- matrix(nrow = nrow(dp_samples), ncol=ncol(dp_samples))
colnames(result) <- colnames(dp_samples)
rownames(result) <- rownames(dp_samples)

for(i in cov_final$Library_name){
  result[,i] <- dp_samples[,i]/cov_final$Coverage[cov_final$Library_name==i]
}

# Write results in a table
#write.table(result, 'male.txt')

# Calculate mean values per individual and per SNP, ommiting any Nas first. 
result <- na.omit(result)
mean_per_indv <- colMeans(result)
mean_per_snp <- rowMeans(result)
write.table(mean_per_indv, 'mean_SD_per_indv.txt')
write.table(mean_per_snp, 'mean_SD_per_snp.txt')

# calculate min and max value
min_ind <- min(mean_per_indv)
max_ind <- max(mean_per_indv)
min_snp <- min(mean_per_snp)
max_snp <- max(mean_per_snp)

# Plotting SD distribution of males
# Plot distribution per individual
jpeg(paste0('Depth_distribution_corrected_by_coverage_in_males_perIndv.jpg'))
par(mfrow=c(1,1))
hist(mean_per_indv, main = 'SD/coverage distribution in males sexual chr', xlab='SD/coverage')
dev.off()
