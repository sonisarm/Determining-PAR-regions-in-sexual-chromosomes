##########################################################
### Author: Sarmiento Cabello, Sonia                   ###
### Version: 1.0.                                      ###
### Objective: Plot coverage distribution males vs     ###
###           females.                                 ###
##########################################################

# Load library
library(vcfR)
library(tidyverse)

##############################
##  Section for homogametic ##
##############################
# Load data
vcf <- read.vcfR('dpf_Males_int_Z_DP.vcf.gz')
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
result_male <- matrix(nrow = nrow(dp_samples), ncol=ncol(dp_samples))
colnames(result_male) <- colnames(dp_samples)
rownames(result_male) <- rownames(dp_samples)

for(i in cov_final$Library_name){
  result_male[,i] <- dp_samples[,i]/cov_final$Coverage[cov_final$Library_name==i]
}

# Write results in a table
#write.table(result_male, 'male.txt')

# Calculate mean values per individual and per SNP, ommiting any Nas first. 
result_male <- na.omit(result_male)
mean_per_indv_males <- colMeans(result_male)
mean_per_snp_males <- rowMeans(result_male)
write.table(mean_per_indv_males, 'mean_SD_per_indv_male.txt', quote =F)
write.table(mean_per_snp_males, 'mean_SD_per_snp_male.txt', quote =F)

# Plotting SD distribution of males
# Plot distribution per individual
jpeg(paste0('Depth_distribution_corrected_by_coverage_in_males_perIndv.jpg'))
par(mfrow=c(1,1))
hist(mean_per_indv_males, main = 'Coverage distribution in sexual chr (males)', xlab='Coverage Sex CHR / Coverage Whole Genome', breaks = 10, xlim = c(0,1))
dev.off()



################################
##  Section for heterogametic ##
################################
# Load data
vcf <- read.vcfR('dpf_Females_int_Z_DP.vcf.gz')
coverage <- read.table('coverage_refpanel.txt', header=T)
coverage <- na.omit(coverage)   #Remove indv with no coverage data (N=7)
female <- read.table('List_females.txt', header=F) # List of females that we analized

#Extract coverages of our samples
cov_final <- coverage %>% filter(Library_name %in% female$V1) #211 females

# Extract read depth
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)

#Visualize the data
dp[1:50, 1:5]
#Extract samples of interest
dp_samples <- dp # only for females. 


# A peak of the read depth of some samples:
#subsample <- dp[1:50, 1:274]
#par(mar=c(12,4,4,2))
#boxplot(dp_samples, col=2:8, las=3)
#title(ylab = "Depth (DP)")


# Divide the individual SD y individual coverage
result_female <- matrix(nrow = nrow(dp_samples), ncol=ncol(dp_samples))
colnames(result_female) <- colnames(dp_samples)
rownames(result_female) <- rownames(dp_samples)

for(i in cov_final$Library_name){
  result_female[,i] <- dp_samples[,i]/cov_final$Coverage[cov_final$Library_name==i]
}

# Write results in a table
#write.table(result_female, 'male.txt')

# Calculate mean values per individual and per SNP, ommiting any Nas first. 
result_female <- na.omit(result_female)
mean_per_indv_females <- colMeans(result_female)
mean_per_snp_females <- rowMeans(result_female)
mean_per_snp_females[1:50]
write.table(mean_per_indv_females, 'mean_SD_per_indv_females.txt')
write.table(mean_per_snp_females, 'mean_SD_per_snp_females.txt')

# Plotting SD distribution of males
# Plot distribution per individual
jpeg(paste0('Depth_distribution_corrected_by_coverage_in_females_perIndv.jpg'))
par(mfrow=c(1,1))
hist(mean_per_indv_females, main = 'Coverage distribution in sexual chr (females)', xlab='Coverage Sex CHR / Coverage Whole Genome', breaks = 5, xlim = c(0,1))
dev.off()
