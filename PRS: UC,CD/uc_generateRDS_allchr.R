# Script for geenrating RDS files for plink bed/bim/fam files for LDPred2
## MUST BE RUN ON CLUSTER ONLY

# FOLLOWING https://privefl.github.io/bigsnpr/articles/LDpred2.html
#Requires R 3.6.3

##install.packages(remotes)
library(remotes)
##remotes::install_github("privefl/bigsnpr")
library(bigsnpr)
library(dplyr)


# RDS file
snp_readBed("/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Data/genetic/UC/uc_ukb_allchr.bed")

# RDS file for each chromosome should be in the indir location
q("no")

