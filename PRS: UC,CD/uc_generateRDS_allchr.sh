#!/usr/bin/env bash
#SBATCH --account=zivlab
#SBATCH --partition=Ziv,common
#SBATCH --job-name=rds
#SBATCH --output=/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Output/uc_rds_allchr.log
#SBATCH --time=1-20:30:00
#SBATCH --mem=100G
#SBATCH --ntasks=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pooja.middha@ucsf.edu

pwd;hostname;date

module load CBI r/4.2.1

Rscript /zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Code/genetic/UC/uc_generateRDS_allchr.R
