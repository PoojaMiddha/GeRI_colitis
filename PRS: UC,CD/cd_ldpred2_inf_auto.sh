#!/usr/bin/env bash
#SBATCH --account=zivlab
#SBATCH --partition=Ziv,common
#SBATCH --job-name=cd_inf_auto
#SBATCH --output=/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Output/cd_ldpred2_inf_auto.log
#SBATCH --time=1-23:55:00
#SBATCH --mem=100G
#SBATCH --ntasks=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pooja.middha@ucsf.edu

pwd;hostname;date

module load CBI r/4.2.1

  Rscript /zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Code/genetic/CD/cd_ldpred2_inf_auto.R
