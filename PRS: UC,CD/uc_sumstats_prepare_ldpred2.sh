#!/usr/bin/env bash
#SBATCH --account=zivlab
#SBATCH --partition=Ziv,common
#SBATCH --job-name=ucsumstats_prep
#SBATCH --output=/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Output/uc_sumstats_hapmapmatch_prepare.log
#SBATCH --time=23:30:00
#SBATCH --mem=100G
#SBATCH --ntasks=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pooja.middha@ucsf.edu

pwd;hostname;date

module load CBI r/3.6.3

Rscript /zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Code/genetic/UC/uc_sumstats_prepare_ldpred2.R
