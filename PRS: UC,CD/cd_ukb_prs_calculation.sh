#!/usr/bin/env bash
#SBATCH --account=zivlab
#SBATCH --partition=Ziv,common
#SBATCH --job-name=cdukb_prs
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pooja.middha@ucsf.edu
#SBATCH --ntasks=4
#SBATCH --mem=8G
#SBATCH --time=23:00:00
#SBATCH --output=/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Output/cd_ukb_prs.log

pwd; hostname; date

module load CBI plink2
cd /zivlab/data3/pmiddha

prspath="/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB"

plink2 \
  --bfile $prspath/Data/genetic/CD/cd_ukb_allchr \
  --score $prspath/Results/CD/cd_ldpred2_auto_model.txt 3 7 8 header list-variants ignore-dup-ids cols=+scoresums \
  --memory 16000 \
  --threads 3 \
  --out $prspath/Results/CD/prs_score/cd_ukb_prs
