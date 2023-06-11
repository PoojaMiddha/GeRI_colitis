#!/usr/bin/env bash
#SBATCH --account=zivlab
#SBATCH --partition=Ziv,common
#SBATCH --job-name=ucgeri_prs_a1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pooja.middha@ucsf.edu
#SBATCH --ntasks=4
#SBATCH --mem=8G
#SBATCH --array=1-22
#SBATCH --time=23:00:00
#SBATCH --output=/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_GeRI/Output/uc_geri_prs_a1_%A_%a.log

pwd;hostname;date

cd /zivlab/data3/pmiddha
echo "Working on chromosome ${SLURM_ARRAY_TASK_ID}"

module load CBI plink2/2.00a3LM
projdir="/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs"
indat="$projdir/PRS_GeRI/Data/genetic/geri_vcf2plink"
outpath="$projdir/PRS_GeRI/Results/UC"
wtpath="$projdir/pre_2023/PRS_UKBB/Results/UC"

plink2 \
  --bfile $indat/chr${SLURM_ARRAY_TASK_ID} \
  --score $wtpath/uc_ldpred2_auto_model.txt 5 7 8 header list-variants ignore-dup-ids cols=+scoresums \
  --memory 16000 \
  --threads 3 \
  --out $outpath/uc_geri_prs_a1_chr${SLURM_ARRAY_TASK_ID}
exit
