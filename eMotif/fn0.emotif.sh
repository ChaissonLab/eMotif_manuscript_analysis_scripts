#!/usr/bin/env bash
#SBATCH --ntasks=16
#SBATCH --time=168:00:00
#SBATCH --mem=70000
#SBATCH --partition=qcb
#SBATCH --account=mchaisso_100
#SBATCH -N 1
#SBATCH --job-name=emotif
#SBATCH --output=slurm.%A_%a.%x.log
###SBATCH --constraint=xeon-2665,avx
###SBATCH --exclude=b10-10
###SBATCH --mail-type=ALL
###SBATCH --mail-user=tsungyul@usc.edu
###SBATCH --array=0,1

#set -eu
set -e
source ~/.bashrc
conda activate nb1
#module load gcc usc

ts=( $(ls -lh /project/mchaisso_100/cmb-17/vntr_genotyping/gtex/eqtl/ResMat/*.ResMat.pickle | awk 'BEGIN {OFS="\t"} {split($9,v1,"/"); split(v1[9],v2,"."); print substr($5,1,length($5)-1), v2[1]}'  | sort -k1,1nr | cut -f 2) )
#ts=( $(cat /project/mchaisso_100/cmb-17/vntr_genotyping/gtex/expression/alltissue.txt) )
t=${ts[$SLURM_ARRAY_TASK_ID]}
sd=/project/mchaisso_100/cmb-17/vntr_genotyping/rpgg2_k21_84k/hgsvc/gtex/emotif/log10/test


date
echo $t
$sd/emotif_single_tissue.gtex.py $t
date
