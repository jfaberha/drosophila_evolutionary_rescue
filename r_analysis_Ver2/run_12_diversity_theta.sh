#!/bin/bash
## 12_diversity_theta — submit from the repository folder (needs only the input files, incl. the grenedalf tables):  sbatch run_12_diversity_theta.sh
#SBATCH --partition=cas
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=0-01:00:00
#SBATCH --job-name="12_diversity_theta"
#SBATCH --output=12_diversity_theta_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rene.shahmohamadloo@wsu.edu

module load r/4.5.1
cd $SLURM_SUBMIT_DIR
echo "Started on $(hostname) at $(date)"
Rscript --vanilla 12_diversity_theta.R && echo "R finished OK - all checks passed at $(date)" || { echo "R FAILED or a check failed at $(date) - see the report in the output folder"; exit 1; }
