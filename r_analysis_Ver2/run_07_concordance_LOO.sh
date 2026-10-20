#!/bin/bash
## 07_concordance_LOO — submit from the repository folder after 02 and 05:  sbatch run_07_concordance_LOO.sh
#SBATCH --partition=cas
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH --time=0-12:00:00
#SBATCH --job-name="07_concordance_LOO"
#SBATCH --output=07_concordance_LOO_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rene.shahmohamadloo@wsu.edu

module load r/4.5.1
cd $SLURM_SUBMIT_DIR
echo "Started on $(hostname) at $(date)"
Rscript --vanilla 07_concordance_LOO.R && echo "R finished OK - all checks passed at $(date)" || { echo "R FAILED or a check failed at $(date) - see the report in the output folder"; exit 1; }
