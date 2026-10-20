#!/bin/bash
## 05_convergent_candidates — submit from the repository folder:  sbatch run_05_convergent_candidates.sh
#SBATCH --partition=cas
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH --time=0-12:00:00
#SBATCH --job-name="05_convergent_candidates"
#SBATCH --output=05_convergent_candidates_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rene.shahmohamadloo@wsu.edu

module load r/4.5.1
cd $SLURM_SUBMIT_DIR
echo "Started on $(hostname) at $(date)"
Rscript --vanilla 05_convergent_candidates.R && echo "R finished OK - all checks passed at $(date)" || { echo "R FAILED or a check failed at $(date) - see the report in the output folder"; exit 1; }
