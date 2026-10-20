#!/bin/bash
## 11_hafpipe_validation — submit from the repository folder (needs only the input files):  sbatch run_11_hafpipe_validation.sh
#SBATCH --partition=cas
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --time=0-02:00:00
#SBATCH --job-name="11_hafpipe_validation"
#SBATCH --output=11_hafpipe_validation_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rene.shahmohamadloo@wsu.edu

module load r/4.5.1
cd $SLURM_SUBMIT_DIR
echo "Started on $(hostname) at $(date)"
Rscript --vanilla 11_hafpipe_validation.R && echo "R finished OK - all checks passed at $(date)" || { echo "R FAILED or a check failed at $(date) - see the report in the output folder"; exit 1; }
