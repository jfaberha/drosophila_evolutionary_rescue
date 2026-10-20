#!/bin/bash
## 02: build the GLM master table once and save it for later modules (~1 h).
## Submit from the repository folder:  sbatch run_02_build_glm_master.sh
#SBATCH --partition=cas
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH --time=0-08:00:00
#SBATCH --job-name="02_master"
#SBATCH --output=02_build_glm_master_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rene.shahmohamadloo@wsu.edu

module load r/4.5.1
cd $SLURM_SUBMIT_DIR
echo "Started on $(hostname) at $(date)"
Rscript --vanilla 02_build_glm_master.R && echo "R finished OK - all checks passed at $(date)" || { echo "R FAILED or a check failed at $(date) - see 02_build_glm_master_report.txt"; exit 1; }
