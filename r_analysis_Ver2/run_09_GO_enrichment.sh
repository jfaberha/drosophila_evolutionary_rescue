#!/bin/bash
## 09_GO_enrichment — submit after 05 (and get_GO_files.sh if the GO files are missing):  sbatch run_09_GO_enrichment.sh
#SBATCH --partition=cas
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=0-02:00:00
#SBATCH --job-name="09_GO_enrichment"
#SBATCH --output=09_GO_enrichment_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rene.shahmohamadloo@wsu.edu

module load r/4.5.1
cd $SLURM_SUBMIT_DIR
echo "Started on $(hostname) at $(date)"
Rscript --vanilla 09_GO_enrichment.R && echo "R finished OK - all checks passed at $(date)" || { echo "R FAILED or a check failed at $(date) - see the report in the output folder"; exit 1; }
