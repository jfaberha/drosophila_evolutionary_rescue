#!/bin/bash
## 11_hafpipe_validation — submit from the repository folder:  sbatch run_11_hafpipe_validation.sh
## Needs the dm3-to-dm6 chain in <output folder>/liftover first (on the login node: bash get_liftover_chain.sh "$RF_OUT_ROOT").
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
Rscript --vanilla liftover_dgrp2_frq.R && Rscript --vanilla 11_hafpipe_validation.R && echo "R finished OK - all checks passed at $(date)" || { echo "R FAILED or a check failed at $(date) - see the report in the output folder"; exit 1; }
