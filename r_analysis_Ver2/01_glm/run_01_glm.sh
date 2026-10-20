#!/bin/bash
## 01_glm: one array task per GLM model (28 models; see `Rscript glm_models.R list`).
## Default MODE=full: refits every SNP and writes the GLM tables read by 02_build_glm_master.R
## (hours to days per model; progress is saved in chunks, so a resubmitted task resumes).
## MODE=validate: fits 3,000 SNPs spread along the genome per model and compares them with the
## tables in RF_GLM_TABLES (a few minutes per model).
## Submit from this folder:
##   sbatch run_01_glm.sh                                     # full refit of all 28 models
##   sbatch --export=ALL,MODE=validate run_01_glm.sh          # validate all 28 models
##   sbatch --array=6 run_01_glm.sh                           # one model only (index from the list)
#SBATCH --partition=cas
#SBATCH --array=0-27
#SBATCH --cpus-per-task=8
#SBATCH --mem=48GB
#SBATCH --time=7-00:00:00
#SBATCH --job-name="01_glm"
#SBATCH --output=01_glm_%a_%A.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rene.shahmohamadloo@wsu.edu

MODE=${MODE:-full}
module load r/4.5.1
cd $SLURM_SUBMIT_DIR
MODELS=($(Rscript --vanilla glm_models.R list))
M=${MODELS[$SLURM_ARRAY_TASK_ID]}
echo "Model $M ($MODE) started on $(hostname) at $(date)"
N=$([ "$MODE" = "validate" ] && echo 3000 || echo 0)
Rscript --vanilla glm_models.R $M $MODE $N $SLURM_CPUS_PER_TASK && echo "R finished OK at $(date)" || { echo "R FAILED at $(date)"; exit 1; }
