#!/bin/bash
## Submit modules 02-12 in dependency order. Run from this folder:  bash submit_all.sh
## 'afterany' is used so that a module whose checks FAIL (exit 1) still lets later modules run.
set -e
j02=$(sbatch --parsable run_02_build_glm_master.sh)
for m in 03_genome_scans 04_rrho_spearman 06_cage_scores_founder_axis 08_candidate_genes 10_cost_temporal; do
	sbatch --parsable --dependency=afterany:$j02 run_$m.sh > /dev/null
done
j05=$(sbatch --parsable --dependency=afterany:$j02 run_05_convergent_candidates.sh)
sbatch --parsable --dependency=afterany:$j05 run_07_concordance_LOO.sh > /dev/null
sbatch --parsable --dependency=afterany:$j05 run_09_GO_enrichment.sh > /dev/null
sbatch --parsable run_11_hafpipe_validation.sh > /dev/null
sbatch --parsable run_12_diversity_theta.sh > /dev/null
echo "Submitted: 02 = job $j02; 03, 04, 06, 08, 10 wait for 02; 05 = job $j05; 07, 09 wait for 05; 11 and 12 run now."
squeue -u $USER
