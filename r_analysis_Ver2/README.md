# rudflies_2023 — analysis code

Code for "Contingency and determinism in rapid adaptation shape evolutionary rescue and
extinction in *Drosophila*" (Shahmohamadloo, Faber-Hammond, Gabidulin, Clare & Rudman).
All analyses are in R. Each module reads `00_config.R` (paths, excluded sample, cage groups,
shared helpers), writes to its own folder under the output folder, and ends by checking its
key results against `expected_values.csv` (see **Checks**).

Sample 90 (unexposed cage 29, 30 September 2023; failed, low-diversity library) is excluded
from every analysis. Phenotypic and census analyses are archived separately on Dryad.

## Modules

| Module | What it does | Paper items |
|---|---|---|
| `00_config.R` | Paths, sample-90 exclusion, cage groups, helpers | — |
| `01_glm/` | All per-SNP GLMs (28 models incl. leave-one-out; quasi-Poisson, founders quasi-binomial) | inputs to everything genomic |
| `02_build_glm_master.R` | Merges GLM tables; FDR; −log10(P); 21-SNP windows; shared objects | SI 2.1.1–2.1.2 counts; inputs to 03–10 |
| `03_genome_scans.R` | Manhattan plots; \|Δp\| at outliers and candidate genes vs matched SNPs; outlier overlap test | Fig. 2D–G, 3D–I, S3, S4, S8; SI 2.1.2–2.1.4 |
| `04_rrho_spearman.R` | Spearman correlations and RRHO of August contrasts | main-text ρ; Fig. S11B–D |
| `05_convergent_candidates.R` | Convergent candidate SNP sets and gene lists; residual Manhattans; trajectories | Fig. 4D, S12; Tables S7, S10, S11; SI 2.3 |
| `06_cage_scores_founder_axis.R` | Population shifts toward the pre-adapted state | Fig. 3J, 4C; SI 2.2.3 |
| `07_concordance_LOO.R` | SNP clusters; leave-one-out and cross-cluster concordance t-tests; Euclidean distances | Tables S3, S4, S6; SI 2.4; Fig. S9 |
| `08_candidate_genes.R` | Spinosad candidate genes (Table S2) vs outliers (hypergeometric, Wilcoxon) | SI 1.8.6; SI 2.1.3–2.1.4 |
| `09_GO_enrichment.R` (+ `get_GO_files.sh`) | GO enrichment of the convergent candidate gene lists | Table S8 |
| `10_cost_temporal.R` | Cost of resistance: temporal change vs pre-adapted outliers; Spearman; top-10K vs matched SNPs | Tables S1, S9, S12, S13; Fig. S13 |
| `11_hafpipe_validation.R` | HAF-pipe and BCFtools frequencies vs DGRP2 founder frequencies | Fig. S14; SI 1.7 |
| `12_diversity_theta.R` (+ `grenedalf_theta.sh`) | Watterson's θ per chromosome and genome-wide at 500x and 3,500x effective coverage (all sampling dates); GLMs among groups on the August samples | Fig. S10; Table S5; SI 1.8.5, 2.2.1 |

Not included yet: the read-processing and HAF-pipe commands (grenepipe, HAF-pipe), the
sequencing statistics (SI 2.1.1) and the inversion analysis (SI 2.2.2).

## Inputs
By default every path is a folder inside this repository; set the `RF_*` environment
variables (or edit the PATHS block of `00_config.R`) to use other locations.

| Folder (default) | Variable | Contents |
|---|---|---|
| `data/input/` | `RF_RAW_INPUT` | HAF-pipe allele-frequency table, sample metadata, VEP/FlyCADD annotation, spinosad candidate gene list, genic-feature table, DGRP2 and BCFtools `.frq` files, grenedalf θ tables (`bash grenedalf_theta.sh`) |
| `data/glm_tables/` | `RF_GLM_TABLES` | per-SNP GLM result tables written by `01_glm/glm_models.R` |
| `data/GO/` | `RF_GO_FILES` | `go-basic.obo` and `fb.gaf.gz` (`bash get_GO_files.sh`) |
| `output/` | `RF_OUT_ROOT` | everything the modules write (one sub-folder per module) |
| `data/published_tables/` | — | published Tables S1, S3, S4, S6–S13, used for the cell-by-cell checks (S7, S10 and S11 also supply the feature annotation text, FlyBase r6.69) |

## Running (Slurm; tested on WSU Kamiak with R 4.5.1)
Submit from this folder (01 from `01_glm/`). Module 01 writes the GLM tables to
`output/01_glm/`; copy them to `data/glm_tables/` (or set `RF_GLM_TABLES`) before running 02.
`sbatch --export=ALL,MODE=validate 01_glm/run_01_glm.sh` instead checks 3,000 SNPs per model
against existing tables in a few minutes. The run scripts email the corresponding author;
change `--mail-user` before running. `bash submit_all.sh` submits 02–12 with the
dependencies below; or submit them one by one:
```
sbatch 01_glm/run_01_glm.sh                # refits all 28 GLMs (array job; hours to days per model)
sbatch run_02_build_glm_master.sh          # first                                   (~9 min)
sbatch run_03_genome_scans.sh              # after 02                                (~10 min)
sbatch run_04_rrho_spearman.sh             # after 02                                (~4 min)
sbatch run_05_convergent_candidates.sh     # after 02                                (~13 min)
sbatch run_06_cage_scores_founder_axis.sh  # after 02                                (~2 min)
sbatch run_08_candidate_genes.sh           # after 02                                (~15 min)
sbatch run_10_cost_temporal.sh             # after 02                                (~43 min)
sbatch run_07_concordance_LOO.sh           # after 05 (Fig. S9 SNP list)             (~28 min)
sbatch run_09_GO_enrichment.sh             # after 05; GO files from get_GO_files.sh (~1 min)
sbatch run_11_hafpipe_validation.sh        # any time (input files only)             (~2 min)
sbatch run_12_diversity_theta.sh           # any time (grenedalf tables)             (~1 min)
```

## Checks
Each module prints its checks and writes `<module>_checks.csv`:
- `PASS` / `FAIL`: a result compared with the value reported in the paper (or every cell of
  a published table in `data/published_tables/`). A FAIL marks the job as failed; all
  outputs are still written.
- `TEXT` / `TEXT?`: a result compared with a reported value that has not been verified;
  `TEXT?` means they differ. This is not a failure.

## Random draws
Matched-SNP and matched-gene draws (modules 03, 07, 08, 10) use fixed seeds (`SEED` in
`00_config.R`, and the seeds set in each module), so reruns give identical results.
