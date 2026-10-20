#!/bin/bash
## Watterson's θ inputs for 12_diversity_theta.R: grenedalf diversity on the HAF-pipe allele-frequency
## table, per chromosome and genome-wide, at the lower (500x) and upper (3,500x) effective coverage.
## Writes diversity_{chrom,genome}_eff{500,3500}.csv; copy them to the input folder (RF_RAW_INPUT).
## Usage:  bash grenedalf_theta.sh <HAF-pipe frequency table, e.g. all.csv.gz> <output folder>
set -e
FREQ=$1; OUT=${2:-.}
mkdir -p "$OUT"
for EFF in 500 3500; do
	for WIN in chromosomes genome; do
		D=$OUT/tmp_${WIN}_${EFF}; mkdir -p "$D"
		grenedalf diversity --frequency-table-path "$FREQ" --window-type $WIN --pool-sizes 240 \
			--filter-sample-min-count 2 --frequency-table-int-factor $EFF --measure theta-watterson \
			--allow-file-overwriting --out-dir "$D"
		NAME=$([ $WIN = chromosomes ] && echo chrom || echo genome)
		mv "$D"/diversity.csv "$OUT/diversity_${NAME}_eff${EFF}.csv" && rmdir "$D"
	done
done
ls -l "$OUT"/diversity_*_eff*.csv
