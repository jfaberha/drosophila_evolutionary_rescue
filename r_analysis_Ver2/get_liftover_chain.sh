#!/bin/bash
## Downloads the UCSC dm3 (FlyBase Release 5) -> dm6 (Release 6) liftOver chain used by liftover_dgrp2_frq.R.
## Run on the login node:  bash get_liftover_chain.sh <output folder (RF_OUT_ROOT)>   (the chain goes to <output folder>/liftover)
set -e
OUT=${1:-${RF_OUT_ROOT:-output}}/liftover
mkdir -p "$OUT"
wget -q -O "$OUT/dm3ToDm6.over.chain.gz" https://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz
ls -l "$OUT/dm3ToDm6.over.chain.gz"
