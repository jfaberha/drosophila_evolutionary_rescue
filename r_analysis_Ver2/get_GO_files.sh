#!/bin/bash
## Download the GO ontology and FlyBase GO annotation used by 09_GO_enrichment.R.
## Run on the Kamiak LOGIN node (compute nodes may not reach the internet):  bash get_GO_files.sh
## Table S8 was made with files downloaded on 12 Oct 2026; newer files can change the table slightly.
DIR=${RF_GO_FILES:-data/GO}
mkdir -p "$DIR" && cd "$DIR"
wget -O go-basic.obo https://purl.obolibrary.org/obo/go/go-basic.obo
wget -O fb.gaf.gz https://current.geneontology.org/annotations/fb.gaf.gz
ls -lh go-basic.obo fb.gaf.gz
echo "Files in $DIR - set RF_GO_FILES=$DIR before running 09 if this is not data/GO"
