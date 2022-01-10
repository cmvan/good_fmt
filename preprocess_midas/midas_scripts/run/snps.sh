# !/bin/bash

# Profile strain-level nucleotide variants of abundant species
# Usage: bash scripts/single/snps.sh outdir sample.gz optional_args

$HOME/MIDAS/scripts/run_midas.py snps --remove_temp "$1" -1 "$2" "${@:3}"
echo "Completed MIDAS snps processing for $2"