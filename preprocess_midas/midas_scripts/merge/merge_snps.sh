# !/bin/bash

# Perform pooled-sample core-genome SNP calling
# Usage: bash scripts/merge/snps.sh outdir comma-separate-samples-paths optional_args

$HOME/MIDAS/scripts/merge_midas.py snps "$2" -i "$3" -t list "${@:4}"
echo "Completed MIDAS merge snps processing for $1"