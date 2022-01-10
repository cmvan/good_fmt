#!/bin/bash

# Prefetches and downloads fastq files from SRA and gzips fastq files
# Usage: bash download_exchange.sh data_path person sample_acc library_layout

# ONLY WORKS FOR SINGLE READS

export data="$1"
export person="$2"
export sample="$3"
export type="$4"

for run_acc in $(python bash_download_exchange.py "$data" print_run_acc "$person" "$sample" "$type"); do
    out="outdir/${type}/${person}/${sample}"
    path="prefetch_loads/${type}/${run_acc}"
    run_file="${out}/${run_acc}.fastq"
    comp_file="${out}/${sample}_compiled.fastq"

    prefetch "$run_acc" -o "$path"
    fasterq-dump "$run_acc" --outdir "$out"
    cat "$run_file" >> "$comp_file"
done
gzip -r "$out"
