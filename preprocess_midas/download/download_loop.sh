#!/bin/bash

# Prefetches and downloads fastq files from SRA and gzips fastq files
# Usage- For downloading files: bash download_loop.sh library_type person sample_acc run_acc out_path
#        For gzipping files   : bash download_loop.sh gz out_path

if [ "$1" = 'gz' ]; then
    echo "Zipping everything ..."
    gzip -r "$2"
    echo "All done!"
else
    export path="prefetch_loads/${1}/${4}"
    export run_file="${5}/${4}.fastq"
    export comp_file="${5}/${3}_compiled.fastq"
    prefetch "$4" -o "$path"
    fasterq-dump "$path" --outdir "$5"
    echo "Concatenating file ..."

    if [ "$1" = 'paired' ]; then
        export run_file="${5}/${4}_1.fastq"
        export run_file2="${5}/${4}_2.fastq"
        export comp_file="${5}/${3}_1_compiled.fastq"
        export comp_file2="${5}/${3}_2_compiled.fastq"
        cat "$run_file2" >> "$comp_file2"
        rm "$run_file2"
    fi

    cat "$run_file" >> "$comp_file"
    echo "Finished concatenating file!"
    rm "$run_file"
fi
