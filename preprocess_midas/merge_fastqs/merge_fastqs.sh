#!/bin/bash

# Merges single and paired reads of sample at specific timepoint
# For use with the Sherlock Research Computer Cluster
# Usage: bash merge_fastqs.sh person sample_acc type day (if type = patient)

export single_read="/scratch/users/cmvan/outdir/single/${1}/${2}/${2}_compiled.fastq.gz"
export first_paired="/scratch/users/cmvan/outdir/paired/${1}/${2}/${2}_1_compiled.fastq.gz"
export second_paired="/scratch/users/cmvan/outdir/paired/${1}/${2}/${2}_2_compiled.fastq.gz"
export out_path="/scratch/users/cmvan/merged_reads/${1}/day_${4}.fastq.gz"

if [ "$3" = "donor" ]; then
    out_path="/scratch/users/cmvan/merged_reads/${1}/${1}.fastq.gz"
elif [ "$3" != "patient" ]; then
    echo "Improper fastq type!"
    exit 1
fi

if [[ -f "$out_path" ]]; then
    cat "$single_read" "$first_paired" "$second_paired" >> "$out_path"
else
    cat "$single_read" "$first_paired" "$second_paired" > "$out_path"
fi

echo "Merge complete!"