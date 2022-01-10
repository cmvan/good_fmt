"""
File: alt_reads_to_freqs.py
---------------------------

"""
import math
import pandas as pd
import os
import sys

gzipped_merged_snps_dir = sys.argv[1]
strain = sys.argv[2]
sample = gzipped_merged_snps_dir.split('/')[-1][:7]
gzipped_alt_reads = os.path.join(gzipped_merged_snps_dir, strain, 'alt_reads.txt.gz')
alt_reads_df = pd.read_csv(gzipped_alt_reads, delimiter='\t', compression='gzip')
freq_points = list(alt_reads_df.columns.values[1:-1])
num_points = len(freq_points)
num_entries = len(alt_reads_df)
alt_freqs = []

for entry in range(num_entries):
    freqs = []
    alt_read_ratios = list(alt_reads_df.iloc[entry, 1:-1])
    if alt_read_ratios.count('0,0') > math.ceil(num_points / 2):
        continue
    for ratio in alt_read_ratios:
        try:
            freq = float(ratio.split(',')[0])/float(ratio.split(',')[1])
            freqs.append(freq)
        except ZeroDivisionError:
            freqs.append(-1)
    freqs.insert(0, str(alt_reads_df.iloc[entry, 0]))
    alt_freqs.append(freqs)

alt_freqs_df = pd.DataFrame(alt_freqs, columns=alt_reads_df.columns.values[:-1])
alt_freqs_file_path = os.path.join(gzipped_merged_snps_dir, strain, 'alt_freqs.txt.gz')
alt_freqs_df.to_csv(alt_freqs_file_path, sep='\t', columns=alt_reads_df.columns.values[:-1], index=False, compression='gzip')