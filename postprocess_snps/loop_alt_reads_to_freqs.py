"""
File: loop_alt_reads_to_freqs.py
--------------------------------

"""
import os
import subprocess
import sys

meta_merged_snps_dir = sys.argv[1]
merge_snps_dirs = [dir for dir in os.listdir(meta_merged_snps_dir) if not dir.startswith('.') 
                   and not dir.endswith('.txt')]

for species in merge_snps_dirs:
    print('\nStarted converting alt reads for %s' % species)
    subprocess.call(['python3', 'postprocess_snps/alt_reads_to_freqs.py', meta_merged_snps_dir, species])
    print('Completed converting alt reads for %s' % species)
print('\nCompleted all conversions for %s' % meta_merged_snps_dir)