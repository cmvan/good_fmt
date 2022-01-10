"""
File: loop_compare_freqs.py
---------------------------

"""
import pandas as pd
import os
import subprocess
import sys

meta_merge_snps_dir = sys.argv[1]
sample = meta_merge_snps_dir.split('/')[-1][:7]
ifp_comparisons_path = os.path.join(meta_merge_snps_dir, sample + '_ifp_comparisons.txt')
ifp_df = pd.read_csv(ifp_comparisons_path, delimiter='\t')
full_ifp_df = ifp_df.query('num_timepoints == 5')
num_species = len(full_ifp_df)

for idx in range(num_species):
    species = full_ifp_df.iloc[idx, 0]
    snps_summary_path = os.path.join(meta_merge_snps_dir, species, 'snps_summary.txt.gz')
    snps_summary_df = pd.read_csv(snps_summary_path, delimiter='\t', compression='gzip')
    """
    print('\nStarting frequency comparison for %s' % species)
    subprocess.call(['python3', 'postprocess_snps/compare_frequencies.py', meta_merge_snps_dir, species])
    print('Completed comparing frequencies for %s' % species)
    """
    # Change as needed to cover other species which need processing but don't satisfy this criteria
    donor = [point for point in snps_summary_df.sample_id.values if 'FAT_DON' in point]
    if 'day_00' in snps_summary_df.sample_id.values and not 'day_84' in snps_summary_df.sample_id.values and donor:
        print('\nStarting frequency comparison for %s' % species)
        subprocess.call(['python3', 'postprocess_snps/compare_frequencies.py', meta_merge_snps_dir, species])
        print('Completed comparing frequencies for %s' % species)
    
print('\nCompleted comparing frequencies for species with full snps data in %s' % meta_merge_snps_dir)


