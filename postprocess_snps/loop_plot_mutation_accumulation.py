"""
File: loop_plot_mutation_accumulation.py
----------------------------------------

"""
import os
import subprocess
import sys

meta_merged_snps_dir = sys.argv[1]
merge_snps_dirs = [dir for dir in os.listdir(meta_merged_snps_dir) if not dir.startswith('.') 
                   and not dir.endswith('.txt')]

for species in merge_snps_dirs:
    sig_freqs_path = os.path.join(meta_merged_snps_dir, species, 'sig_freqs.txt')
    unsig_freqs_path = os.path.join(meta_merged_snps_dir, species, 'unsig_freqs.txt')
    if os.path.exists(sig_freqs_path) and os.path.exists(unsig_freqs_path):
        print('\nStarted plotting frequency trajectories for %s' % species)
        subprocess.call(['python3', 'postprocess_snps/plot_mutation_accumulation.py', sig_freqs_path, unsig_freqs_path])
        print('Completed plotting frequency trajectories for %s' % species)
print('Completed plotting frequency trajectories for all viable species in %s' % meta_merged_snps_dir)