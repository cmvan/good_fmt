"""
File: private_snvs.py
---------------------

Pseudocode:
    Find all samples w/ strain of interest
    Merge all initial frequency columns into new dataframe
    Iterate through rows and checking to see which have high frequency only in one sample
        Can use a dictionary to tie samples to their own private marker snps
    Compare those private marker snps to human microbiome cohort 
        Where repolarization comes in (from minor/major to ref/alt)
        If site not listed or altfreq about 0, then private snv!
    Write file with private marker snvs for each sample

    Additional Notes:
        - Use dictionary with site_ids as keys and sample frequencies as values; fill with -1 for failed lookups
        - If fewer than half of the sites have valid frequencies, throw out
        - Consider making boolean numpy array for whether given site is a private snp, can pipe to snp frequencies dataframe

Questions:
    What to do with unique sites (i.e. sites not found in other samples)
"""
import numpy as np
import pandas as pd
import sys
import os


meta_merged_snps_dir = os.path.join(os.curdir, 'midas_out/snps_merges') #sys.argv[1]
strain = sys.argv[1]
alt_freqs_identifier = 'new_alt_freqs.txt.gz' #sys.argv[3]
strain_alt_freqs_paths = []
initial_freqs_df = pd.DataFrame()
priv_sites = dict()

sample_merges = [dir for dir in os.listdir(meta_merged_snps_dir) if not dir.startswith('.')]
for sample_dir in sample_merges:
    strain_dir = os.path.join(meta_merged_snps_dir, sample_dir, strain, alt_freqs_identifier)
    if os.path.exists(strain_dir):
        strain_alt_freqs_paths.append(strain_dir)
num_samples = len(strain_alt_freqs_paths)

for i in range(num_samples):
    alt_freqs_df = pd.read_csv(strain_alt_freqs_paths[i], delimiter='\t', compression='gzip')
    filtered_cols = [col for col in alt_freqs_df.columns if not col.startswith('FAT_')]
    alt_freqs_df = alt_freqs_df.iloc[:, [0, 2]]

    sample = strain_alt_freqs_paths[i].split('/')[3][:7]
    initial_rename = sample + '_' + str(alt_freqs_df.columns[1])
    alt_freqs_df = alt_freqs_df.rename(columns={alt_freqs_df.columns[1]:initial_rename})
    if initial_freqs_df.empty:
        initial_freqs_df = alt_freqs_df
    else:
        initial_freqs_df = pd.merge(initial_freqs_df, alt_freqs_df, on='site_id', how='outer')
num_sites = len(initial_freqs_df)

for i in range(num_samples):
    sample = initial_freqs_df.columns[i + 1][:7] + '_' + strain
    priv_sites[sample] = []

maxes = np.nanmax(np.array(initial_freqs_df.iloc[:, 1:]), axis=1)
max_idxs = np.nanargmax(np.array(initial_freqs_df.iloc[:, 1:]), axis=1)
for i in range(10000):
    max = maxes[i]
    max_idx = max_idxs[i]
    freqs = np.array(initial_freqs_df.iloc[i, 1:])
    np.delete(freqs, max_idx)
    if max >= 0.5 and pd.Series(freqs).isna().sum() < num_samples / 2:
        sample = initial_freqs_df.columns[max_idx + 1][:7] + '_' + strain
        site_id = initial_freqs_df.iloc[i, 0]
        priv_sites[sample].append((site_id, max))

for sample in priv_sites:
    sample_df = pd.DataFrame(priv_sites[sample], columns=['site_id', 'initial_freq'])
    file_name = sample + '.txt'
    test_path = os.path.join(os.getcwd(), 'midas_out', file_name)
    sample_df.to_csv(test_path, sep = '\t', index=False)