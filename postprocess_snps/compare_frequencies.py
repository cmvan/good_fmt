"""
File: compare_frequencies.py
----------------------------

"""
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import sys

def clip(arr):
    return np.clip(arr, 0, 2)

gzipped_merged_snps_dir = sys.argv[1]
strain = sys.argv[2]
sample = gzipped_merged_snps_dir.split('/')[-1][:7]
gzipped_freqs = os.path.join(gzipped_merged_snps_dir, strain, 'new_alt_freqs.txt.gz')
freqs_df = pd.read_csv(gzipped_freqs, delimiter='\t', compression='gzip')
donor = freqs_df.columns.values[1]

like_freqs_df = freqs_df[((freqs_df[donor] > 0.8) & (freqs_df['day_00'] > 0.8)) | 
                         ((freqs_df[donor] < 0.2) & (freqs_df['day_00'] < 0.2))]
dif_freqs_df = freqs_df[((freqs_df[donor] > 0.8) & (freqs_df['day_00'] < 0.2)) | 
                        ((freqs_df[donor] < 0.2) & (freqs_df['day_00'] > 0.8))]
intermediate_freqs_df = freqs_df[(freqs_df[donor].between(0.2, 0.8)) | (freqs_df['day_00'].between(0.2, 0.8))]
unsig_freqs_df = pd.concat([like_freqs_df, intermediate_freqs_df], axis=0)
unsig_freqs_df.sort_values(by=['site_id'])
breakthrough_freqs_df = like_freqs_df[(like_freqs_df[donor] == 0) & (like_freqs_df['day_00'] == 0) & 
                                      (like_freqs_df['day_84'] > 0)]

dif_freqs_filepath = os.path.join(gzipped_merged_snps_dir, strain, 'sig_freqs.txt')
dif_freqs_df.to_csv(dif_freqs_filepath, sep='\t', columns=dif_freqs_df.columns.values, index=False)
unsig_freqs_filepath = os.path.join(gzipped_merged_snps_dir, strain, 'unsig_freqs.txt')
unsig_freqs_df.to_csv(unsig_freqs_filepath, sep='\t', columns=unsig_freqs_df.columns.values, index=False)
breakthrough_freqs_filepath = os.path.join(gzipped_merged_snps_dir, strain, 'breakthrough_freqs.txt')
breakthrough_freqs_df.to_csv(breakthrough_freqs_filepath, sep='\t', columns=breakthrough_freqs_df.columns.values,
                             index=False)

fig, ax = plt.subplots(figsize=(15, 8))
ax.scatter(clip(like_freqs_df['day_00']), clip(like_freqs_df[donor]), color='b', alpha=0.1)
ax.scatter(clip(dif_freqs_df['day_00']), clip(dif_freqs_df[donor]), color='r', alpha=0.1)
ax.scatter(clip(intermediate_freqs_df['day_00']), clip(intermediate_freqs_df[donor]), color='gray', alpha=0.1)

ax.set_xlabel('%s_00 %s Frequencies' % (sample, strain), labelpad=15)
ax.set_ylabel('%s %s Frequencies' % (donor, strain), labelpad=15)
ax.set_title('%s Frequencies Between %s and %s at day 0' % (strain, donor, sample), pad=15)

#plt.show()
figure_path = os.path.join(os.getcwd(), 'figures/compare_frequencies', sample, '%s_%s_%s_Frequencies_8_26_21' 
                           % (donor, sample, strain))
plt.savefig(figure_path, dpi=720)