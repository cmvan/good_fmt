"""
File: plot_private_snvs.py
--------------------------

"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


# Do sample, strain as arguments???
priv_snps_file = 'midas_out/FAT_015_Collinsella_aerofaciens_61484.txt' 
sample_strain = priv_snps_file.split('/')[1] # TO CHANGE
priv_snps_df = pd.read_csv(priv_snps_file, delimiter='\t')

alt_freqs_file = 'midas_out/snps_merges/FAT_015_merged_snps/Collinsella_aerofaciens_61484/new_alt_freqs.txt.gz'
alt_freqs_df = pd.read_csv(alt_freqs_file, delimiter='\t', compression='gzip')

priv_freqs_df = alt_freqs_df[alt_freqs_df['site_id'].isin(priv_snps_df['site_id'])]
num_freqs = len(priv_freqs_df)
timepoints = np.array(priv_freqs_df.columns[2:])
index = np.where(timepoints.reshape(timepoints.size, 1) == np.array(timepoints))[1]
num_timepoints = len(timepoints)

fig, ax = plt.subplots()
for idx in range(num_freqs):
    freqs = priv_freqs_df.iloc[idx, 2:]
    good_idx = (freqs != -1)
    plt.plot(index[good_idx], freqs[good_idx], color='blue', linewidth=0.3)
    
ax.set_xticks(range(len(timepoints)))
ax.set_xticklabels(timepoints)
plt.xticks(fontsize=6)
plt.xlabel('Timepoints', size=8, labelpad=10)
plt.yticks(fontsize=6)
plt.ylabel('Allele Frequency', size=8, labelpad=10)
plt.title('Private SNP Trajectories for %s' % (sample_strain.split('.')[0]), size=10, pad=10)
plt.show()