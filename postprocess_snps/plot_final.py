"""
File: plot_final.py
-------------------

"""
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

def plot_freqs(df, query, color, return_df=False):
    freqs_df = df.query(query)
    num_freqs = len(freqs_df)
    for idx in range(num_freqs):
        freqs = freqs_df.iloc[idx, 2:]
        good_idx = (freqs != -1)
        plt.plot([time[-2:] for time in timepoints[good_idx]], freqs[good_idx], color=color, linewidth=0.3)
    if return_df:
        return freqs_df


sig_freqs_file = 'midas_out/snps_merges/FAT_012_merged_snps/Ruminococcus_bromii_62047/sig_freqs.txt'
sig_freqs_df = pd.read_csv(sig_freqs_file, delimiter='\t')

donor = sig_freqs_df.columns.values[1]
sample = sig_freqs_file.split('/')[2][:7]
strain = sig_freqs_file.split('/')[3]
timepoints = np.array(sig_freqs_df.columns.values[2:])

fig = plt.figure(figsize=(15, 8))
rows = 6
columns = 5
grid = plt.GridSpec(rows, columns, wspace = 0.7, hspace = 2)
plt.subplot(grid[:3, :4])

query = '%s > %s' % (donor, timepoints[0])
g_don_freqs_df = plot_freqs(sig_freqs_df, query, '#FC6F6F', True)
#plot_freqs(g_don_freqs_df, 'day_42 == 0', '#FFC300')
plot_freqs(g_don_freqs_df, 'day_42 == 1', '#8ED8F1')

plt.xticks(fontsize=6)
plt.xlabel('Days After FMT', size=8, labelpad=15)
plt.yticks(fontsize=6)
plt.ylabel('Frequency of Mutation', size=8, labelpad=15)
plt.title('Frequency Trajectories for %s in %s (Donor-Dominated Sites)' % (strain, sample), 
           size=10, pad=10, wrap=True)
plt.subplots_adjust(top=0.92)

#plt.show()
fig_path = os.path.join(os.getcwd(), 'figures/final', '%s_%s_Final_42_1' % (sample, strain))
plt.savefig(fig_path, dpi=720)