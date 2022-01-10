"""
File: plot_mutation_accumulation.py
-----------------------------------

MAY NEED TO FIX LABEL ORDERING, CHECK ON IT
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from textwrap import wrap


def plot_sig_freqs(df, hg, query, cluster_label, color, return_df=False):
    freqs_df = df.query(query)
    num_freqs = len(freqs_df)
    hg.append((cluster_label, num_freqs))
    for idx in range(num_freqs):
        freqs = freqs_df.iloc[idx, 2:]
        good_idx = (freqs != -1)
        plt.plot(timepoints[good_idx], freqs[good_idx], color=color, linewidth=0.3)
    if return_df:
        return freqs_df


sig_freqs_file = sys.argv[1] #'midas_out/snps_merges/FAT_012_merged_snps/Ruminococcus_bromii_62047/sig_freqs.txt' 
sig_freqs_df = pd.read_csv(sig_freqs_file, delimiter='\t')

donor = sig_freqs_df.columns.values[1]
sample = sig_freqs_file.split('/')[2][:7]
strain = sig_freqs_file.split('/')[3]
timepoints = np.array(sig_freqs_df.columns.values[2:])
num_timepoints = len(timepoints)
g_don_hg = []
g_rec_hg = []
colors = ['#FC6F6F', '#FFC300', '#8ED8F1']

# For graphing components of plot
#colors = ['#FC6F6F']
#colors = ['#FC6F6F', '#FFC300']
#colors = ['#FC6F6F', '#8ED8F1']


fig = plt.figure(figsize=(15, 8))
rows = 6
columns = 5
grid = plt.GridSpec(rows, columns, wspace = 0.7, hspace = 2)
plt.subplot(grid[:2, :3])

query = '%s > %s' % (donor, timepoints[0])
g_don_freqs_df = plot_sig_freqs(sig_freqs_df, g_don_hg, query, 'donor_freqs', '#FC6F6F', True)
plot_sig_freqs(g_don_freqs_df, g_don_hg, 'day_42 == 0', 'day_42_0', '#FFC300')
plot_sig_freqs(g_don_freqs_df, g_don_hg, 'day_42 == 1', 'day_42_1', '#8ED8F1')

plt.xticks(fontsize=6)
plt.xlabel('Timepoints', size=8, labelpad=10)
plt.yticks(fontsize=6)
plt.ylabel('Allele Frequency', size=8, labelpad=10)
plt.title('Minor Allele Frequency Changes in %s for %s (Donor > Recipient)' % (sample, strain), 
           size=10, pad=10, wrap=True)
plt.subplots_adjust(top=0.92)

plt.subplot(grid[:2, 3:])
plt.bar([df[0] for df in g_don_hg], [df[1] for df in g_don_hg], color=colors)
plt.yticks(fontsize=6)
plt.ylabel('# of Frequencies', size=8, labelpad=10)
plt.title('Higher Donor Cluster Distributions', size=10, pad=5)
plt.subplots_adjust(top=0.92)

for idx in range(num_timepoints):
    plt.subplot(grid[2, idx])
    good_idx = (g_don_freqs_df[timepoints[idx]] != -1)
    plt.hist(g_don_freqs_df[timepoints[idx]][good_idx], bins=50)
    plt.xticks(fontsize=8)
    plt.xlabel('Frequencies', size=7, labelpad=5)
    plt.yticks(fontsize=6)
    plt.ylabel('# of Frequencies', size=7, labelpad=5)
    plt.yscale('log')
    title = 'Frequency Distributions at %s' % timepoints[idx]
    plt.title('\n'.join(wrap(title, 23)), size=9, wrap=True)
    plt.subplots_adjust(left=0.07, bottom=0.07, right=0.95)


plt.subplot(grid[3:5, :3])

query = '%s < %s' % (donor, timepoints[0])
g_rec_freqs_df = plot_sig_freqs(sig_freqs_df, g_rec_hg, query, 'rec_freqs', '#FC6F6F', True)
plot_sig_freqs(g_rec_freqs_df, g_rec_hg, 'day_42 == 0', 'day_42_0', '#FFC300')
plot_sig_freqs(g_rec_freqs_df, g_rec_hg, 'day_42 == 1', 'day_42_1', '#8ED8F1')

plt.xticks(fontsize=6)
plt.xlabel('Timepoints', size=8, labelpad=10)
plt.yticks(fontsize=6)
plt.ylabel('Allele Frequency', size=8, labelpad=10)
plt.title('Minor Allele Frequency Changes in %s for %s (Donor < Recipient)' % (sample, strain), 
           size=10, pad=10, wrap=True)
plt.subplots_adjust(top=0.92)

plt.subplot(grid[3:5, 3:])
plt.bar([df[0] for df in g_rec_hg], [df[1] for df in g_rec_hg], color=colors)
plt.yticks(fontsize=6)
plt.ylabel('# of Frequencies', size=8, labelpad=10)
plt.title('Higher Recipient Cluster Distributions', size=10, pad=5)
plt.subplots_adjust(top=0.92)

for idx in range(num_timepoints):
    plt.subplot(grid[5, idx])
    good_idx = (g_rec_freqs_df[timepoints[idx]] != -1)
    plt.hist(g_rec_freqs_df[timepoints[idx]][good_idx], bins=50)
    plt.xticks(fontsize=8)
    plt.xlabel('Frequencies', size=7, labelpad=5)
    plt.yticks(fontsize=6)
    plt.ylabel('# of Frequencies', size=7, labelpad=5)
    plt.yscale('log')
    title = 'Frequency Distributions at %s' % timepoints[idx]
    plt.title('\n'.join(wrap(title, 23)), size=9, wrap=True)
    plt.subplots_adjust(left=0.07, bottom=0.07, right=0.95)

plt.show()
#fig_path = os.path.join(os.getcwd(), 'figures/mutation_frequencies', sample, '%s_%s_Draft' % (sample, strain)) # Change filename as needed
#plt.savefig(fig_path, dpi=720)