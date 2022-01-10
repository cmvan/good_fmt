"""
File: plot_abundances_stacked_bar.py
------------------------------------
Given a species profile, this program creates a stacked bar plot of the relative abundances of the
bacterial families for each of the samples in the provided directory with the largest abundances 
at the bottom of the stacked bar plot.

This program assumes the command line argument given is a directory containing all the species
profiles for a person across all timepoints and that the species profiles were produced from the
MIDAS metagenomic pipeline.

Specifically for graphing all species profiles for a given patient from the Li et al. 2016 study 
"""
from distinctipy import distinctipy as dp
import matplotlib.pyplot as plt
import os
import sys
import plot_abundances_util as plot_util


if __name__=='__main__':

    # Obtain and prepare data for graphing
    dir = sys.argv[1]
    if dir[-1] != '/':
        dir = dir + '/'
    files = sorted(os.listdir(sys.argv[1]))
    family_to_abundance_list, collective_families, sample_names = plot_util.categorize_profile_pooled_data(dir, files)
    meta_abundances = plot_util.make_graph_abundance_lists(family_to_abundance_list, collective_families)
    num_abundances = len(meta_abundances)
    bar_colors = dp.get_colors(num_abundances, pastel_factor=0.7)

    # Graph families and relative abundances into a stacked bar chart
    fig, ax = plt.subplots()
    bottom = len(files) * [0]
    for idx in range(num_abundances):
        ax.bar(sample_names, meta_abundances[idx][1], label=meta_abundances[idx][0], bottom=bottom,
               color=bar_colors[idx])
        bottom = [bottom[i] + meta_abundances[idx][1][i] for i in range(5)]

    # Adding Legend, Titles, and Labels
    legend = ax.legend(bbox_to_anchor=(1.2,1), loc="upper left", ncol=4, fontsize='x-small', 
                        title="Gut Bacterial Families")
    legend._legend_box.sep = 5
    ax.set_ylabel('Relative Abundance', labelpad=12)
    ax.set_xlabel('Timepoint', labelpad=12)
    person = dir.split(sep='/')[1]
    ax.set_title('Relative Abundances for ' + person + ' (By Family)', pad=15)
    fig.subplots_adjust(right=0.45)
    plt.show()