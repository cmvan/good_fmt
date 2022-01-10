"""
File: plot_abundances_at_timepoint_stacked_bar.py
-------------------------------------------------
Given a species profile, this program creates a stacked bar plot of the relative abundances of the
bacterial families of the sample with the largest abundances at the bottom of the stacked bar plot.

This program assumes the command line argument given is a species profile created from the MIDAS
metagenomic pipeline

Specifically for graphing individual species profiles at a timepoint from Li et al. 2016 study 
"""
from distinctipy import distinctipy as dp
import matplotlib.pyplot as plt
import pandas as pd
import sys
import plot_abundances_util as plot_util


if __name__=='__main__':

    # Obtain and organize species data into list of tuples with families and pooled abundances
    person, sorted_family_pairs = plot_util.parse_profile_data(sys.argv[1])
    num_abundances = len(sorted_family_pairs)

    # Graph families and relative abundances into a stacked bar chart
    fig, ax = plt.subplots()
    bar_colors = dp.get_colors(num_abundances, pastel_factor=0.7)
    bottom = 0
    for idx in range(num_abundances):
        ax.bar(person, sorted_family_pairs[idx][1], label=sorted_family_pairs[idx][0], bottom=bottom,
               color=bar_colors[idx])
        bottom = bottom + sorted_family_pairs[idx][1]

    # Adding Legend, Titles, and Labels
    legend = ax.legend(bbox_to_anchor=(1.2,1), loc="upper left", ncol=3, fontsize='medium', 
                       title="Gut Bacterial Families")
    legend._legend_box.sep = 10
    ax.set_ylabel('Relative Abundance', labelpad=12)
    ax.set_title('Relative Abundances for ' + person + ' (By Family)', pad=15)
    fig.subplots_adjust(right=0.4)
    plt.show()