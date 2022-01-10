"""
File: plot_abundances_at_timepoint_piechart.py
----------------------------------------------
Given a species profile, this program creates a pie chart of the relative abundances of the
bacterial families of the sample with the largest abundances at the bottom of the stacked bar plot.

This program assumes the command line argument given is a species profile created from the MIDAS
metagenomic pipeline

Specifically for graphing individual species profiles at a timepoint from Li et al. 2016 study 
Reference for donut: https://medium.com/@kvnamipara/a-better-visualisation-of-pie-charts-by-matplotlib-935b7667d77f
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
    families = [family for family, abundances in sorted_family_pairs]
    abundances = [abundances for family, abundances in sorted_family_pairs]

    # Plot pie chart
    fig, ax = plt.subplots()
    pie_colors = dp.get_colors(num_abundances, pastel_factor=0.5)
    patches, texts = ax.pie(abundances, colors=pie_colors, normalize=False)

    # Draw circle for donut-shape
    center_circle = plt.Circle((0,0),0.40,fc='white')
    fig = plt.gcf()
    fig.gca().add_artist(center_circle)
    
    # Adding Legend, Titles, and Labels
    legend = plt.legend(patches, families, bbox_to_anchor=(1.1,1), loc="upper left", ncol=3,
                        fontsize='small', title="Gut Bacterial Families")
    legend._legend_box.sep = 10
    ax.set_title('Relative Abundances for ' + person + ' (By Family)', pad=10)
    fig.subplots_adjust(right=0.5)

    plt.show()