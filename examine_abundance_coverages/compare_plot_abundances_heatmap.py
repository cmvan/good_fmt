"""
File: plot_abundances_stacked_bar.py
------------------------------------
Given a species profile, this program creates a heatmap of the relative abundances of the
bacterial families for each of the samples in the provided directory with the largest abundances 
at the bottom of the stacked bar plot.

This program assumes the command line argument given is a directory containing all the species
profiles for a person across all timepoints and that the species profiles were produced from the
MIDAS metagenomic pipeline.

Specifically for graphing all species profiles for a given patient from the Li et al. 2016 study 
"""
import matplotlib.pyplot as plt
import numpy as np
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

    # Create matrix and labels
    family_labels = []
    abundance_matrix = []
    for i in range(30):
        family_labels.append(meta_abundances[i][0])
        
        # To graph abundances relative to each other
        abundances_sum = np.array(meta_abundances[i][1]).sum()
        proportions = [abundance / abundances_sum for abundance in meta_abundances[i][1]]
        abundance_matrix.append(proportions) 
        
        # To graph actual abundances
        #abundance_matrix.append(meta_abundances[i][1]) 
    
    # Graph heatmap
    hm = plt.imshow(abundance_matrix, cmap='coolwarm', interpolation="nearest")
    plt.colorbar(hm, label="Abundances")
    plt.xticks(ticks=np.arange(len(sample_names)),labels=sample_names, rotation=90, fontsize='small')
    plt.yticks(ticks=np.arange(len(family_labels)),labels=family_labels)
    person = dir.split(sep='/')[1]
    plt.title('Relative Abundances for ' + person + ' (15 Most Abundant Families)', pad=15)

    plt.show()