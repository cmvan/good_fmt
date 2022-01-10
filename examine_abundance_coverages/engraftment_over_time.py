"""
File: engraftment_over_time.py
-------------------------------------
This program graphs changes in relative abundance for species which reflect a high donor abundance
(> 1e-03) and a low recipient abundance at day 00 (< 2e-04) in a given sample from the 
Li et al. 2016 study. It also writes a file containing the relative abundances for all such
species in the donor and across all recipient timepoints

The program takes the relative abundance file for the donor and the recipient as command-line 
arguments.
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys


def get_sample_name(rel_abundances_filepath):
    sample_merged_species_dir = rel_abundances_filepath.split('/')[-2]
    if sample_merged_species_dir.find('DON') != -1:
        donor_identifer = sample_merged_species_dir.split('_')[2]
        return 'FAT_DON_' + donor_identifer
    else:
        return sample_merged_species_dir[:7]


if __name__=='__main__':
    donor_rel_abundances_file = sys.argv[1]
    donor_name = get_sample_name(donor_rel_abundances_file)
    donor_abundance_data = pd.read_csv(donor_rel_abundances_file, delimiter='\t', compression='gzip')
    recipient_rel_abundances_file = sys.argv[2]
    recipient_name = get_sample_name(recipient_rel_abundances_file)
    recipient_abundance_data = pd.read_csv(recipient_rel_abundances_file, delimiter='\t', compression='gzip')
    num_species = len(donor_abundance_data) # num of species for donor and recipient always equal
    timepoints = list(recipient_abundance_data.columns)[1:]
    num_timepoints = len(timepoints)
    meta_donor_rec_abundances = []

    # Group relative abundance for donor, all recipient timepoints together for each abundant species
    for idx in range(int(num_species)):
        rec_abundances = list(recipient_abundance_data.iloc[idx])
        species = rec_abundances.pop(0)
        donor_abundance = float(donor_abundance_data[donor_abundance_data.species_id == species][donor_name]) # timepoints[0] with autologous treatment
        
        # Determine whether species is adequately abundant; if so, then group
        if donor_abundance > 1e-03 and rec_abundances[0] < 2e-04 and donor_abundance > 10 * rec_abundances[0]:
            clipped_rec_abundances = np.zeros(num_timepoints)
            np.clip(rec_abundances, 10e-5, 1, out=clipped_rec_abundances)
            meta_donor_rec_abundances.append((species, donor_abundance, clipped_rec_abundances))
    meta_donor_rec_abundances = sorted(meta_donor_rec_abundances, key=lambda x: abs(x[1] - x[2][0]), reverse=True)
    species = [abundance[0] for abundance in meta_donor_rec_abundances]
    num_species = len(species)


    # Write donor, recipient relative abundances to output file for all abundant species
    output_file_path = os.path.join(os.getcwd(), 'midas_out/species_merges', '%s_merged_species' % 
                                    recipient_name, 'engraftment.txt')
    output_file = open(output_file_path, 'w')
    headers = ['species_id', donor_name]
    headers.extend(timepoints)
    headers.extend(['initial_abundance_difference'])
    output_file.write('\t'.join(headers))
    output_file.write('\n')
    for idx in range(num_species):
        inputs = [meta_donor_rec_abundances[idx][0], str(meta_donor_rec_abundances[idx][1])]
        inputs.extend([str(abundance) for abundance in meta_donor_rec_abundances[idx][2]])
        inputs.extend([str(meta_donor_rec_abundances[idx][1] - meta_donor_rec_abundances[idx][2][0])])
        output_file.write('\t'.join(inputs))
        output_file.write('\n')
    output_file.close()

    # Plot changes in relative abundances
    fig, ax = plt.subplots(figsize=(15, 8)) #30, 10 for wide; 15, 8 for proportioned
    colors = ['black', 'red', 'blue', 'orange', 'purple']
    clipped_donor = np.zeros(num_species)
    xticks = []
    xticklabels = []
    np.clip([abundance[1] for abundance in meta_donor_rec_abundances], 10e-5, 1, out=clipped_donor)

    current_position = 0
    for idx in range(num_species):
        positions = []
        donor_abundance = clipped_donor[idx]
        initial_rec_abundance = meta_donor_rec_abundances[idx][2][0]
        xticks.append(current_position)
        xticklabels.append(species)
        ax.plot([current_position,current_position],[donor_abundance, initial_rec_abundance],'k-')

        current_position += 1
        
        for time_idx in range(num_timepoints):
            if time_idx == 0:
                continue
            positions.append(current_position) 
            current_position += 1
            
        ax.plot(positions, meta_donor_rec_abundances[idx][2][1:],'r-')
    
    all_positions = list(np.arange(stop=num_species * num_timepoints))
    ax.plot([i for i in all_positions if i % num_timepoints == 0], clipped_donor, marker='^', markersize=5, 
            ls='', color='k', label='Donor')
    for idx in range(num_timepoints):
        marker = 'v' if idx == 0 else '.'
        ax.plot([i for i in all_positions if i % num_timepoints == idx], [abundances[2][idx] for abundances in meta_donor_rec_abundances], 
                marker=marker, markersize=5, ls='', color=colors[idx], label=timepoints[idx])

    # Set legend, labels, title
    ax.set_xticks(xticks)
    legend = ax.legend(bbox_to_anchor=(1.01,1), title='Timepoints')        
    ax.set_xlabel('Species', labelpad=15)
    ax.set_xticklabels(species, size=8, rotation=45, ha='right', rotation_mode='anchor')
    ax.set_ylabel('Relative Abundances', labelpad=15)
    ax.set_yscale('log')
    ax.set_title(donor_name + ' vs. ' + recipient_name + ' Relative Abundances Across All Timepoints')
    fig.subplots_adjust(bottom=0.25)

    #plt.show()
    # Save the figure
    fig_path = os.path.join(os.getcwd(), 'figures', 'potential_engraftment', 
                            '%s_vs_%s_Potential_Engraftment' % (donor_name, recipient_name))
    plt.savefig(fig_path, dpi=720)