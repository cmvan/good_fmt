"""
File: engraft_compare_freqs.py
------------------------------
This program compares allele frequencies between the donor and the recipient at day 84 for a given
species to look for/confirm signs of engraftment. 
"""
import matplotlib.pyplot as plt
import os
import pandas as pd
import sys

"""
priv_snps_file = sys.argv[1]
sample_strain = priv_snps_file.split('/')[1] # TO CHANGE
priv_snps_df = pd.read_csv(priv_snps_file, delimiter="\t")

sample = sys.argv[2]
donor = sys.argv[3]
merged_snps_dir = os.path.join(os.curdir, 'midas_out/snps_merges', sample + '_merged_snps') #sys.argv[1]
engrafted_abundances_file = 'midas_out/species_merges/%s_merged_species/engraftment.txt' % sample
engraft_abund_df = pd.read_csv(engrafted_abundances_file, delimiter='\t')

species = [dir for dir in os.listdir(merged_snps_dir) if not dir.startswith('.') and not dir.contains('.')]
high_freq_species = [s for s in engraft_abund_df['species_id'] if s in species]
"""

############################################################################################################################

alt_freqs_file = 'midas_out/snps_merges/FAT_015_merged_snps/Collinsella_aerofaciens_61484/new_alt_freqs.txt.gz' #sys.argv[2]
alt_freqs_df = pd.read_csv(alt_freqs_file, delimiter='\t', compression='gzip')
donor = alt_freqs_df.columns[1]
sample = alt_freqs_file.split('/')[2][:7]
strain = alt_freqs_file.split('/')[3]
alt_freqs_df = alt_freqs_df[(alt_freqs_df[donor] != -1) & (alt_freqs_df['day_84'] != -1)]

fig, ax = plt.subplots(figsize=(15, 8))
ax.scatter(alt_freqs_df[donor], alt_freqs_df['day_84'], alpha=0.1)
ax.set_xlabel('%s Frequencies' % (donor), labelpad=15)
ax.set_ylabel('%s_84 %s Frequencies' % (sample, strain), labelpad=15)
ax.set_title('%s Frequencies Between %s and %s at day 84' % (strain, donor, sample), pad=15)

plt.show()