"""
File: plot_abundances_util.py
-----------------------------
Contains several helper functions for plotting abundances

Specifically for graphing species profiles for a given recipient from the Li et al. 2016 study 
"""
import pandas as pd
import operator


"""
Aggregates species and their abundances into families with pooled abundances and creates
a list of tuples mapping families to abundances sorted by decreasing pooled abundance.
"""
def map_abundances_to_families(species_list, relative_abundance):
    family_to_abundance = {}
    num_species = len(species_list)

    for idx in range(num_species):
        family = species_list[idx].split(sep='_', maxsplit=1)[0]
        if family not in family_to_abundance:
            family_to_abundance[family] = relative_abundance[idx]
        else:
            family_to_abundance[family] += relative_abundance[idx]
    family_pairs = [(family,family_to_abundance[family]) for family in family_to_abundance]
    sorted_family_pairs= list(sorted(family_pairs, key=operator.itemgetter(1), reverse=True))

    return sorted_family_pairs


"""
Obtains raw species profiling data from file and parses the patient's alias and a list of tuples
containing the families of the species and their pooled abundances if families bool is True
"""
def parse_profile_data(file):
    person_species_data = pd.read_csv(file, delimiter='\t')
    person_species_data = person_species_data[person_species_data.count_reads > 0]
    person = file.split(sep='/')[1] + file.split(sep='/')[2][-3:]
    species_list = person_species_data.species_id
    relative_abundance = person_species_data.relative_abundance
    sorted_family_pairs = []

    sorted_family_pairs = map_abundances_to_families(species_list, relative_abundance)
        
    return person, sorted_family_pairs


"""
Groups species from profile into families with their pooled abundances and creates a collective
list of families from across all samples in dir. Also provides a list of sample names to use as
labels for the graph
"""
def categorize_profile_pooled_data(dir, files):

    family_to_abundance_list = []
    collective_families = []
    sample_names = []

    for file in files:
        # Obtain data (species list and their relative abundances especially)
        person_species_data = pd.read_csv(dir + file + '/species/species_profile.txt', delimiter='\t')
        person_species_data = person_species_data[person_species_data.count_reads > 0]
        sample_names.append(file.split(sep="_p")[0])
        species_list = person_species_data.species_id
        relative_abundance = person_species_data.relative_abundance

        # Group families and their pooled abundances into dictionary sorted by decreasing relative abundance
        family_pairs = map_abundances_to_families(species_list, relative_abundance)
        family_list = [family_pair[0] for family_pair in family_pairs]

        # Append to grouped lists
        family_to_abundance_list.append(family_pairs)
        collective_families = list(set(collective_families + family_list))

    return family_to_abundance_list, collective_families, sample_names


"""
Creates a list of tuples tying families to each of the pooled abundances across the
samples at the different timepoints. Resulting structure used for graphing
"""
def make_graph_abundance_lists(family_to_abundance_list, collective_families):
    meta_abundances = []

    for family in collective_families:
        abundances_to_graph = []
        sum_abundances = 0
        for idx in range(len(family_to_abundance_list)):
            species_abundance_pair = [(x,y) for x,y in family_to_abundance_list[idx] if x == family]
            if species_abundance_pair:
                abundances_to_graph.append(species_abundance_pair[0][1])
                sum_abundances += species_abundance_pair[0][1]
            else:
                abundances_to_graph.append(0)
        meta_abundances.append((family, abundances_to_graph, sum_abundances))
    meta_abundances = sorted(meta_abundances, key=operator.itemgetter(2), reverse=True)

    return meta_abundances