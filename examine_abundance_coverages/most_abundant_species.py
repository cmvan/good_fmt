"""
File: most_abundant_species.py
------------------------------
Finds all species with coverage greater than 1 at at least 1 timepoint for donor/patient
Adaptable to find all species satisfying a certain coverage, relative abundance, number of reads

Specifically for parsing species profiles from the Li et al. 2016 study 
"""
import pandas as pd
import os


def most_abundant_at_timepoint(file):
    species_profile = pd.read_csv(file, delimiter='\t', usecols=['species_id', 'coverage'])

    # Can change this specification as needed
    return species_profile[species_profile.coverage > 5]['species_id']


def most_abundant(path, type):
    collective_abundant_species = []
    files = ''

    if type == 'dir':
        files = os.listdir(path)
    elif type == 'file':
        files = [path]
    else:
        raise TypeError

    for file in files:
        file_path = path
        if type == 'dir':
            if path[-1] != '/':
                path = path + '/'
            file_path = path + file
        abundant_species = most_abundant_at_timepoint(file_path)
        collective_abundant_species = list(set(collective_abundant_species + list(abundant_species)))

    return collective_abundant_species