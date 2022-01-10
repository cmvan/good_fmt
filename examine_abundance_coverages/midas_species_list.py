"""
File: midas_species_list.py
---------------------------
This program prints a list of the most abundant species in a given sample from the Li et al.
2016 study
"""
from most_abundant_species import most_abundant
import sys

species = most_abundant(sys.argv[1], sys.argv[2])
num_species = len(species)
species_string = ""

# To get species list printed out
print(species)
print(len(species))

for i in range(num_species - 1):
    species_string += species[i] + ","
species_string += species[num_species - 1]

# To get format string for MIDAS_merge printed out
print(species_string)