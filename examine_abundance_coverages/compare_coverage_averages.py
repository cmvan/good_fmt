"""
File: compare_coverage_averages.py
----------------------------------
Given the species profiles of two samples, this program writes a file comparing coverages of
species where at least one has nonzero coverage. It also includes the average coverage between
the two samples. The file is sorted by decreasing average coverage
"""
import config
from most_abundant_species import most_abundant
import operator
import os
import pandas as pd
import sys


"""
This function returns the sample identifier from the path to the sample's species profile
"""
def get_person_name(sample_path):
    name = sample_path.split('/')[1]
    if name.find('DON') < 0:
        name = name + sample_path.split('/')[2][-3:]
    return name


if __name__=='__main__':
    person1 = pd.read_csv(sys.argv[1], delimiter="\t")
    person2 = pd.read_csv(sys.argv[2], delimiter="\t")
    person1_name = get_person_name(sys.argv[1])
    person2_name = get_person_name(sys.argv[2])

    person1_species = most_abundant(sys.argv[1], 'file')
    person2_species = most_abundant(sys.argv[2], 'file')
    most_abundant_compiled = list(set(person1_species + list(person2_species)))

    species_coverage_matrix = []
    species_string = ""

    # Group coverages and average for a species and place into greater list 
    for species in most_abundant_compiled:
        person1_coverage = float(person1[person1.species_id == species]['coverage'])
        person2_coverage = float(person2[person2.species_id == species]['coverage'])
        avg_coverage = (person1_coverage + person2_coverage) / 2
        species_to_coverages = (species, person1_coverage, person2_coverage, avg_coverage)
        species_coverage_matrix.append(species_to_coverages)
        species_string += species + ","
    species_coverage_matrix = sorted(species_coverage_matrix, key=operator.itemgetter(3), reverse=True)

    # Species string derived and printed to complete exploratory MIDAS merges
    species_string = species_string[:-1]
    print(species_string)

    # Write species coverages and average to file
    output_file_name = person1_name + '_vs_' + person2_name + '_coverage_averages.txt'
    output_file_path = os.path.join(config.cov_avg_dir, output_file_name)
    output_file = open(output_file_path, "w")
    output_file.write('\t'.join(['species_id', person1_name, person2_name, 'avg_coverage']))
    for species, person1, person2, avg in species_coverage_matrix:
        output_file.write("\n")
        output_file.write("\t".join([species, str(person1), str(person2), str(avg)]))
    output_file.close()

    
    """
    Using midas species merge and marker coverages
    
    prevalence = pd.read_csv(sys.argv[1], delimiter="\t")
    prevalence = prevalence[prevalence.prevalence > 0]
    prevalence = prevalence.sort_values(by ='mean_coverage', ascending = 0)
    print(prevalence)
    """
