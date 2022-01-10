"""
File: compare_coverage_differences_scatterplot.py
-------------------------------------------------
Given the species profiles of two samples, this program writes a file comparing coverages of
species where at least one has nonzero coverage and includes the difference of coverages. The file
is sorted by decreasing coverage difference. The program also plots a scatterplot to visually
the coverages.
"""
import config
import matplotlib.pyplot as plt
import os
import pandas as pd
from pandas.core.frame import DataFrame
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
    coverage_directory = config.cov_diff_dir

    person1_species = person1.species_id
    person1_coverages = []
    person2_coverages = []
    person1_person2_meta = []

    # Group coverages, difference, species together and place into larger list
    for species in person1_species:
        person1_coverage = float(person1[person1.species_id == species]['coverage'])
        person2_coverage = float(person2[person2.species_id == species]['coverage'])

        if person1_coverage != 0 or person2_coverage != 0:
            person1_coverages.append(person1_coverage)
            person2_coverages.append(person2_coverage)
            species_meta = (species, person1_coverage, person2_coverage, 
                            person2_coverage - person1_coverage)
            person1_person2_meta.append(species_meta)

    # Create output file
    person1_person2_meta_df = DataFrame(person1_person2_meta,
                                        columns=['species_id', person1_name + "_coverages",
                                        person2_name + '_coverages', 'coverage_difference'])
    person1_person2_meta_df = person1_person2_meta_df.sort_values(by=['coverage_difference'], key=abs,
                                                                ascending=False)
    output_file_path = os.path.join(coverage_directory, person1_name + '_vs_' +
                                    person2_name + '_coverage_differences.txt')
    person1_person2_meta_df.to_csv(output_file_path, sep='\t', index=False)

    # Plot scatterplot comparing coverages
    fig, ax = plt.subplots()
    ax.scatter(person1_coverages, person2_coverages, color='r', edgecolors='black')
    ax.set_xlabel(person1_name + ' Coverages', labelpad=15)
    ax.set_ylabel(person2_name + ' Coverages', labelpad=15)
    ax.set_title(person1_name + ' Coverages vs. ' + person2_name + ' Coverages', pad=10)
    plt.show()