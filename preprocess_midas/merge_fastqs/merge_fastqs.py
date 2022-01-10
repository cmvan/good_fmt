"""
File: merge_fastqs.py
---------------------
This program individual sample runs into one larger file for each sample in the Li et al. 2016
study. For use with the Sherlock Research Computing Cluster

Time and Memory Specs (for Sherlock):
6 minutes * (10 samples and 1 for donors) = 66 minutes --> 1.5 hours
266 GB (same as outdir)
"""
from parse_metadata import process_data
import subprocess
import sys

people_to_sample = process_data(sys.argv[1], 'people')
people = list(people_to_sample.keys())
days = ['00', '02', '14', '42', '84'] # timepoints from Li et al. study
subprocess.call(['mkdir', '/scratch/users/cmvan/merged_reads'])

for person in people:
    type = 'patient'
    sample_accs = people_to_sample[person]
    num_accs = len(sample_accs)
    merged_dir = '/scratch/users/cmvan/merged_reads/' + person
    subprocess.call(['mkdir', merged_dir])

    if num_accs != 5:
        type = 'donor'
    for idx in range(num_accs):
        merge_args = ['bash', 'merge_fastqs.sh', person, sample_accs[idx], type, days[idx]]
        subprocess.call(merge_args)
    print('Merges for ' + person + ' are all complete!')
