"""
File: merge_midas_fastqs.py
---------------------------
This program completes snps merges for all recipient samples in the Li et al. 2016 study

NOTE: merge_midas_by_increment.py is the updated version of this program, but modified for
        optimization with Sherlock
"""
import os
from parse_metadata import process_data
import subprocess
import sys

people_to_sample = process_data(sys.argv[1], 'people')
people = list(people_to_sample.keys())
days = ['00', '02', '14', '42', '84'] # timepoints from Li et al. study
merge_dir = '/scratch/users/cmvan/midas_out/merges/'
if not os.path.isdir(merge_dir):
    os.mkdir(merge_dir)

for person in people:
    if person.find('DON') != -1:
        continue
    outdir = '/scratch/users/cmvan/midas_out/merges/' + person + '_merged_snps'
    os.mkdir(outdir)
    snps_list = ''
    for day in days:
        snps = '/scratch/users/cmvan/midas_out/' + person + '/day_' + day
        snps_list = snps_list + snps + ','
    snps_list = snps_list[:-1]
    subprocess.call(['bash', sys.argv[2], person, outdir, snps_list])
        