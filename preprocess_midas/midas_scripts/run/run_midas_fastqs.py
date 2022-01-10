"""
File: run_midas_fastqs.py
-------------------------
This program completes MIDAS processing for all samples in the Li et al. 2016 study
Can complete MIDAS species, genes or snps processing

UPDATE: run_midas_by_increment.py is the updated version of this program, but modified for
        optimization with Sherlock

Time and Memory Specs for species:
About 42 minutes per timepoint * (50 for patients and 5 for donors) = 1 day 14.5 hours
1 M per sample * (10 samples + 1 combined for donors) = 11 M (very small)
"""
from parse_metadata import process_data
import subprocess
import sys

def call_midas(script, outdir, fastq_path):
    midas_args = ['bash', script, outdir, fastq_path]
    subprocess.call(midas_args)


if __name__=='__main__':
    people_to_sample = process_data(sys.argv[1], 'people')
    people = list(people_to_sample.keys())
    days = ['00', '02', '14', '42', '84'] # timepoints from Li et al. study

    for person in people:
        if len(people_to_sample[person]) == 5:
            for day in days:
                outdir = '/scratch/users/cmvan/midas_out/' + person + '/day_' + day
                fastq_path = '/scratch/users/cmvan/merged_reads/' + person + '/day_' + day + '.fastq.gz'
                call_midas(sys.argv[2], outdir, fastq_path)
        else:
            outdir = '/scratch/users/cmvan/midas_out/' + person
            fastq_path = '/scratch/users/cmvan/merged_reads/' + person + '/' + person + '.fastq.gz'
            call_midas(sys.argv[2], outdir, fastq_path)
        print("MIDAS processing for " + person + " is complete!")