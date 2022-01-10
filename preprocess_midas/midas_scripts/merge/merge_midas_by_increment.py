"""
File: merge_midas_by_increment.py
---------------------------------
This program writes SBATCH job files to complete snps merges for each sample from the Li et al.
2016 study. After it submits each job file, Sherlock calls upon the program to complete the snps
merge on the given sample.
"""
import os
from parse_metadata import process_data
import subprocess
import sys


if os.path.exists(sys.argv[1]):
    people_to_sample = process_data(sys.argv[1], 'people')
    people = list(people_to_sample.keys())

    # Create directory for merge job files
    job_dir = os.path.join(os.getcwd(), 'merge_midas_snps_jobs')
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)

    # Write each job file and execute it (in Sherlock)
    for person in people:
        job_name = person + '_merge_snps_process'
        job_filename = person + '_merge_snps.job'
        job_file_path = os.path.join(job_dir, job_filename)

        job_file = open(job_file_path, 'w')
        job_file.write('#!/bin/bash\n')
        job_file.write('#SBATCH --nodes=2\n')
        job_file.write('#SBATCH --mem=2G\n')
        job_file.write('#SBATCH --job-name=%s\n' % job_name)
        job_file.write('#SBATCH --error=sherlock_out/%s_merge_snps.err\n' % job_name)
        job_file.write('#SBATCH --output=sherlock_out/%s_merge_snps.out\n' % job_name)
        job_file.write('#SBATCH --time=12:00:00\n')
        job_file.write('#SBATCH --qos=normal\n\n')

        job_file.write('ml py-numpy/1.14.3_py27\n')
        job_file.write('ml py-pandas/0.23.0_py27\n')
        job_file.write('ml biology\n')
        job_file.write('ml py-biopython/1.70_py27\n')
        job_file.write('ml py-pysam/0.14.1_py27\n')
        job_file.write('ml bowtie2\n\n')

        # $HOME/snps.sh hardcoded as program was written in mind of running midas_snps
        # Could generalize by modifying to include another command-line arg
        job_file.write('python $HOME/merge_midas_by_increment.py ' + person + ' ' +
                        sys.argv[1] + ' $HOME/merge_snps.sh')
        
        job_file.close()
        
        subprocess.call(['sbatch', '-c', '2', job_file_path])

else:
    person = sys.argv[1]
    people_to_sample = process_data(sys.argv[2], 'people')
    days = ['00', '02', '14', '42', '84'] # timepoints from Li et al. study

    # Completes snps merge if sample is a recipient sample
    if len(people_to_sample[person]) == 5:
        outdir = '/scratch/users/cmvan/midas_out/merges/' + person + '_merged_snps'
        os.mkdir(outdir)
        snps_list = ''
        for day in days:
            snps = '/scratch/users/cmvan/midas_out/' + person + '/day_' + day
            snps_list = snps_list + snps + ','
        snps_list = snps_list[:-1]
        subprocess.call(['bash', sys.argv[3], person, outdir, snps_list, '--all_sites', '--min_samples', '4'])
        print("MIDAS snps merge for " + person + " is complete!")