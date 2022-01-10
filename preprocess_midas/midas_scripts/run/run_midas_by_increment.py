"""
File: run_midas_snps.py
-----------------------
This program writes SBATCH job files to complete MIDAS processing for each sample from the Li et al.
2016 study. After it submits each job file, Sherlock calls upon the program again, where it then
calls for the MIDAS processing on the given sample.

Can complete MIDAS species, genes, or snps processing; but snps processing currently hard-coded, as
program was written with that processing in mind

Time and Memory Specs for snps: (With 1 CPU)
FAT_DON_8 - 131m5.940s (real), 2.57GB max memory, 36 genomes
FAT_012_00 - 138m13.472s (real), 2.9GB max memory, 40 genomes
FAT_012_02 - 151m10.684s (real), 2.7GB max memory, 43 genomes
FAT_012_14 - 129m41.292s (real), 2.54GB max memory, 41 genomes
FAT_012_42 - 73m44.060s (real), 1.71GB max memory, 26 genomes
FAT_012_84 - 119m50.740s (real), 1.73GB max memory, 27 genomes
"""
import os
from parse_metadata import process_data
import subprocess
import sys


"""
This function calls the given midas script with the given output directory and path to fastqs files
"""
def call_midas(script, outdir, fastq_path):
    midas_args = ['bash', script, outdir, fastq_path]
    subprocess.call(midas_args)


if os.path.exists(sys.argv[1]):
    people_to_sample = process_data(sys.argv[1], 'people')
    people = list(people_to_sample.keys())

    # Create directory for MIDAS processing job files
    job_dir = os.path.join(os.getcwd(), 'midas_snps_jobs')
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)

    # Write each job file and execute it (in Sherlock)
    for person in people:
        job_name = person + '_snps_process'
        job_filename = person + '_snps.job'
        job_file_path = os.path.join(job_dir, job_filename)

        job_file = open(job_file_path, 'w')
        job_file.write('#!/bin/bash\n')
        job_file.write('#SBATCH --nodes=2\n')
        job_file.write('#SBATCH --mem=5G\n')
        job_file.write('#SBATCH --job-name=%s\n' % job_name)
        job_file.write('#SBATCH --error=sherlock_out/%s_snps.err\n' % job_name)
        job_file.write('#SBATCH --output=sherlock_out/%s_snps.out\n' % job_name)
        job_file.write('#SBATCH --time=18:00:00\n')
        job_file.write('#SBATCH --qos=normal\n\n')

        job_file.write('ml py-numpy/1.14.3_py27\n')
        job_file.write('ml py-pandas/0.23.0_py27\n')
        job_file.write('ml biology\n')
        job_file.write('ml py-biopython/1.70_py27\n')
        job_file.write('ml py-pysam/0.14.1_py27\n')
        job_file.write('ml bowtie2\n\n')

        # $HOME/snps.sh hardcoded as program was written in mind of running midas_snps
        # Can generalize using another command-line arg
        job_file.write('time python $HOME/run_midas_by_increment.py ' + person + ' ' +
                        sys.argv[1] + ' $HOME/snps.sh')
        job_file.close()
        
        subprocess.call(['sbatch', '-c', '2', job_file_path])

else:
    person = sys.argv[1]
    people_to_sample = process_data(sys.argv[2], 'people')
    days = ['00', '02', '14', '42', '84'] # timepoints from Li et al. study

    # If recipient sample, complete snps processing for each timepoint
    if len(people_to_sample[person]) == 5:
        for day in days:
            outdir = '/scratch/users/cmvan/midas_out/' + person + '/day_' + day
            fastq_path = '/scratch/users/cmvan/merged_reads/' + person + '/day_' + day + '.fastq.gz'
            call_midas(sys.argv[3], outdir, fastq_path)
    # Else complete one bout of snps processing for the donor sample
    elif len(people_to_sample[person]) == 3 or len(people_to_sample[person]) == 1:
        outdir = '/scratch/users/cmvan/midas_out/' + person
        fastq_path = '/scratch/users/cmvan/merged_reads/' + person + '/' + person + '.fastq.gz'
        call_midas(sys.argv[3], outdir, fastq_path)
        
    print("MIDAS processing for " + person + " is complete!")