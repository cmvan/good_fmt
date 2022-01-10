"""
File: bash_download_exchange.py
-------------------------------
This program reads data from a tsv file containing the run accessions, aliases, sample accessions,
and library layouts from the Li et al. 2016 study and downloads the appropriate fastq files from the
NCBI Sequence Read Archive (SRA) using the SRA Tool-Kit

Accompanied by download_exchange.sh to download and gzip fastq files

This program assumes bash_download_exchange.py and download_exchange.sh are located in the current 
directory from which the program is called.

NOT GENERALIZED FOR PAIRED READS
"""
import subprocess
import sys
from parse_metadata import process_data


if __name__=='__main__':

    single_samples, paired_samples, people_to_sample = process_data(sys.argv[1], 'samples')

    if len(sys.argv) < 3:
        for person in people_to_sample:
            for sample in people_to_sample[person]:
                subprocess.call(['bash', 'download_exchange.sh', sys.argv[1], person, sample, 'single'])
                subprocess.call(['bash', 'download_exchange.sh', person, sample, 'paired'])
    
    else:
        sample_dict = single_samples
        if sys.argv[5] == 'paired':
            sample_dict = paired_samples
        
        for run in sample_dict[sys.argv[4]]:
            print(run)