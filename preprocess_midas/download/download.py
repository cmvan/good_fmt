"""
File: download.py
-----------------
This program reads data from a tsv file containing the run accessions, aliases, sample accessions,
and library layouts from the Li et al. 2016 study and downloads the appropriate fastq files from
the NCBI Sequence Read Archive (SRA) using the SRA Tool-Kit

The program assumes download.py is located in the current directory from which the program is called

USES OS.SYSTEM --> SECURITY CONCERN; BE AWARE OF CODE INJECTION
"""
import os
import sys
from parse_metadata import process_data


def download(people_dict, sample_dict, type):
    """
    This function downloads the appropriate fastq files from the SRA by order of patient in the
    people-to-sample dictionary. When all files are downloaded, everything is organized by type
    of read, patient, and sample accession, in that order. Each sample accession folder will
    include the original 3-4 fastq files from the run accessions under that sample and one 
    concatenated fastq file for the whole sample

    This function will make a directory called prefetch_loads to place the prefetchs and all
    runs will have the path outdir/$type/$person-from-people-dict/$sample_acc/$run.fastq.gz 
    (where $ denotes variables). The compiled sample fastq file will have the path
    outdir/$type/$person-from-people-dict/$sample_acc/$sample_acc_compiled.fastq.gz 
    """
    for person in people_dict:
        for sample_acc in people_dict[person]:
            out = 'outdir/' + type + '/' + person + '/' + sample_acc
            for run in sample_dict[sample_acc]:
                path = 'prefetch_loads/' + type + '/' + run
                prefetch_cmd = "prefetch " + run + " -o " + path
                fasterq_cmd = "fasterq-dump " + run + " --outdir " + out
                cat_cmd = 'cat ' + out + '/' + run + '.fastq' + ' >> ' + out + '/' + sample_acc + '_compiled.fastq'
                os.system(prefetch_cmd)
                os.system(fasterq_cmd)
                os.system(cat_cmd) 
            os.system("gzip -r " + out)


if __name__=='__main__':

    single_samples, paired_samples, people_to_sample = process_data(sys.argv[1], 'samples')

    download(people_to_sample, single_samples, 'single')
    download(people_to_sample, paired_samples, 'paired')