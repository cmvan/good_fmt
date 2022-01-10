"""
File: parse_metadata.py
-----------------
This program parses the metadata from the Li et al. 2016 study and organizes them for download.py
and bash_download.py.

The program assumes the tsv file provided contains the run accessions, sample accessions, aliases,
and library layouts.
"""
import pandas as pd 
import numpy as np


def tie_runs_to_samples(sorted_data, sample_accs):
    """
    This function takes the data sorted by sample accession and iterates through the list of sample
    accessions to obtain the list of run accessions associated with each sample accession. For each
    sample accession, it organizes the runs by single and paired reads before adding the subcategories
    of runs with the sample accession itself to dictionaries for single and paired reads respectively.

    It returns a dictionary which associates each sample accession to its single read run accessions
    and another dictionary which associates each sample accession to its paired read run accessions
    """
    single_samples = {}
    paired_samples = {}

    for sample_acc in sample_accs:
        single = []
        paired = []
        runs = sorted_data[sorted_data.sample_accession == sample_acc]['run_accession']

        for run in runs:
            ppp = (sorted_data[sorted_data.run_accession == run]['library_layout']).to_string()
            if ppp.find('SINGLE') > 0:
                single.append(run)
            else:
                paired.append(run)

        single_samples[sample_acc] = single
        paired_samples[sample_acc] = paired
    
    return single_samples, paired_samples


def tie_people_to_sample(sorted_data, people_list):
    """
    This function creates and returns a dictionary that aligns the prefixed aliases of patients and
    donors to alist of their respective sample accessions by using the prefixes of the aliases to 
    differentiate between the patients and donors. Particularly, it ensures that samples from the 
    same patient but different timepoints are associated together in a list.
    """
    people_sample_acc = {}

    for person in people_list:
        # Isolates the first seven characters of the alias. Ensures sample accessions of
        # patients from different time points are grouped together
        cut = person[:7]
        if cut == 'FAT_DON':
            cut = person[:10]
            if cut[9] == '-':
                cut = cut[:9]

        sample_acc = np.unique(sorted_data[sorted_data.alias == person]['sample_accession'])[0]
        if cut not in people_sample_acc:
            people_sample_acc[cut] = []
        people_sample_acc[cut].append(sample_acc)
    
    return people_sample_acc


def process_data(tsv, type):
    """
    This program reads the tsv file containing the run and sample accessions, aliases, and library
    layouts from the Li et al. 2016 study and sorts the data by sample accession in increasing 
    alphabetical-numerical order. It then creates two sample-to-run accession dictionaries for single 
    and paired reads respectively and a people to sample accession dictionary before returning the 3
    dictionaries.
    """
    data = pd.read_csv(tsv, sep='\t')
    sorted_data = data.sort_values(by='sample_accession')
    people = np.unique(sorted_data.alias.tolist())
    people_to_sample = tie_people_to_sample(sorted_data, people)
    sample_accs = np.unique(sorted_data.sample_accession.tolist())

    if type == 'samples':
        single_samples, paired_samples = tie_runs_to_samples(sorted_data, sample_accs)
        return single_samples, paired_samples, people_to_sample
    elif type == 'people':
        return people_to_sample
    else:
        raise TypeError