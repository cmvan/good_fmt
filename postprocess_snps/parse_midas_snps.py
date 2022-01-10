"""
File: parse_midas_snps.py
-------------------------
Largely adapted from Roodger, Good et al. 2021 
"""
import numpy as np
import os
import pandas as pd
import process_snps_util as snps_util
import sys


"""
Creates a coverage_distribution.txt file in the species snp directory
 
In this file, rows are samples, columns are D,count pairs
samples are guaranteed to be same order as snps_depth.txtfile

Also creates a gene_coverage.txt file in the species snp directory
In this file, rows are genes, columns are samples, 
entries are avg coverage of that gene for that sample 
"""
def calculate_coverage_distribution(big_merged_snps_dir, species_name):
    merged_snps_dir = big_merged_snps_dir + '/' + species_name
    sample = big_merged_snps_dir.split('/')[-1][:7]

    sys.stderr.write("Calculating coverage distribution for %s in %s ...\n" % (species_name, sample))

    # These are the default MIDAS parameters. Used to ensure consistency
    # Only sites that pass this prevalence threshold in terms of having at least the prevalence_min_coverage can pass. 
    # These sites are stored in the sample_depth_histograms
    prevalence_threshold = 0.95 
    prevalence_min_coverage = 3

    allowed_variant_types = set(["1D","2D","3D","4D"]) # use all types of sites to include most information

    depth_file = pd.read_csv("%s/snps_depth.txt" % merged_snps_dir, delimiter="\t")
    file_entries = len(depth_file)
    info_file = pd.read_csv("%s/snps_info.txt" % merged_snps_dir, delimiter="\t")

    samples = list(depth_file.columns.values)[1:]
    sample_depth_histograms = {sample: {} for sample in samples} # stores only sites passng prevalence threshold.
    full_sample_depth_histograms = {sample: {} for sample in samples} # stores all sites. 

    gene_total_depths = {}
    gene_total_sites = {}

    num_sites_processed = 0

    for idx in range(file_entries):

        # make sure it is an allowed site
        variant_type = info_file['site_type'][idx]
        """
        For biallelic sites only
        
        snp_type = info_file['snp_type'][idx]
        if not variant_type in allowed_variant_types or snp_type != 'bi':
            continue
        """
        if not variant_type in allowed_variant_types:
            continue

        gene_name = info_file['gene_id'][idx]
        items = depth_file.iloc[idx, 1:]
        depths = [int(item) for item in items]

        # Manual prevalence filter
        if sum(depth >= prevalence_min_coverage for depth in depths) * 1.0 / len(depths) >= prevalence_threshold:
            # Add to genome-wide depth distribution
            for sample, D in zip(samples,depths):
                if D not in sample_depth_histograms[sample]:
                    sample_depth_histograms[sample][D]=0
                sample_depth_histograms[sample][D] += 1
        
        for sample,D in zip(samples,depths):
            if D not in full_sample_depth_histograms[sample]:
                full_sample_depth_histograms[sample][D]=0
            full_sample_depth_histograms[sample][D] += 1
        
        # Add to gene-specific avg depth
        if gene_name not in gene_total_depths:
            gene_total_depths[gene_name] = np.zeros(len(samples))*1.0
            gene_total_sites[gene_name] = 0.0
            
        gene_total_depths[gene_name]+=depths
        gene_total_sites[gene_name]+=1
        
        num_sites_processed+=1
            
        if num_sites_processed%100000==0:
            sys.stderr.write("Processed %dk sites!\n" % (num_sites_processed / 1000))

    # First write (filtered) genome-wide coverage distribution. This is filtered by prevalence
    output_file = open("%s/coverage_distribution.txt" % merged_snps_dir, "w")
    output_file.write("SampleID\tD,n(D) ...")
    for sample in samples:
        output_file.write("\n")
        sorted_histogram_keys = sorted(full_sample_depth_histograms[sample].keys())
        output_file.write("\t".join([sample]+["%d,%d" % (D, full_sample_depth_histograms[sample][D]) for D in sorted_histogram_keys]))
    output_file.close()

    # Write unfiltered genome-wide coverage distribution
    output_file = open("%s/full_coverage_distribution.txt" % merged_snps_dir, "w")
    output_file.write("SampleID\tD,n(D) ...")
    for sample in samples:
        output_file.write("\n")
        sorted_histogram_keys = sorted(full_sample_depth_histograms[sample].keys())
        output_file.write("\t".join([sample]+["%d,%d" % (D, full_sample_depth_histograms[sample][D]) for D in sorted_histogram_keys]))
    output_file.close()

    # Then write gene-specific coverages
    output_file = open("%s/gene_coverage.txt" % merged_snps_dir, "w")
    output_file.write("\t".join(["Gene"]+samples)) # Header line
    for gene_name in sorted(gene_total_depths.keys()):
        avg_depths = gene_total_depths[gene_name]/(gene_total_sites[gene_name]+(gene_total_sites[gene_name]==0))
        output_file.write("\n")
        output_file.write("\t".join([gene_name]+["%0.1f" % D for D in avg_depths]))
    output_file.close()

    # Done
    print("Finished calculating coverage distributions!\n")


def pipe_midas_snps(big_merged_snps_dir, species_name):
    merged_snps_dir = big_merged_snps_dir + '/' + species_name
    sample = big_merged_snps_dir.split('/')[-1][:7]

    # lower_factor = 0.3 is the default to be consistent with MIDAS gene presence criterion
    # upper factor = 3 is the default for (logarithmic) symmetry 
    min_samples = 1

    # Load genomic coverage distributions
    sample_coverage_histograms, sample_list = snps_util.parse_coverage_distribution(merged_snps_dir)

    # depth threshold map returns the lower and upper depth values that are 0.3 * median and 3 * median depth in the data. 
    # min_nonzero_median_coverage, lower_factor, upper_factor as args for below
    depth_threshold_map = snps_util.calculate_relative_depth_threshold_map(sample_coverage_histograms, sample_list)

    # Open MIDAS output files
    depth_file = pd.read_csv("%s/snps_depth.txt" % merged_snps_dir, delimiter="\t")
    file_entries = len(depth_file)
    info_file = pd.read_csv("%s/snps_info.txt" % merged_snps_dir, delimiter="\t")
    minor_freq_file = pd.read_csv("%s/snps_freq.txt" % merged_snps_dir, delimiter="\t")
    output_file = open("%s/alt_reads.txt" % merged_snps_dir, "w")

    # samples
    prevalence_threshold = min([min_samples * 1.0 / len(sample_list), 0.5])

    # create depth threshold vector from depth threshold map
    lower_depth_threshold_vector = []
    upper_depth_threshold_vector = []
    for sample in sample_list:
        lower_depth_threshold_vector.append(depth_threshold_map[sample][0])
        upper_depth_threshold_vector.append(depth_threshold_map[sample][1])
        
    lower_depth_threshold_vector = np.array(lower_depth_threshold_vector)
    upper_depth_threshold_vector = np.array(upper_depth_threshold_vector)

    # Figure out which samples passed our avg_depth_threshold
    # 1e09 comes from the calculate_relative_depth_threshold_map definition above, which is a code for a
    # bad sample. A bad sample has median depth less than 5 or greater than 0.6 fraction of the genome
    # is outside the acceptable range of good depths. 
    passed_samples = (lower_depth_threshold_vector < 1e09) 
    total_passed_samples = passed_samples.sum()

    # Let's focus on those from now on
    samples = [sample_list[idx] for idx in range(len(sample_list)) if passed_samples[idx]]
    lower_depth_threshold_vector = lower_depth_threshold_vector[passed_samples]
    upper_depth_threshold_vector = upper_depth_threshold_vector[passed_samples]

    # print header
    print_str = "\t".join(["site_id"] + samples + ["Total Ratios"])
    output_file.write(print_str)
    output_file.write("\n")

    # Only going to look at 1D, 2D, 3D, and 4D sites
    allowed_variant_types = set(['1D','2D','3D','4D'])

    num_sites_processed = 0

    for idx in range(file_entries):
        variant_type = info_file['site_type'][idx]
        """
        # For biallelic sites only

        snp_type = info_file['snp_type'][idx]

        # make sure it is a biallelic coding site
        if not variant_type in allowed_variant_types or snp_type != 'bi':
            continue
        """
        if not variant_type in allowed_variant_types:
            continue

        gene_name = info_file['gene_id'][idx]
        contig = info_file['ref_id'][idx]
        location = str(info_file['ref_pos'][idx])

        ref_allele = info_file['ref_allele'][idx]
        minor_allele = info_file['minor_allele'][idx]
        depths = np.array([float(item) for item in depth_file.iloc[idx, 1:]])
        ref_freqs = np.array([float(item) for item in minor_freq_file.iloc[idx, 1:]])

        if ref_allele != minor_allele:
            ref_freqs = np.array([1 - float(item) for item in minor_freq_file.iloc[idx, 1:]])
        
        depths = depths[passed_samples]
        ref_freqs = ref_freqs[passed_samples]
        
        refs = np.round(ref_freqs * depths)  
        alts = depths - refs

        passed_sites = (depths >= lower_depth_threshold_vector) * 1.0   
        passed_sites *= (depths <= upper_depth_threshold_vector)
        
        # make sure the site is prevalent in enough samples to count as "core"
        if ((passed_sites).sum() * 1.0 / total_passed_samples) < prevalence_threshold:
            continue
            #passed_sites *= 0
            
        refs = refs * passed_sites
        alts = alts * passed_sites
        depths = depths * passed_sites
        
        total_alts = alts.sum()
        total_refs = depths.sum()
        total_depths = total_alts + total_refs

        alts = np.append(alts, total_alts)
        refs = np.append(refs, total_refs)

        new_site_id_str = "|".join([contig, location, gene_name, variant_type])
        
        # format row into output file
        read_strs = ["%g,%g" % (A, A + R) for A,R in zip(alts, refs)]
        print_str = "\t".join([new_site_id_str] + read_strs)
        output_file.write(print_str)
        output_file.write("\n")
        
        num_sites_processed += 1
        if num_sites_processed % 100000 == 0:
            sys.stderr.write("%dk sites processed ...\n" % (num_sites_processed / 1000))   
            #if debug:
                #break

    output_file.close() 
    print("Finished piping snps!\n")


def calculate_ifp(big_merged_snps_dir, species_name):
    merged_snps_dir = big_merged_snps_dir + '/' + species_name
    sample = big_merged_snps_dir.split('/')[-1][:7]
    merged_snps_dir = os.path.join(big_merged_snps_dir, species_name)
    alt_reads = pd.read_csv(merged_snps_dir + '/alt_reads.txt', delimiter='\t')
    num_alt_read_entries = len(alt_reads)
    samples = alt_reads.columns.values[1:-1]
    ifp_fractions = []

    for sample in samples:
        num_ifp_sites = 0
        sites_processed = 0
        for entry in range(num_alt_read_entries):
            alt_read_pair = alt_reads[sample][entry]
            alt_count = float(alt_read_pair.split(",")[0])
            read_count = float(alt_read_pair.split(",")[1])
            if alt_count == 0 and read_count == 0:
                continue
            alt_read_percentage = alt_count / read_count
            if alt_read_percentage > 0.2 and alt_read_percentage < 0.8:
                num_ifp_sites += 1
            sites_processed += 1
        ifp_fractions.append((sample, str(num_ifp_sites), str(sites_processed),
                            str(num_ifp_sites / sites_processed)))

    output_file = open("%s/intermediate_freq_polymorphism.txt" % merged_snps_dir, "w") # GENERALIZE
    output_file.write("\t".join(["sample", "num_ifp_sites", "sites_processed", "ifp_fraction"]))
    for sample, num_ifp_sites, sites_processed, ifp_frac in ifp_fractions:
        output_file.write("\n")
        output_file.write("\t".join([sample, num_ifp_sites, sites_processed, ifp_frac]))
    output_file.close()

    print("Finished calculating intermediate frequency polymorphism!\n")


def order_ifps(big_merge_snps_dir):
    sample = big_merge_snps_dir.split('/')[-1][:7]
    species_list = sorted([species for species in os.listdir(big_merge_snps_dir) 
                        if not species.startswith('.') and not species.endswith('.txt')])
    species_to_ifp_averages = []

    for species in species_list:
        species_ifp_file = os.path.join(big_merge_snps_dir, species,'intermediate_freq_polymorphism.txt')
        species_ifp_data = pd.read_csv(species_ifp_file, delimiter='\t')
        species_ifps = np.array(species_ifp_data.ifp_fraction, dtype=np.float64)
        ifp_average = np.average(species_ifps)
        species_to_ifp_averages.append((species, len(species_ifps), ifp_average))

    ifp_comparisons_filepath = os.path.join(big_merge_snps_dir, sample + '_ifp_comparisons.txt')
    species_to_ifp_averages_df = pd.DataFrame(species_to_ifp_averages, columns=['species', 'num_timepoints',
                                                                                'ifp_average'])
    species_to_ifp_averages_df = species_to_ifp_averages_df.sort_values(by=['num_timepoints', 'ifp_average'], 
                                                                        ascending=[False, True])
    species_to_ifp_averages_df.to_csv(ifp_comparisons_filepath, sep='\t', float_format=np.float64,
                                    index=False)