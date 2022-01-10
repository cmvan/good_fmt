"""
File: process_snps_util.py
--------------------------

"""
import os
from math import floor, ceil
import stats_utils

def parse_coverage_distribution(merged_snps_dir): # species_name as arg
    """
    *** prevalence_filter = True as fn arg ***
    if prevalence_filter:
        full_str = ""
    else:
        full_str = "full_"
    """

    filename = "%s/coverage_distribution.txt" % (merged_snps_dir)

    if not os.path.isfile(filename):
        return [], []

    coverage_distribution_file = open(filename, "r")

    line = coverage_distribution_file.readline() # header
    samples = []
    sample_coverage_histograms = []
    for line in coverage_distribution_file:
        items = line.split()
        sample_coverage_histogram = {}
        for item in items[1:]:
            subitems = item.split(",")
            sample_coverage_histogram[float(subitems[0])] = float(subitems[1])
        sample_coverage_histograms.append(sample_coverage_histogram)
        samples.append(items[0])
    
    return sample_coverage_histograms, samples


def calculate_relative_depth_threshold_map(sample_coverage_histograms, samples, min_nonzero_median_coverage=3, lower_factor=0.3, upper_factor=3):
    # returns map of sample name: coverage threshold
    # essentially filtering out samples whose marker depth coverage
    # does not exceed the average coverage threshold
    
    num_samples = len(samples)
    depth_threshold_map = {}
    for i in range(num_samples):
        
        # Check if coverage distribution meets certain requirements
        is_bad_coverage_distribution = False
        
        # First check if passes median coverage requirement
        nonzero_median_coverage = stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histograms[i])
        if round(nonzero_median_coverage) < min_nonzero_median_coverage:
            is_bad_coverage_distribution=True
    
        # Passed median coverage requirement
        # Now check whether a significant number of sites fall between lower and upper factor. 
        lower_depth_threshold = floor(nonzero_median_coverage * lower_factor) - 0.5 # why is 0.5 being added/subtracted? NRG
        upper_depth_threshold = ceil(nonzero_median_coverage * upper_factor) + 0.5
    
        depths, depth_CDF = stats_utils.calculate_CDF_from_histogram(sample_coverage_histograms[i])
        # remove zeros
        if depths[0] < 0.5:
            depth_CDF -= depth_CDF[0]
            depth_CDF /= depth_CDF[-1]
        
        fraction_in_good_range = depth_CDF[(depths > lower_depth_threshold) * (depths < upper_depth_threshold)].sum()
    
        if fraction_in_good_range < 0.6: # where does 0.6 come from? NRG
            is_bad_coverage_distribution = True
            
        if is_bad_coverage_distribution:
            lower_depth_threshold = 1000000001
            upper_depth_threshold = 1000000001
        
        depth_threshold_map[samples[i]] = (lower_depth_threshold, upper_depth_threshold)
        
    return depth_threshold_map
