"""
File: postprocess_pipeline.py
-----------------------------

"""
import parse_midas_snps as parse_snps
import os
import sys


big_merge_snps_dir = sys.argv[1]
sample = big_merge_snps_dir.split('/')[-1][:7]
# species = sys.argv[2]
species_list = sorted([species for species in os.listdir(big_merge_snps_dir) if not species.startswith('.')])

for species in species_list:
    parse_snps.calculate_coverage_distribution(big_merge_snps_dir, species)
    parse_snps.pipe_midas_snps(big_merge_snps_dir, species)
    parse_snps.calculate_ifp(big_merge_snps_dir, species)
parse_snps.order_ifps(big_merge_snps_dir)
print("Completed pipeline for %s in %s\n\n" % (species, sample))