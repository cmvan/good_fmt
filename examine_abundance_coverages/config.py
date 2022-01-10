"""
File: config.py
---------------
Directory of important variables used across programs
"""
#!/usr/bin/env python
import os

eac_cur_dir = os.getcwd() + '/examine_abundances_coverages'
cov_avg_dir = os.path.join(eac_cur_dir, 'compare_coverage_averages')
cov_diff_dir = os.path.join(eac_cur_dir, 'compare_coverage_differences')

