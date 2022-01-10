"""
File: stats_utils.py
--------------------

"""
import numpy as np

####
#
# Calculates CDF from histogram
# histogram is map of value: counts
#
####
def calculate_CDF_from_histogram(histogram):

    xs = sorted(histogram.keys())
    ns = np.array([histogram[x] for x in xs]) * 1.0
    xs = np.array(xs)*1.0
    CDF = ns.cumsum()/ns.sum()
    return xs, CDF


####
#
# Calculates median from histogram
# histogram is map of value: counts
#
####
def calculate_nonzero_median_from_histogram(histogram):

    xs, CDF = calculate_CDF_from_histogram(histogram)
    
    if len(xs) < 2:
        return xs[0]
    
    if xs[0] < 0.5:
        if CDF[0] > 0.8:
            return xs[0]
        CDF -= CDF[0]
        CDF /= CDF[-1]
        
    median_idx = np.nonzero(CDF >= 0.5)[0][0]
    return xs[median_idx]