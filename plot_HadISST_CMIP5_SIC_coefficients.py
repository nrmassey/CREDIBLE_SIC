#! /usr/bin/env python  
#############################################################################
#
# Program : plot_HadISST_CMIP5_SIC_coefficients.py
# Author  : Neil Massey
# Purpose : PLot the coefficients used in the reconstruction of the sea-ice 
#           concentration file (SIC) from the synthetic SST
#           To create the mapping run ./create_HadISST_SST_SIC_mapping.py
#                                  or ./create_CMIP5_SST_SIC_mapping.py
# Inputs  : same as ./create_HadISST_CMIP5_syn_SSTs.py
# Output  : in the output directory:
#           output_name
#            
# Date    : 22/07/15
#
#############################################################################

from create_HadISST_CMIP5_SIC_from_SST import load_sst_sic_mapping
from create_HadISST_SST_SIC_mapping import get_year_intervals, calc_polynomial
from create_SIC_from_mapping import create_coefficient_splines
import matplotlib.pyplot as plt
import numpy
from scipy import interpolate

if __name__ == "__main__":
    sst_sic_map = load_sst_sic_mapping("rcp45", True, 1)
    lon = 175
    m = 0
    
    years = get_year_intervals()
    sp0 = plt.subplot(111)
    sp1 = sp0.twinx()
    sy = years[0][0]
    ey = years[-1][1]
    mv = -1e30
    
    for lat in range(0,20):
        coeff_splines = create_coefficient_splines(sst_sic_map, m, lat, lon, mv)
        y_new = numpy.arange(sy, ey+1)
        sp0.plot(y_s, coeffs[0], 'k')
        sp0.plot(y_new, coeff_splines[0], 'r')

    plt.show()
