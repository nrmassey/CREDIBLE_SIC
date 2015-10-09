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

from create_HadISST_SST_SIC_mapping import *
from create_SIC_from_mapping import create_coefficient_interpolants
import matplotlib.pyplot as plt
import numpy
from netcdf_file import *
from scipy import interpolate

#############################################################################

def get_CMIP5_SST_SIC_mapping_fname(rcp):
    out_dir = "/Users/Neil/Coding/CREDIBLE_output/output/"
    map_file = out_dir + rcp + "_2006_2100/cmip5_polyfit_" + rcp + "_2006_2100_1_anoms.nc"
    return map_file

#############################################################################

if __name__ == "__main__":
    lon = 175
    m = 0
    
    hadisst_sy = 1899
    hadisst_ey = 2010
    
    ref_start = 1986
    ref_end = 2005
    rn = 400
    hadisst_deg = 1
    RCP = "rcp45"
    
    hadisst_fit_file = get_HadISST_SST_SIC_mapping_fname(hadisst_sy, hadisst_ey, rn, hadisst_deg, anoms=True)
    cmip5_fit_file = get_CMIP5_SST_SIC_mapping_fname(RCP)
    
    fh = netcdf_file(hadisst_fit_file)
    hadisst_fit_data = fh.variables["polyfit"][:]
    fh.close()

    fh = netcdf_file(cmip5_fit_file)
    cmip5_fit_data = fh.variables["polyfit"][:]
    fh.close()
    
    print hadisst_fit_data.shape, cmip5_fit_data.shape
    
    years = get_year_intervals()
    sp0 = plt.subplot(111)
    sp1 = sp0.twinx()
    sy = years[0][0]
    ey = years[-1][1]
    mv = -1e30
    m = 0
    y_s = numpy.arange(sy,ey,10)
    
    for lat in range(160,180):
        coeff_splines = create_coefficient_interpolants(hadisst_fit_data, cmip5_fit_data, m, lat, lon, mv)
        y_new = numpy.arange(sy, ey+1)
#        sp0.plot(y_s, coeffs[0], 'k')
        sp0.plot(y_new, coeff_splines[0], 'r')

    plt.show()
