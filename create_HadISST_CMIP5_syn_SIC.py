#! /usr/bin/env python  
#############################################################################
#
# Program : create_HadISST_CMIP5_syn_SIC.py
# Author  : Neil Massey
# Purpose : Create a sea-ice concentration file (SIC) from an SST file and
#           a mapping file between SST and SIC.
#           To create the mappings run ./create_HadISST_SST_SIC_mapping.py
#                                  and ./create_CMIP5_SST_SIC_mapping.py
# Inputs  : input   : file of sst values
#           sic_fit : file of fitted sst to sic values
#           output  : output file name
# Output  : in the output directory:
#           output_name
#            
# Date    : 11/06/15
#
#############################################################################

import os, sys, getopt

import numpy

sys.path.append("../CREDIBLE_SST")
from cmip5_functions import load_data

sys.path.append("/Users/Neil/python_lib")
from create_HadISST_SST_SIC_mapping import get_HadISST_SST_SIC_mapping_fname, get_HadISST_monthly_ref_filenames
from sst_sic_mapping import *
from netcdf_file import *

# smoothing functions
import pyximport
pyximport.install(setup_args={'include_dirs':[numpy.get_include()]})
from window_smooth import window_smooth_3D
from fill_ice import *
from remove_isolated_ice import *

#############################################################################

def get_lsm():
    lsm_fname = "/Users/Neil/ClimateData/HadISST2/HadISST.2.1.0.0_sst_lsm.nc"
    fh = netcdf_file(lsm_fname)
    lsm = fh.variables["sst"][:]
    lsm_mv = fh.variables["sst"]._attributes["_FillValue"]
    fh.close()
    return lsm, lsm_mv

#############################################################################

def get_syn_SST_lon_lat_time_vars(input):
    fh = netcdf_file(input)
    lon_var = fh.variables["longitude"]
    lat_var = fh.variables["latitude"]
    time_var = fh.variables["time"]
    
    return lon_var, lat_var, time_var

#############################################################################

def create_SIC_from_mapping_file(input, output, hadisst_fit_fname, cmip5_fit_fname, rn):
    # load each mapping filename
    hadisst_mapping, mv = load_mapping(hadisst_fit_fname)
    cmip5_mapping, mv2  = load_mapping(cmip5_fit_fname)

    # load the hadisst monthly reference values
    sst_ref_fname, sic_ref_fname = get_HadISST_monthly_ref_filenames(1899, 2010, rn)
    sst_ref_data = load_data(sst_ref_fname, "sst")
    sic_ref_data = load_data(sic_ref_fname, "sic")
    
    # concatenate (along the t-axis) the two mapping files
    all_mapping = numpy.concatenate((hadisst_mapping, cmip5_mapping), axis=0)
    
    # load the sst file in
    sst_input = load_data(input, "sst")
    
    # truncate to 200 years worth of mapping data
    sub_mapping = all_mapping[5:]
    sub_sst = sst_input[0:-12]
    
    # reconstruct
    print "Constructing sea-ice"
    # remove the sst reference to produce the anomalies
    n_rpts = sub_sst.shape[0] / sst_ref_data.shape[0]
    sub_sst_anoms = sub_sst - numpy.tile(sst_ref_data, [n_rpts,1,1])

    syn_sic_anoms = calc_sic_from_sst(sub_sst_anoms, sub_mapping, mv, all_mapping.shape[1]-1)
    
    # add the reference back on
    syn_sic = syn_sic_anoms + numpy.tile(sic_ref_data, [n_rpts,1,1])
    
    print "Filling ice holes"
#    sic_filled = fill_ice(syn_sic, mv)
    
    print "Removing isolated ice"
#    sic_removed = remove_isolated_ice(sic_filled, mv)

    # smooth the ice with a 3x1 smoothing window
    weights = numpy.ones([1,3,1], 'f')
    print "Smoothing"
#    sic_smooth = window_smooth_3D(sic_removed, weights, mv, smooth_zero=True)
    # fix range
    sic_smooth = syn_sic
    
    sic_smooth[sic_smooth < mv] = mv
    sic_smooth[(sic_smooth < 0.0) & (sic_smooth != mv)] = 0.0
    sic_smooth[sic_smooth > 1.0] = 1.0
    
    # remove the block of sea ice in the Baltic sea caused by using the
    # 1986->2005 mean in months 04 to 11
    for m in range(4,12):
        sic_smooth[m::12,28:31,210:213] = 0.0
    
    # save the output
    lon_var, lat_var, time_var = get_syn_SST_lon_lat_time_vars(input)
    save_sic(output, sic_smooth, lon_var, lat_var, time_var, mv)
    print output

#############################################################################

if __name__ == "__main__":
    in_file = ""
    out_file = ""
    cmip5_fit_fname = ""
    hadisst_rn = 400
    hadisst_deg = 1
    hadisst_sy = 1850
    hadisst_ey = 2010
    opts, args = getopt.getopt(sys.argv[1:], 'i:f:o:d:',
                               ['input=', 'sic_fit=', 'output=', 'deg='])

    for opt, val in opts:
        if opt in ['--input', '-i']:
            in_file = val
        if opt in ['--output', '-o']:
            out_file = val
        if opt in ['--sic_fit', '-f']:
            cmip5_fit_fname = val
        if opt in ['--deg', '-d']:
            hadisst_deg = int(val)

    hadisst_fit_fname = get_HadISST_SST_SIC_mapping_fname(hadisst_sy, hadisst_ey, hadisst_rn, hadisst_deg)
    print hadisst_fit_fname
    print cmip5_fit_fname
    create_SIC_from_mapping_file(in_file, out_file, hadisst_fit_fname, cmip5_fit_fname, hadisst_rn)