#! /usr/bin/env python  
#############################################################################
#
# Program : create_MaRIUS_future_SIC.py
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
from create_HadISST_sst_anoms import get_HadISST_monthly_reference_fname, get_HadISST_input_filename
from create_HadISST_CMIP5_syn_SSTs import save_3d_file
from create_CMIP5_sst_anoms import get_start_end_periods
from sst_sic_mapping import *
from calc_sea_ice_extent import calc_grid_areas
from netcdf_file import *
from cdo import *
from create_HadISST_CMIP5_syn_SIC import *
from create_MaRIUS_future_SSTs import *

# smoothing functions
import pyximport
pyximport.install(setup_args={'include_dirs':[numpy.get_include()]})
from window_smooth import window_smooth_3D
from fill_ice import *
from remove_isolated_ice import *

#############################################################################

def create_Ma_SIC(input, output, cmip5_fit_fname, sy, ey, rn):
    # load each mapping filename
    cmip5_mapping, mv = load_mapping(cmip5_fit_fname)

    # load the hadisst monthly reference values
    sst_ref_fname, sic_ref_fname = get_HadISST_monthly_ref_filenames(1899, 2010, rn)
    sic_ref_data = load_data(sic_ref_fname, "sic")
    sst_ref_fname = get_HadISST_monthly_reference_fname(1899, 2010, 1986, 2005, rn)
    sst_ref_data = load_data(sst_ref_fname, "sst")
    lon_var, lat_var, time_var = get_syn_SST_lon_lat_time_vars(input)
    
    # load the sst file in
    sst_input = load_data(input, "sst")
    
    # truncate to sy->ey years worth of mapping data
    sub_mapping = cmip5_mapping[(sy-2015)/10:(ey-2015)/10]
       
    # reconstruct
    print "Constructing sea-ice"
    # remove the sst reference to produce the anomalies
    n_rpts = sst_input.shape[0] / sst_ref_data.shape[0]
    sub_sst_anoms = sst_input - numpy.tile(sst_ref_data, [n_rpts,1,1])

    # create the SIC from the SST
    syn_sic_anoms = calc_sic_from_sst(sub_sst_anoms, sub_mapping, mv, sub_mapping.shape[1]-1)
    # add the reference back on
    n_rpts = syn_sic_anoms.shape[0] / sic_ref_data.shape[0]
    syn_sic = syn_sic_anoms + numpy.tile(sic_ref_data, [n_rpts,1,1])
    syn_sic[syn_sic < -1] = mv
    # fix range        
    #
    print "Filling ice holes"
    sic_filled = fill_ice(syn_sic, mv)
    
    print "Removing isolated ice"
    sic_removed = remove_isolated_ice(sic_filled, mv)

    # smooth the ice with a 3x1 smoothing window
    weights = numpy.ones([1,3,1], 'f')
    print "Smoothing"
    syn_sic[syn_sic==mv] = 0.0
    sic_removed = syn_sic
    sic_smooth = window_smooth_3D(sic_removed, weights, mv, smooth_zero=True)
    sic_smooth[(sic_smooth > 1.0) & (sic_smooth != mv)] = 1.0
    for m in [0,1,2,3,4,10,11]:
        syn_sic[m::12][syn_sic[m::12] != mv] = numpy.abs(syn_sic[m::12][syn_sic[m::12] != mv])
    sic_smooth[(sic_smooth < 0.0) & (sic_smooth != mv)] = 0.0

    # remove the block of sea ice in the Baltic sea caused by using the
    # 1986->2005 mean in months 04 to 11
    for m in range(4,12):
        sic_smooth[m::12,28:31,210:213] = 0.0
    
    if ey == 2101:
        # create 2101
        sic_smooth[-24:-12] = sic_smooth[-36:-24]
        sic_smooth[-12:] = sic_smooth[-36:-24]
    
    # restore the LSM
    sic_smooth[sst_input==mv] = 0.0
    # save the output
    save_sic(output, sic_smooth, lon_var, lat_var, time_var, mv)
    print output

#############################################################################

if __name__ == "__main__":
    in_file = ""
    out_file = ""
    cmip5_fit_fname = "/Users/Neil/Coding/CREDIBLE_output/output/rcp45_2006_2100/cmip5_polyfit_rcp45_2006_2100_1_anoms.nc"
    hadisst_deg = 1
    opts, args = getopt.getopt(sys.argv[1:], 'y:z:p:n:',
                               ['start_year=', 'end_year=', 'percentile=', 'number='])

    for opt, val in opts:
        if opt in ['--start_year', '-y']:
            sy = int(val)
        if opt in ['--end_year', '-z']:
            ey = int(val)
        if opt in ['--percentile', '-p']:
            p = float(val)
        if opt in ['--number', '-n']:
            n = int(val)
            
    rt = "rcp85"
    rfs = 1986
    rfe = 2005
    path = get_Ma_output_directory(rt, rfs, rfe, sy, ey)+"/"
    in_file = get_Ma_output_name(rt, rfs, rfe, sy, ey, p, n)
    out_file = in_file.replace("SST", "SIC")
    create_Ma_SIC(path+in_file, path+out_file, cmip5_fit_fname, sy, ey, 400)