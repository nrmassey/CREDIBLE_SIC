#! /usr/bin/env python  
#############################################################################
#
# Program : create_HAPPI_future_SIC.py
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
# Date    : 24/05/16
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

def get_HAPPI_output_directory(run_type, sy, ey):
    histo_sy, histo_ey, rcp_sy, rcp_ey = get_start_end_periods()
    out_dir = "HAPPI_SST_"+run_type+"_"+str(sy)+"_"+str(ey)
    path = "/Users/Neil/Coding/HAPPI_output/"
    if not os.path.exists(path+out_dir):
        os.mkdir(path+out_dir)
    return path + out_dir

#############################################################################

def get_HAPPI_output_name(run_type, ref_start, ref_end, sy, ey, model):
    out_name = "HAPPI15_"+run_type+"_"+str(sy)+"01_"+str(ey)+\
              "12_added_to_OSTIA_"+str(ref_start)+"01_"+str(ref_end)+"12_1DEG_"+model+".nc"
    return out_name
    
#############################################################################

def create_HAPPI_SIC(input, output, cmip5_fit_fname, sy, ey, rn):
    # load each mapping filename
    cmip5_mapping, mv = load_mapping(cmip5_fit_fname)

    # load the hadisst monthly reference values
    sst_ref_fname, sic_ref_fname = get_HadISST_monthly_ref_filenames(1899, 2010, rn)
    sic_ref_data = load_data(sic_ref_fname, "sic")
    sst_ref_fname = get_HadISST_monthly_reference_fname(1899, 2010, 1986, 2005, rn)
    sst_ref_data = load_data(sst_ref_fname, "sst")
    lon_var, lat_var, time_var = get_syn_SST_lon_lat_time_vars(input)
    
    # load the sst file in
    sst_input = load_data(input, "tos")
    # sst data is in different order in file
    lenX2 = sst_input.shape[2] / 2
    sst_right = numpy.array(sst_input[:,:,lenX2:])
    sst_left  = numpy.array(sst_input[:,:,:lenX2])
    # recombine
    sst_input[:,:,:lenX2] = sst_right
    sst_input[:,:,lenX2:] = sst_left
    
    sst_input[sst_input > 1000] = mv
    
    # truncate to sy->ey years worth of mapping data
    ref_yr = 2015
    rcp_offset = -10 # fudge to use RCP2.6 with RCP4.5 data
    sub_mapping = cmip5_mapping[(sy-rcp_offset-ref_yr)/10:(ey-rcp_offset-ref_yr)/10]
       
    # reconstruct
    print "Constructing sea-ice"
    # remove the sst reference to produce the anomalies
    n_rpts = sst_input.shape[0] / sst_ref_data.shape[0]
    sub_sst_anoms = sst_input - numpy.tile(sst_ref_data, [n_rpts,1,1])

    # create the SIC from the SST
    syn_sic_anoms = calc_sic_from_sst(sub_sst_anoms, sub_mapping, mv, sub_mapping.shape[1]-1)
    # add the reference back on
#    n_rpts = syn_sic_anoms.shape[0] / sic_ref_data.shape[0]
#    syn_sic = syn_sic_anoms + numpy.tile(sic_ref_data, [n_rpts,1,1])
    syn_sic = syn_sic_anoms
    # apply a sigmoid, write out with multiple widths (0,6,12)
    W=0
    if (W != 0):
        syn_sic = 1.0 / (1.0 + numpy.exp(-W*(syn_sic-0.5)))

    # restore the LSM
    syn_sic[sst_input<-1000] = mv
    
    print "Filling ice holes"
    sic_filled = fill_ice(syn_sic, mv)
    
    print "Removing isolated ice"
    sic_removed = remove_isolated_ice(sic_filled, mv)

    # smooth the ice with a 3x1 smoothing window
    weights = numpy.ones([1,3,1], 'f')
    print "Smoothing"
    sic_removed[sst_input<-1000] = mv
    sic_smooth = window_smooth_3D(sic_removed, weights, mv, smooth_zero=False)
    sic_smooth[(sic_smooth > 1.0) & (sic_smooth != mv)] = 1.0
    for m in [0,1,2,3,4,10,11]:
        sic_smooth[m::12][sic_smooth[m::12] != mv] = numpy.abs(sic_smooth[m::12][sic_smooth[m::12] != mv])

    # remove the block of sea ice in the Baltic sea caused by using the
    # 1986->2005 mean in months 04 to 11
    for m in range(4,12):
        sic_smooth[m::12,28:31,210:213] = 0.0
    
    if ey == 2101:
        # create 2101
        sic_smooth[-24:-12] = sic_smooth[-36:-24]
        sic_smooth[-12:] = sic_smooth[-36:-24]
    
    # clamp to 0.0 to 1.0
    sic_smooth[sic_smooth > 1.0] = 1.0
    sic_smooth[(sic_smooth < 0.0) & (sic_smooth != mv)] = 0.0
    # restore the LSM
    sic_smooth[sst_input<-1000] = mv
    # switch the values back
    sic_right = numpy.array(sic_smooth[:,:,lenX2:])
    sic_left  = numpy.array(sic_smooth[:,:,:lenX2])
    # recombine
    sic_smooth[:,:,:lenX2] = sic_right
    sic_smooth[:,:,lenX2:] = sic_left

    # save the output
    save_sic(output, sic_smooth, lon_var, lat_var, time_var, mv)
    print output

#############################################################################

if __name__ == "__main__":
    in_file = ""
    out_file = ""
    cmip5_fit_fname = "/Users/Neil/Coding/CREDIBLE_output/output/rcp45_2006_2100/cmip5_polyfit_rcp45_2006_2100_1_anoms.nc"
    hadisst_deg = 1
    opts, args = getopt.getopt(sys.argv[1:], 'y:z:m:',
                               ['start_year=', 'end_year=', 'model='])

    model = ""
    for opt, val in opts:
        if opt in ['--start_year', '-y']:
            sy = int(val)
        if opt in ['--end_year', '-z']:
            ey = int(val)
        if opt in ['--model', '-m']:
            model = val
            
    rt = "RCP26"
    rfs = 1996
    rfe = 2016
    path = get_HAPPI_output_directory(rt, sy, ey)+"/"
    in_file = get_HAPPI_output_name(rt, rfs, rfe, sy, ey, model)
    out_file = in_file.replace("HAPPI15_", "HAPPI15_SIC_")
    out_file = out_file.replace(".nc", "_anom.nc")
    create_HAPPI_SIC(path+in_file, path+out_file, cmip5_fit_fname, sy, ey, 400)