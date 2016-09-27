#! /usr/bin/env python  
#############################################################################
#
# Program : create_CMIP5_SST_SIC_mapping.py
# Author  : Neil Massey
# Purpose : Create a mapping between the sea-surface temperature (SST) and
#           sea-ice concentration (SIC) for the HadISST2 dataset
#           The mapping is created with the following algorithm:
#
#           For every grid box
#               For every month
#                   For sample pairs (SST,SIC) where SST < -1.79, SIC = 1.0
#                   Remove sample pairs (SST,SIC) where SIC == 0.0
#                   If there are more than 5 sample pairs
#                       Fit a polynomial of degree d to the sample pairs, where
#                           x = SST and y = SIC
#                           y = p[0]*pow(x,d) + p[1]*pow(x,d-1) + p[d]
#
# Inputs  : start : start year to fit the polynomial to
#           end   : end year to fit the polynomial to
#           deg   : degree of the polynomial to fit (recommend 4)
#           rn    : run number
# Output  : in the output directory:
#           hadisst_polyfit_<start>_<end>_<deg>_<rn>.nc
#           a netcdf file containing the polynomial fit for every month and every grid point
#           array format is: [d,m,lat,lon], where d is position in the polynomial above
#               m is month, lat and lon are the grid box
#            
# Date    : 10/06/15
#
#############################################################################

import os, sys, getopt
sys.path.append("../CREDIBLE_SST")
from cmip5_functions import load_data, get_missing_value
from sst_sic_mapping import *

import numpy

import pyximport
pyximport.install(setup_args={'include_dirs':[numpy.get_include()]})
sys.path.append("/Users/Neil/python_lib")
from running_gradient_filter import *

#############################################################################

def get_CMIP5_ens_mean_directory(run_type, ref_start, ref_end):
    parent_dir = "/Users/Neil/Coding/CREDIBLE_output/output/"
    cmip5_dir = parent_dir + run_type + "_" + str(ref_start) + "_" + str(ref_end)
    return cmip5_dir

#############################################################################

def get_CMIP5_ens_mean_anom_filename(run_type, ref_start, ref_end, var, hemi):
    # input path
    parent_dir = get_CMIP5_ens_mean_directory(run_type, ref_start, ref_end)+"/"
    if var == "tos":
        model = "Omon"
        parent_dir += "concat_sst_anoms"
    elif var == "sic":
        model = "OImon"
        parent_dir += "concat_sic_anoms"
    sst_start = 2006
    sst_end = 2100
    cmip5_ens_mean_name = "atlas_"+var+"_"+model+"_"+hemi+"_"+run_type+"_ens_mean_"+str(sst_start)+"01-"+str(sst_end)+"12_1x1"
    cmip5_ens_mean_name += "_anoms.nc"
    
    return parent_dir + "/" + cmip5_ens_mean_name

#############################################################################

def load_CMIP5_anom_data(run_type):
    # load the HadISST file - get the name from the run number
    ref_start = 1986
    ref_end = 2005
    
    # get the filenames
    cmip5_sic_arctic_anom_name = get_CMIP5_ens_mean_anom_filename(run_type, ref_start, ref_end, "sic", "arctic")
    cmip5_tos_arctic_anom_name = get_CMIP5_ens_mean_anom_filename(run_type, ref_start, ref_end, "tos", "arctic")
    
    cmip5_sic_antarctic_anom_name = get_CMIP5_ens_mean_anom_filename(run_type, ref_start, ref_end, "sic", "antarctic")
    cmip5_tos_antarctic_anom_name = get_CMIP5_ens_mean_anom_filename(run_type, ref_start, ref_end, "tos", "antarctic")
    
    # load the data
    cmip5_sic_arctic_anoms = load_data(cmip5_sic_arctic_anom_name, "sic")
    cmip5_tos_arctic_anoms = load_data(cmip5_tos_arctic_anom_name, "tos")
    time = load_data(cmip5_tos_arctic_anom_name, "time")

    cmip5_sic_antarctic_anoms = load_data(cmip5_sic_antarctic_anom_name, "sic")
    cmip5_tos_antarctic_anoms = load_data(cmip5_tos_antarctic_anom_name, "tos")
    
    # get the missing value
    mv = get_missing_value(cmip5_tos_arctic_anom_name, "tos")
    
    # amalgamate the data into one array, splitting at the equator and copying the
    # arctic data into the NH and the antarctic data into the SH
    ant_s = cmip5_sic_arctic_anoms.shape[1] / 2

    cmip5_sic_arctic_anoms[:,ant_s:,:] = cmip5_sic_antarctic_anoms[:,ant_s:,:]
    cmip5_tos_arctic_anoms[:,ant_s:,:] = cmip5_tos_antarctic_anoms[:,ant_s:,:]
    
    # trim the first 4 years - we only want 2010 to 2100
#    S = 12 * (2010-2006)
    S=0
    cmip5_tos_anom_data = cmip5_tos_arctic_anoms[S:]
    cmip5_sic_anom_data = cmip5_sic_arctic_anoms[S:]
    
    return cmip5_tos_anom_data, cmip5_sic_anom_data, mv

#############################################################################

def get_CMIP5_SST_SIC_mapping_fname(start, end, run_type, deg):
    parent_dir = get_CMIP5_ens_mean_directory(run_type, start, end)
    if not os.path.exists(parent_dir):
        os.mkdir(parent_dir)
    out_name = "cmip5_polyfit_"+str(run_type)+"_"+str(start)+"_"+str(end)+"_"+str(deg)
    out_name += "_anoms.nc"
    return parent_dir + "/" + out_name

#############################################################################

def get_CMIP5_lon_lat_vars(run_type):
    ref_start = 1986
    ref_end = 2005
    cmip5_sic_arctic_anom_name = get_CMIP5_ens_mean_anom_filename(run_type, ref_start, ref_end, "sic", "arctic")
    fh = netcdf_file(cmip5_sic_arctic_anom_name)
    
    lon_var = fh.variables["longitude"]
    lat_var = fh.variables["latitude"]
    return lon_var, lat_var

#############################################################################

def create_CMIP5_SST_SIC_mapping(run_type, deg):
    # load the data
    cmip5_tos_anoms, cmip5_sic_anoms, mv = load_CMIP5_anom_data(run_type)
            
    # year intervals
    years = [[y, y+10] for y in range(2010,2100,10)]
    
    # calculate the mapping
    mapping = calc_sic_mapping(cmip5_tos_anoms, cmip5_sic_anoms, mv, deg)
    
    # get the output filename
    out_name = get_CMIP5_SST_SIC_mapping_fname(2006, years[-1][1], run_type, deg)
    
    # save the polyfits
    lon_var, lat_var = get_CMIP5_lon_lat_vars(run_type)
    save_mapping(out_name, mapping, lat_var, lon_var, years, mv, deg)
    
#############################################################################

if __name__ == "__main__":
    run_type = "rcp45"
    deg = 2
    opts, args = getopt.getopt(sys.argv[1:], 'r:d:',
                               ['run_type=', 'deg='])

    for opt, val in opts:
        if opt in ['--run_type', '-r']:
            run_type = val
        if opt in ['--deg', '-d']:
            deg = int(val)
    create_CMIP5_SST_SIC_mapping(run_type, deg)