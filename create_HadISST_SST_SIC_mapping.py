#! /usr/bin/env python  
#########################################################################################
#
# Program : create_HadISST_SST_SIC_mapping.py
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
#########################################################################################

import os, sys, getopt
sys.path.append("../CREDIBLE_SST")
from cmip5_functions import load_data, get_missing_value
from create_HadISST_sst_anoms import get_HadISST_input_filename, get_HadISST_output_directory 
from sst_sic_mapping import *
from cdo import *

import numpy
from netcdf_file import *

#########################################################################################

def get_HadISST_SST_SIC_mapping_fname(start, end, rn, deg):
    out_dir = get_HadISST_output_directory(start, end, rn)
    out_name = "hadisst_polyfit_"+str(start)+"_"+str(end)+"_"+str(deg)+"_"+str(rn)
    out_name += "_anoms.nc"
    return out_dir + "/" + out_name

#########################################################################################

def get_HadISST_monthly_anomaly_filenames(start, end, rn):
    out_dir = get_HadISST_output_directory(start, end, rn)
    sic_fname = out_dir + "/HadISST_" + str(rn) + "_mon_sic_anoms.nc"
    sst_fname = out_dir + "/HadISST_" + str(rn) + "_mon_sst_anoms.nc"
    return sst_fname, sic_fname

#########################################################################################

def get_HadISST_monthly_ref_filenames(start, end, rn):
    out_dir = get_HadISST_output_directory(start, end, rn)
    sic_fname = out_dir + "/HadISST_" + str(rn) + "_mon_sic_ref.nc"
    sst_fname = out_dir + "/HadISST_" + str(rn) + "_mon_sst_ref.nc"
    return sst_fname, sic_fname

#########################################################################################

def load_HadISST_anom_data(rn):
    start = 1899
    end = 2010

    sst_anom_fname, sic_anom_fname = get_HadISST_monthly_anomaly_filenames(start, end, rn)
    sst_anom_data = load_data(sst_anom_fname, "sst")
    sic_anom_data = load_data(sic_anom_fname, "sic")
    mv = get_missing_value(sst_anom_fname, "sst")
    
    return sst_anom_data, sic_anom_data, mv

#########################################################################################

def load_HadISST_ref_data(rn):
    start = 1899
    end = 2010

    sst_ref_fname, sic_ref_fname = get_HadISST_monthly_ref_filenames(start, end, rn)
    sic_ref_data = load_data(sic_ref_fname, "sic")
    sst_ref_data = load_data(sst_ref_fname, "sst")
    
    return sst_ref_data, sic_ref_data

########################################################################################

def get_HadISST_lon_lat_vars(rn):
    start = 1899
    end = 2010
    sst_fname, sic_fname = get_HadISST_monthly_anomaly_filenames(start, end, rn)

    fh = netcdf_file(sst_fname)
    lon_var = fh.variables["longitude"]
    lat_var = fh.variables["latitude"]
    return lon_var, lat_var

########################################################################################

def create_HadISST_SST_SIC_monthly_anoms(rn):
    start = 1899
    end = 2010

    hadisst_fname = get_HadISST_input_filename(rn)
    sst_anom_fname, sic_anom_fname = get_HadISST_monthly_anomaly_filenames(start, end, rn)
    # do the sst anoms
    cdo = Cdo()
    sst_cdo_string=" -ymonmean -selvar,sst -selyear,"+str(1986)+"/"+str(2005)+" " +hadisst_fname
    cdo.ymonsub(input=" -selvar,sst "+ hadisst_fname + sst_cdo_string, output=sst_anom_fname)

    sic_cdo_string=" -ymonmean -selvar,sic -selyear,"+str(1986)+"/"+str(2005)+" " +hadisst_fname
    cdo.ymonsub(input=" -selvar,sic "+ hadisst_fname + sic_cdo_string, output=sic_anom_fname)

########################################################################################

def create_HadISST_SST_SIC_mapping(rn, deg):
    # load the anomalies from the yearly month means of 1986->2010
    sst_anom_data, sic_anom_data, mv = load_HadISST_anom_data(rn)
        
    # year intervals
    years = [[y, y+10] for y in range(1850,2010,10)]
    # get the output filename
    out_name = get_HadISST_SST_SIC_mapping_fname(years[0][0], years[-1][1], rn, deg)
    print out_name

    # calculate the mapping between the anomalies
    mapping = calc_sic_mapping(sst_anom_data, sic_anom_data, mv, deg)
    
    # save the polyfits
    lon_var, lat_var = get_HadISST_lon_lat_vars(rn)
    save_mapping(out_name, mapping, lat_var, lon_var, years, mv, deg)

#############################################################################

if __name__ == "__main__":
    run_n = 400
    deg   = 2
    opts, args = getopt.getopt(sys.argv[1:], 'n:d:',
                               ['runn=', 'deg='])

    for opt, val in opts:
        if opt in ['--runn', '-n']:
            run_n = val
        if opt in ['--deg', '-d']:
            deg = int(val)
    
    create_HadISST_SST_SIC_monthly_anoms(run_n)
    create_HadISST_SST_SIC_mapping(run_n, deg)