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
from create_HadISST_SST_SIC_mapping import find_polyfit, calc_polynomial, save_SST_SIC_mapping, calc_polyfits
from netcdf_file import *
import numpy

#############################################################################

def get_CMIP5_ens_mean_directory(run_type, ref_start, ref_end):
    parent_dir = "/Users/Neil/Coding/CREDIBLE_output/output/"
    cmip5_dir = parent_dir + run_type + "_" + str(ref_start) + "_" + str(ref_end)
    return cmip5_dir

#############################################################################

def get_CMIP5_ens_mean_filename(run_type, ref_start, ref_end, var, hemi, anoms=False):
    # input path
    parent_dir = get_CMIP5_ens_mean_directory(run_type, ref_start, ref_end)+"/"+var
    if var == "tos":
        model = "Omon"
    elif var == "sic":
        model = "OImon"
        
    cmip5_ens_mean_name = "atlas_"+var+"_"+model+"_"+hemi+"_"+run_type+"_ens_mean_"+str(ref_start)+"01-"+str(ref_end)+"12_1x1"
    if anoms:
        cmip5_ens_mean_name += "_anoms.nc"
    else:
        cmip5_ens_mean_name += ".nc"
    
    return parent_dir + "/" + cmip5_ens_mean_name

#############################################################################

def get_CMIP5_SST_SIC_mapping_fname(start, end, run_type, deg, anoms=True):
    parent_dir = get_CMIP5_ens_mean_directory(run_type, start, end)
    if not os.path.exists(parent_dir):
        os.mkdir(parent_dir)
    out_name = "cmip5_polyfit_"+str(run_type)+"_"+str(start)+"_"+str(end)+"_"+str(deg)
    if anoms:
        out_name += "_anoms.nc"
    else:
        out_name += ".nc"
    return parent_dir + "/" + out_name

#############################################################################

def create_CMIP5_SST_SIC_mapping(start, end, run_type, max_deg, anoms=False):
    # load the HadISST file - get the name from the run number
    ref_start = 2006
    ref_end = 2100
    cmip5_sic_arctic_name = get_CMIP5_ens_mean_filename(run_type, ref_start, ref_end, "sic", "arctic", anoms)
    cmip5_tos_arctic_name = get_CMIP5_ens_mean_filename(run_type, ref_start, ref_end, "tos", "arctic", anoms)
    
    cmip5_sic_antarctic_name = get_CMIP5_ens_mean_filename(run_type, ref_start, ref_end, "sic", "antarctic", anoms)
    cmip5_tos_antarctic_name = get_CMIP5_ens_mean_filename(run_type, ref_start, ref_end, "tos", "antarctic", anoms)
    
    print cmip5_sic_arctic_name
    
    nc_arctic_sic_fh = netcdf_file(cmip5_sic_arctic_name)
    nc_arctic_tos_fh = netcdf_file(cmip5_tos_arctic_name)
    nc_antarctic_sic_fh = netcdf_file(cmip5_sic_arctic_name)
    nc_antarctic_tos_fh = netcdf_file(cmip5_tos_arctic_name)
    
    # load the tos and sic
    tos_arctic_var = nc_arctic_tos_fh.variables["tos"]
    sic_arctic_var = nc_arctic_sic_fh.variables["sic"]
    tos_antarctic_var = nc_antarctic_tos_fh.variables["tos"]
    sic_antarctic_var = nc_antarctic_sic_fh.variables["sic"]
    
    # load the lat and lon
    lon_var = nc_arctic_tos_fh.variables["longitude"]
    lat_var = nc_arctic_tos_fh.variables["latitude"]
        
    # get the start and end points in the time axis
    base = 2006    # start of CMIP5 data
    nm=12
    s = (start-base)*nm
    e = (end-base)*nm

    # get the output filename
    out_name = get_CMIP5_SST_SIC_mapping_fname(start, end, run_type, max_deg, anoms)

    # do the calculation
    ant_s = lat_var.shape[0] / 2
    ant_e = lat_var.shape[0]

    # calculate the arctic and antarctic separately
    arctic_polyfits, arctic_degrees, mv = calc_polyfits(tos_arctic_var, sic_arctic_var, max_deg, s, e, 0, ant_s, anoms)
    antarctic_polyfits, antarctic_degrees, mv = calc_polyfits(tos_antarctic_var, sic_antarctic_var, max_deg, s, e, ant_s, ant_e, anoms)

    # put the antarctic fits into the arctic fits
    arctic_polyfits[:,:,ant_s:ant_e,:] = antarctic_polyfits[:,:,ant_s:ant_e,:]
    arctic_degrees[:,ant_s:ant_e,:] = antarctic_degrees[:,ant_s:ant_e,:]
    
    # save the polyfits
    save_SST_SIC_mapping(arctic_polyfits, arctic_degrees, out_name, lat_var, lon_var, mv)
    nc_arctic_sic_fh.close()
    nc_arctic_tos_fh.close()
    nc_antarctic_sic_fh.close()
    nc_antarctic_tos_fh.close()


#############################################################################

if __name__ == "__main__":
    start = 1978
    end   = 2010
    run_type = "rcp45"
    deg   = 4
    opts, args = getopt.getopt(sys.argv[1:], 's:e:r:d:',
                               ['start=', 'end=', 'run_type=', 'deg='])

    for opt, val in opts:
        if opt in ['--start', '-s']:
            start = int(val)
        if opt in ['--end', '-e']:
            end = int(val)
        if opt in ['--run_type', '-r']:
            run_type = val
        if opt in ['--deg', '-d']:
            deg = int(val)
    anoms = True
    create_CMIP5_SST_SIC_mapping(start, end, run_type, deg, anoms)