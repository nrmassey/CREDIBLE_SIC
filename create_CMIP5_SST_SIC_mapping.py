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
from create_HadISST_SST_SIC_mapping import find_polyfit, calc_polynomial, save_SST_SIC_mapping, calc_polyfits, get_year_intervals
from create_CMIP5_sst_anoms import get_concat_anom_sst_ens_mean_fname
from netcdf_file import *
from cmip5_functions import load_sst_data
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

def create_CMIP5_SST_SIC_mapping(run_type, max_deg, anoms=True):
    # load the HadISST file - get the name from the run number
    ref_start = 2006
    ref_end = 2100
    sst_load_anoms = False
    cmip5_sic_arctic_name = get_CMIP5_ens_mean_filename(run_type, ref_start, ref_end, "sic", "arctic", anoms)
    cmip5_tos_arctic_name = get_CMIP5_ens_mean_filename(run_type, ref_start, ref_end, "tos", "arctic", sst_load_anoms)
    
    cmip5_sic_antarctic_name = get_CMIP5_ens_mean_filename(run_type, ref_start, ref_end, "sic", "antarctic", anoms)
    cmip5_tos_antarctic_name = get_CMIP5_ens_mean_filename(run_type, ref_start, ref_end, "tos", "antarctic", sst_load_anoms)
    
    nc_arctic_sic_fh = netcdf_file(cmip5_sic_arctic_name)
    nc_arctic_tos_fh = netcdf_file(cmip5_tos_arctic_name)
    nc_antarctic_sic_fh = netcdf_file(cmip5_sic_arctic_name)
    nc_antarctic_tos_fh = netcdf_file(cmip5_tos_arctic_name)
    
    # load the tos and sic
    tos_arctic_var = nc_arctic_tos_fh.variables["tos"]
    sic_arctic_var = nc_arctic_sic_fh.variables["sic"]
    tos_antarctic_var = nc_antarctic_tos_fh.variables["tos"]
    sic_antarctic_var = nc_antarctic_sic_fh.variables["sic"]
    mv = sic_antarctic_var._attributes["_FillValue"]
    
    # load the lat and lon
    lon_var = nc_arctic_tos_fh.variables["longitude"]
    lat_var = nc_arctic_tos_fh.variables["latitude"]

    # do the calculation
    ant_s = lat_var.shape[0] / 2
    ant_e = lat_var.shape[0]

    # byte swap if smoothing
    smooth = False
    if smooth:
        tos_arctic_data = tos_arctic_var[:].byteswap().newbyteorder()
        sic_arctic_data = sic_arctic_var[:].byteswap().newbyteorder()
        tos_antarctic_data = tos_antarctic_var[:].byteswap().newbyteorder()
        sic_antarctic_data = sic_antarctic_var[:].byteswap().newbyteorder()
    
        # smooth with 40 year filter
        P = 40
        sm_tos_arctic = running_gradient_3D(tos_arctic_data, P, mv)
        sm_sic_arctic = running_gradient_3D(sic_arctic_data, P, mv)
        sm_tos_antarctic = running_gradient_3D(tos_antarctic_data, P, mv)
        sm_sic_antarctic = running_gradient_3D(sic_antarctic_data, P, mv)
    else:
        sm_tos_arctic    = tos_arctic_var[:]
        sm_sic_arctic    = sic_arctic_var[:]
        sm_tos_antarctic = tos_antarctic_var[:]
        sm_sic_antarctic = sic_antarctic_var[:]
        
    # create the output - one array for degrees, one for the polyfit coefficients
    # size of degrees array is [years, months, lats, lons]
    # size of polyfits array is [degrees+1, years, months, lats, lons]
    years = get_year_intervals()
    cmip5_base_year = 2006
    
    n_years = 0
    for y in years:
        if y[0] >= cmip5_base_year:
            n_years += 1
    output_degrees  = numpy.zeros([n_years, 12, lat_var.shape[0], lon_var.shape[0]], 'f')
    output_polyfits = numpy.zeros([max_deg+1, n_years, 12, lat_var.shape[0], lon_var.shape[0]], 'f')
    
    # loop over the year intervals
    yi = 0
    for year in years:
        start = year[0]
        end = year[1]
        if start < cmip5_base_year:
            continue
        # calculate the start index and end index
        si = (start - cmip5_base_year)*12
        ei = (end - cmip5_base_year)*12
        
        # extrapolate beyond the 10 years to ensure a good continuity
        if si > 120:
            si -= 120
        if ei < sm_tos_arctic.shape[0] - 120:
            ei += 120

        # subset the smoothed data and add the ensemble mean trajectory
        sm_tos_arctic_local = sm_tos_arctic[si:ei]
        sm_sic_arctic_local = sm_sic_arctic[si:ei]
        sm_tos_antarctic_local = sm_tos_antarctic[si:ei]
        sm_sic_antarctic_local = sm_sic_antarctic[si:ei]

        # calculate the arctic and antarctic separately
        arctic_polyfits, arctic_degrees = calc_polyfits(sm_tos_arctic_local, sm_sic_arctic_local, max_deg, mv, 0, ant_s)
        antarctic_polyfits, antarctic_degrees, = calc_polyfits(sm_tos_antarctic_local, sm_sic_antarctic_local, max_deg, mv, ant_s, ant_e)

        # put the antarctic fits into the arctic fits
        arctic_polyfits[:,:,ant_s:ant_e,:] = antarctic_polyfits[:,:,ant_s:ant_e,:]
        arctic_degrees[:,ant_s:ant_e,:] = antarctic_degrees[:,ant_s:ant_e,:]
        
        # place in the outputs
        output_polyfits[:,yi,:,:,:] = arctic_polyfits
        output_degrees[yi] = arctic_degrees
        yi+=1

    # get the output filename
    out_name = get_CMIP5_SST_SIC_mapping_fname(cmip5_base_year, years[-1][1], run_type, max_deg, anoms)

    # save the polyfits
    save_SST_SIC_mapping(output_polyfits, output_degrees, out_name, years, cmip5_base_year, years[-1][1], lat_var, lon_var, mv)
        
    # close all the files
    nc_arctic_sic_fh.close()
    nc_arctic_tos_fh.close()
    nc_antarctic_sic_fh.close()
    nc_antarctic_tos_fh.close()

#############################################################################

if __name__ == "__main__":
    run_type = "rcp45"
    deg   = 4
    opts, args = getopt.getopt(sys.argv[1:], 'r:d:',
                               ['run_type=', 'deg='])

    for opt, val in opts:
        if opt in ['--run_type', '-r']:
            run_type = val
        if opt in ['--deg', '-d']:
            deg = int(val)
    anoms = True
    create_CMIP5_SST_SIC_mapping(run_type, deg, anoms)