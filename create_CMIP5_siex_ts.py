#! /usr/bin/env python  
#############################################################################
#
# Program : create_CMIP5_siex_ts.py
# Author  : Neil Massey
# Purpose : Create timeseries of yearly or monthly sea ice extent
# Inputs  : run_type  : rcp4.5 | rc8.5 | histo
#           ref_start : year to start reference period, 1850->2005
#           ref_end   : year to end reference period, 1850->2005
#           run_type  : historical | rcp45 | rcp85
# Notes   : all reference values are calculated from the historical run_type
#           CMIP5 ensemble members are only included if their historical run 
#           includes the reference period
# Output  : in the output/ directory filename is:
#            
# Date    : 05/05/15
#
#############################################################################

import os, sys, getopt

sys.path.append("../CREDIBLE_SST")
sys.path.append("../python_lib")
from create_CMIP5_sic_anoms import get_concat_anom_sic_output_fname
from create_CMIP5_sst_anoms import get_start_end_periods
from cmip5_functions import get_output_directory
from filter_cmip5_members import read_cmip5_index_file
import numpy, numpy.ma
from Geodetic import gA
from netcdf_file import *

#############################################################################

def get_siex_anom_ts_fname(run_type, ref_start, ref_end, monthly=False):
    out_dir = get_output_directory(run_type, ref_start, ref_end)
    out_name = out_dir + "/" + out_dir.split("/")[-1] + "_siex_anom_ts"
    if monthly:
        out_name += "_mon"
    out_name += ".nc"
    return out_name

#############################################################################

def create_siex_anoms_ts(run_type, ref_start, ref_end, monthly):
    # get the filtered list of CMIP5 ensemble members
    cmip5_rcp_idx = read_cmip5_index_file(run_type, ref_start, ref_end)
    histo_sy, histo_ey, rcp_sy, rcp_ey = get_start_end_periods()
    
    # create the storage
    n_ens = len(cmip5_rcp_idx)
    if monthly:
        f = 12
    else:
        f = 1
    n_t = (rcp_ey - histo_sy + 1)*f
    nh_siex = numpy.zeros([n_ens, n_t], 'f')
    sh_siex = numpy.zeros([n_ens, n_t], 'f')
    
    # create the area - load the lat and lon coords in
    fn = get_concat_anom_sic_output_fname(cmip5_rcp_idx[0][0], 
                                          cmip5_rcp_idx[0][1], run_type,
                                          ref_start, ref_end)
    fh = netcdf_file(fn, 'r')
    lon = fh.variables['longitude'][:]
    lat = fh.variables['latitude'][:]
    t_attr = fh.variables['time']._attributes
    t_vals = fh.variables['time'][:]
    fh.close()
    
    areas = numpy.zeros([lat.shape[0], lon.shape[0]], 'f')
    # loop through each latitude
    lat_d = lat[1] - lat[0]
    lon_d = lon[1] - lon[0]
    for y in range(0, lat.shape[0]):
        areas[y,:] = gA(lat[y]-lat_d, 0.0, lat[y]+lat_d, lon_d) / (1000**2)
    eq = lat.shape[0]*0.5
    
    for idx in range(0, n_ens):
        print cmip5_rcp_idx[idx][0] + ", " + cmip5_rcp_idx[idx][1] + ", " +str(idx)

        fn = get_concat_anom_sic_output_fname(cmip5_rcp_idx[idx][0], 
                                              cmip5_rcp_idx[idx][1], run_type,
                                              ref_start, ref_end)
        fh = netcdf_file(fn, 'r')
        sic_var = fh.variables['sic']
        mv = sic_var._attributes["_FillValue"]
        sic = numpy.ma.masked_equal(sic_var[:], mv)
        # sic expressed in percentages, hence 0.01
        ext = sic * 0.01 * areas
        nh_siex[idx] = numpy.sum(numpy.sum(ext[:,0:eq], axis=1), axis=1)
        sh_siex[idx] = numpy.sum(numpy.sum(ext[:,eq:], axis=1), axis=1)
        fh.close()
        
    # save the northern and southern hemisphere sea ice extent anomaly timeseries
    out_name = get_siex_anom_ts_fname(run_type, ref_start, ref_end, monthly)
    out_fh = netcdf_file(out_name, "w")
    # create dimensions and variables
    time_out_dim = out_fh.createDimension("time", t_vals.shape[0])
    time_out_var = out_fh.createVariable("time", t_vals.dtype, ("time",))
    ens_out_dim = out_fh.createDimension("ens", n_ens)
    ens_out_var = out_fh.createVariable("ens", 'f', ("ens",))
    nh_out_var = out_fh.createVariable("nh_siex", nh_siex.dtype, ("ens", "time",))
    sh_out_var = out_fh.createVariable("sh_siex", sh_siex.dtype, ("ens", "time",))
    # write out variables
    time_out_var._attributes = t_attr
    time_out_var[:] = t_vals[:]
    ens_out_var[:] = numpy.arange(0, n_ens)
    # write out data
    nh_out_var[:] = sh_siex[:]
    sh_out_var[:] = nh_siex[:]
    
    out_fh.close()


#############################################################################

if __name__ == "__main__":
    ref_start = -1
    ref_end = -1
    run_type = ""
    start_idx = 0
    end_idx = 0
    monthly = False
    opts, args = getopt.getopt(sys.argv[1:], 'r:s:e:m',
                               ['run_type=', 'ref_start=', 'ref_end=', 'monthly'])

    for opt, val in opts:
        if opt in ['--run_type', '-r']:
            run_type = val
        if opt in ['--ref_start', '-s']:
            ref_start = int(val)
        if opt in ['--ref_end', '-e']:
            ref_end = int(val)
        if opt in ['--monthly', '-m']:
            monthly = True
            
    create_siex_anoms_ts(run_type, ref_start, ref_end, monthly)