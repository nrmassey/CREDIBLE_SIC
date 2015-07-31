#! /usr/bin/env python  
#############################################################################
#
# Program : create_HadISST_CMIP5_SIC_from_SST.py
# Author  : Neil Massey
# Purpose : Create a sea-ice concentration file (SIC) from the synthetic SST
#           file created by ./create_HadISST_CMIP5_syn_SSTs.py
#           To create the mapping run ./create_HadISST_SST_SIC_mapping.py
#                                  or ./create_CMIP5_SST_SIC_mapping.py
# Inputs  : same as ./create_HadISST_CMIP5_syn_SSTs.py
# Output  : in the output directory:
#           output_name
#            
# Date    : 11/06/15
#
#############################################################################

import os, sys, getopt
sys.path.append("../CREDIBLE_SST")
from create_HadISST_CMIP5_syn_SSTs import get_syn_sst_filename, get_output_directory, create_hadisst_annual_cycle
from create_CMIP5_SST_SIC_mapping import get_CMIP5_SST_SIC_mapping_fname
from create_HadISST_SST_SIC_mapping import get_HadISST_SST_SIC_mapping_fname, get_year_intervals
from create_SIC_from_mapping import create_SIC_from_mapping, write_SIC
from create_HadISST_sst_anoms import get_HadISST_reference_fname
from create_CMIP5_sst_anoms import get_start_end_periods, save_3d_file, get_concat_anom_sst_ens_mean_smooth_fname
from cmip5_functions import load_sst_data
from netcdf_file import netcdf_file
import numpy
import pyximport
pyximport.install(setup_args={'include_dirs':[numpy.get_include()]})
from fill_ice import *
from zonal_smoother import *
from remove_isolated_ice import *

#############################################################################

def load_sst_sic_mapping(run_type, use_anoms=True, max_deg=4):
    deg = max_deg
    rn = "400"

    # get the hadisst and cmip5 mapping filenames    
    cmip5_map_fname = get_CMIP5_SST_SIC_mapping_fname(2006, 2100, run_type, deg, use_anoms)
    hadisst_map_fname = get_HadISST_SST_SIC_mapping_fname(1899, 2010, rn, deg, use_anoms)
        
    # load the CMIP5 mapping
    cmip5_map_fh = netcdf_file(cmip5_map_fname, 'r')
    cmip5_polyfit_data = cmip5_map_fh.variables["polyfit"][:]
    cmip5_degree_data  = cmip5_map_fh.variables["degrees"][:]
    cmip5_map_fh.close()
    
    # load the HadISST mapping
    hadisst_map_fh = netcdf_file(hadisst_map_fname, 'r')
    hadisst_polyfit_data = hadisst_map_fh.variables["polyfit"][:]
    hadisst_degree_data  = hadisst_map_fh.variables["degrees"][:]
    hadisst_map_fh.close()
    
    return [hadisst_polyfit_data, hadisst_degree_data, cmip5_polyfit_data, cmip5_degree_data]

#############################################################################

def impose_lsm(data, mv):
    lsm_path = "/Users/Neil/ClimateData/HadISST2/HadISST.2.1.0.0_sst_lsm.nc"
    fh = netcdf_file(lsm_path)
    lsm_var = fh.variables["sst"]
    lsm_mv = lsm_var._attributes["_FillValue"]
    lsm_data = lsm_var[:]
    data[:,lsm_data[0]==lsm_mv] = mv
    return data

#############################################################################

def load_syn_SSTs(run_type, ref_start, ref_end, neofs, eof_year, sample, intvarmode, monthly):
    sst_fname = get_syn_sst_filename(run_type, ref_start, ref_end, neofs, eof_year, sample, intvarmode, monthly)
    sst_fh = netcdf_file(sst_fname)
     
    sst_var = sst_fh.variables["sst"]
    lat_var = sst_fh.variables["latitude"]
    mv = sst_var._attributes["_FillValue"]
    sst_data =  numpy.array(sst_var[:]).byteswap().newbyteorder()
    lats = numpy.array(lat_var[:]).byteswap().newbyteorder()

    return sst_data, lats, mv

#############################################################################

def calc_sea_ice(run_type, ref_start, ref_end, neofs, eof_year, sample, intvarmode, monthly, use_anoms=True, max_deg=4):

    # calculate the sea ice using the following method.
    # (all SIC is now expressed as anomalies wrt to anomalies of SST)
    # 1. Load in the timeseries of SST data
    # 2. Subtract the 1986->2005 reference period
    # 3. Subtract the annual cycle residuals (after tiling to n_years / 12)
    # 4. Use this (the monthly SST anomalies from the reference SST) to predict the SIC anomalies
    # 5. Add on the (monthly-varying) SIC reference, after tiling to n_years/12

    # get the time periods
    histo_sy, histo_ey, rcp_sy, rcp_ey = get_start_end_periods()
    hadisst_ey = 2010
    run_n = 400

    # read the sst var in
    sst_data, lats, mv = load_syn_SSTs(run_type, ref_start, ref_end, neofs, eof_year, sample, intvarmode, monthly)    
    
    # read the hadisst sst reference (1 timesteps) in
    hadisst_sst_ref_fname = get_HadISST_reference_fname(histo_sy, hadisst_ey, ref_start, ref_end, run_n)
    hadisst_sst_ref = load_sst_data(hadisst_sst_ref_fname, "sst")
    
    # read the hadisst sst annual cycle anomalies from the reference (12 timesteps) and tile it over the sst_data timeseries
    n_repeats = sst_data.shape[0] / 12
    hadisst_sst_ac_tile = create_hadisst_annual_cycle(run_type, ref_start, ref_end, n_repeats, run_n=400)

    cmip5_ens_mean_anoms_fname = get_concat_anom_sst_ens_mean_smooth_fname(run_type, ref_start, ref_end, monthly=monthly)
    cmip5_ens_mean_anoms = load_sst_data(cmip5_ens_mean_anoms_fname, "sst")

    # calculate the anomalies as sst_data - hadisst_ref_sst - hadisst_ref_sst_ac
    sst_anoms = sst_data - hadisst_sst_ac_tile - hadisst_sst_ref
    
    # read in and tile the SIC annual cycle
    hadisst_sic_ref_fname = hadisst_sst_ref_fname[:-3] + "_sic.nc"
    hadisst_sic_ac_ref = load_sst_data(hadisst_sic_ref_fname, "sic")
    hadisst_sic_ac_ts = numpy.tile(hadisst_sic_ac_ref, [n_repeats,1,1])
    
    # get the year intervals
    years = get_year_intervals()
    #
    sy = years[0][0]
    # create the output values
    sic_out = numpy.zeros(sst_data.shape, 'f')
    yi = 0
    
    # create the output filename
    sst_fname = get_syn_sst_filename(run_type, ref_start, ref_end, neofs, eof_year, sample, intvarmode, monthly)
    out_fname = sst_fname.replace("ssts", "sic")
    out_fname = out_fname[:-3] + "_" + str(max_deg) + ".nc"

    # load the sst / sic mapping   
    sst_sic_map = load_sst_sic_mapping(run_type, use_anoms, max_deg)
    
    fh = netcdf_file(sst_fname, 'r')
    lat_var = fh.variables["latitude"]
    lon_var = fh.variables["longitude"]
    time_var = fh.variables["time"]
    attrs = fh.variables["sst"]._attributes

    # save out the SST anoms as a test to compare 
#    out_sst_anoms_fname = sst_fname[:-3] + "_anoms.nc"
#    save_3d_file(out_sst_anoms_fname, sst_anoms, lon_var, lat_var, attrs, time_var, "sst_anoms")
#    print out_sst_anoms_fname

    # create the SIC from the SST->SIC mapping
    sic_out = create_SIC_from_mapping(sst_anoms, time_var, sst_sic_map, mv, use_anoms)

    # Add the annual cycle back on
    sic_out += hadisst_sic_ac_ts
    # enforce the limits
    sic_idx = numpy.where((sic_out<0.0) & (sic_out > -1000.0))
    sic_out[sic_idx] = 0.0
    sic_idx = numpy.where((sic_out>1.0) & (sic_out < 1e6))
    sic_out[sic_idx] = 1.0

    # set very small amounts of SIC to zero
    sic_out[(sic_out < 0.1) & (sic_out != mv)] = 0.0
    # impose the HadISST sst lsm
    sic_out = impose_lsm(sic_out, mv)
    # fill in any holes in the ice
    sic_fill = fill_ice(sic_out, mv)
    # remove isolated pixels in the ice
    sic_remove = remove_isolated_ice(sic_fill, mv)
    # apply a zonal smoother
    sic_zonsm = zonal_smoother(sic_remove, lats, 64, 15, mv)
    sic_out = sic_zonsm
    # write out
    write_SIC(out_fname, sic_out, lon_var, lat_var, time_var, mv)
    fh.close()
    print out_fname

#############################################################################

if __name__ == "__main__":
    ref_start = -1
    ref_end = -1
    run_type = ""       # run type can be rcp45, rcp85, rcp26 or likely -
                        # fits ensemble mean to AR5 ch11 likely scenario for 2016->2035
    neofs = 0
    eof_year = 2050
    sample = 100
    intvarmode = 0      # internal variability mode - 0 = none, 1 = yearly, 2 = monthly
    monthly = False     # use the monthly EOFs / PCs ?
    max_deg = 4
    opts, args = getopt.getopt(sys.argv[1:], 'r:s:e:n:f:a:i:v:d:m',
                               ['run_type=', 'ref_start=', 'ref_end=', 'neofs=', 'eof_year=', 'sample=', 'intvarmode=',
                                'varneofs=', 'monthly', 'deg='])

    for opt, val in opts:
        if opt in ['--run_type', '-r']:
            run_type = val
        if opt in ['--ref_start', '-s']:
            ref_start = int(val)
        if opt in ['--ref_end', '-e']:
            ref_end = int(val)
        if opt in ['--neofs', '-n']:
            neofs = int(val)
        if opt in ['--eof_year', '-f']:
            eof_year = int(val)
        if opt in ['--sample', '-a']:
            sample = int(val)
        if opt in ['--intvar', '-i']:
            intvarmode = int(val)
        if opt in ['--monthly', '-m']:
            monthly = True
        if opt in ['--deg', '-d']:
            max_deg = int(val)

    if max_deg > 1:
        print "Fitting using degs > 1 not currently implemented"
        sys.exit()

    use_anoms=True
    calc_sea_ice(run_type, ref_start, ref_end, neofs, eof_year, sample, intvarmode, monthly, use_anoms, max_deg)
