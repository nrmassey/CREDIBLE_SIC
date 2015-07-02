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
from create_HadISST_CMIP5_syn_SSTs import get_syn_sst_filename, get_output_directory
from create_CMIP5_SST_SIC_mapping import get_CMIP5_SST_SIC_mapping_fname
from create_HadISST_SST_SIC_mapping import get_HadISST_SST_SIC_mapping_fname
from create_SIC_from_mapping import create_SIC_from_mapping, write_SIC
from create_HadISST_sst_anoms import get_HadISST_reference_fname
from create_CMIP5_sst_anoms import get_start_end_periods
from netcdf_file import netcdf_file
import numpy
import pyximport
pyximport.install(setup_args={'include_dirs':[numpy.get_include()]})
from fill_ice import *
from zonal_smoother import *
from remove_isolated_ice import *

#############################################################################

def get_year_intervals():
#    years = [[1899,1910], [1909,1920], [1919,1930], [1929,1940], [1939,1950], 
#             [1949,1960], [1959,1970], [1969,1980], [1979,1990], [1989,2000],
#             [1999,2010]]
    years = [[1899,1940], [1929,1970], [1959,2000], [1989,2010],
             [2006,2040], [2039,2070], [2059,2100]]

    return years
    
#############################################################################

def load_sst_sic_mapping(run_type, use_anoms=False, max_deg=4):
    deg = max_deg
    rn = "400"
    
    years = get_year_intervals()
    maps = []
    
    for year in years:
        anoms = False
        if year[0] > 2000:
            # load the mappings for CMIP5
            anoms = use_anoms
            map_fname = get_CMIP5_SST_SIC_mapping_fname(year[0], year[1], run_type, deg, anoms)
        else:
            # load the mappings for HadISST
            map_fname = get_HadISST_SST_SIC_mapping_fname(year[0], year[1], rn, deg)
        
        # load the mapping
        map_fh = netcdf_file(map_fname, 'r')
        polyfit_data = map_fh.variables["polyfit"][:]
        degree_data = map_fh.variables["degrees"][:]
        maps.append([polyfit_data, degree_data, map_fname, anoms])
        map_fh.close()
    return maps

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

def calc_sea_ice(sst_fname, sst_sic_map, ref_fname, sic_ref_fname, use_anoms=False, max_deg=4):
    # read the sst var in
    sst_fh = netcdf_file(sst_fname)
    
    sst_var = sst_fh.variables["sst"]
    lon_var = sst_fh.variables["longitude"]
    lat_var = sst_fh.variables["latitude"]
    time_var = sst_fh.variables["time"]
    mv = sst_var._attributes["_FillValue"]
    
    lats = numpy.array(lat_var[:]).byteswap().newbyteorder()

    # read the reference file in
    ref_fh = netcdf_file(ref_fname)
    ref_var = ref_fh.variables["sst"]
    ref_data = numpy.array(ref_var[:])
    ref_fh.close()

    # read the sic reference file in
    sic_ref_fh = netcdf_file(sic_ref_fname)
    sic_ref_var = sic_ref_fh.variables["sic"]
    sic_ref_data = numpy.array(sic_ref_var[:])
    sic_ref_fh.close()

    # get the year intervals
    years = get_year_intervals()
    #
    sy = years[0][0]
    # create the output values
    sic_out = numpy.zeros(sst_var.shape, 'f')
    yi = 0
    
    # create the output filename
    out_fname = sst_fname.replace("ssts", "sic")
    out_fname = out_fname[:-3] + "_" + str(max_deg) + ".nc"
    print out_fname
    # rejig the first CMIP5 period to start in 1999
    years[4][0] = 1999
    
    for yeari in years:
        # get the indices into the sst data
        si = (yeari[0]-sy)*12
        ei = (yeari[1]-sy+1)*12
        print yeari[0], yeari[1], sst_sic_map[yi][2], sst_sic_map[yi][3]
        sst_data = sst_var[si:ei]
        # should we express this as an anomaly?
        if sst_sic_map[yi][3] == True:
            sst_data = sst_data - ref_data
        sic_data = create_SIC_from_mapping(sst_data, time_var[si:ei], sst_sic_map[yi][0], sst_sic_map[yi][1], mv, use_anoms)
        # add on the sic reference if this is an anomaly
        if sst_sic_map[yi][3] == True:
            n_rpts = (yeari[1]-yeari[0]+1)
            add_sic_ref = numpy.tile(sic_ref_data, (n_rpts,1,1))
            sic_data += add_sic_ref
        
        # enforce the limits
        sic_idx = numpy.where((sic_data<0.0) & (sic_data > -1000.0))
        sic_data[sic_idx] = 0.0
        sic_idx = numpy.where((sic_data>1.0) & (sic_data < 1000.0))
        sic_data[sic_idx] = 1.0
        
        sic_out[si:ei] = sic_data
        yi += 1

    sic_out[(sic_out < 0.1) & (sic_out != mv)] = 0.0
    sic_out = impose_lsm(sic_out, mv)
    sic_fill = fill_ice(sic_out, mv)
    sic_remove = remove_isolated_ice(sic_fill, mv)
    sic_zonsm = zonal_smoother(sic_remove, lats, 64, 15, mv)
    sic_out = sic_zonsm
    # write out
    write_SIC(out_fname, sic_out, lon_var, lat_var, time_var, mv)

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
    
    histo_sy, histo_ey, rcp_sy, rcp_ey = get_start_end_periods()
    hadisst_ey = 2010
    sst_fname = get_syn_sst_filename(run_type, ref_start, ref_end, neofs, eof_year, sample, intvarmode, monthly)
    ref_fname = get_HadISST_reference_fname(histo_sy, hadisst_ey, ref_start, ref_end, 400)
    sic_ref_fname = ref_fname[:-3] + "_sic.nc"
    sst_sic_map = load_sst_sic_mapping(run_type, use_anoms=True, max_deg=max_deg)
    calc_sea_ice(sst_fname, sst_sic_map, ref_fname, sic_ref_fname, use_anoms=True, max_deg=max_deg)
