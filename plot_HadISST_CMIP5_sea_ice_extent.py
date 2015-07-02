#! /usr/bin/env python  
#############################################################################
#
# Program : create_HadISST_CMIP5_SIC_from_SST.py
# Author  : Neil Massey
# Purpose : 
# Inputs  : 
# Output  : in the output directory:
#           
#            
# Date    : 11/06/15
#
#############################################################################

import os, sys, getopt
sys.path.append("../CREDIBLE_SST")
import numpy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from netcdf_file import *
from create_HadISST_CMIP5_syn_SSTs import get_syn_sst_filename
from create_HadISST_sst_anoms import get_HadISST_input_filename, get_HadISST_smooth_fname, get_HadISST_month_smooth_filename
from create_CMIP5_sst_anoms import get_start_end_periods 
from calc_sea_ice_extent import *

#############################################################################

def plot_SST_SIC_extent(sst_fname, sic_fname, hadisst_fname):

    # load the data
    sst_fh = netcdf_file(sst_fname, 'r')
    sic_fh = netcdf_file(sic_fname, 'r')
    had_fh = netcdf_file(hadisst_fname, 'r')
    
    sst_var = sst_fh.variables["sst"]
    sic_var = sic_fh.variables["sic"]
    lat_var = sst_fh.variables["latitude"]
    lon_var = sst_fh.variables["longitude"]
    had_sic_var = had_fh.variables["sic"]
    had_sst_var = had_fh.variables["sst"]
    
    sst_data = sst_var[:]
    sic_data = sic_var[:]
    had_sic_data = had_sic_var[:]
    had_sst_data = had_sst_var[:]
    lat_data = lat_var[:]
    d_lon = lon_var[1] - lon_var[0]
    mv = sic_var._attributes["_FillValue"]
    
    # calculate the sea-ice extent
    sic_arctic, sic_antarctic = calc_sea_ice_extent(sic_data, lat_data, d_lon, mv, 1e-6)
    had_sic_arctic, had_sic_antarctic = calc_sea_ice_extent(had_sic_data, lat_data, d_lon, mv, 1e-6)
    
    sst_arctic, sst_antarctic = calc_sst_arctic_means(sst_data, lat_data, d_lon, mv)
    had_sst_arctic, had_sst_antarctic = calc_sst_arctic_means(had_sst_data, lat_data, d_lon, mv)
    
    gs = gridspec.GridSpec(2,10)
    
    sp0 = plt.subplot(gs[0,:])
    sp1 = sp0.twinx()
    sp2 = plt.subplot(gs[1,:])
    sp3 = sp2.twinx()
    x = [1899 + 1.0/12*i for i in range(0, sst_data.shape[0])]
    x2 =[1850 + 1.0/12*i for i in range(0, had_sst_data.shape[0])]
    print had_sst_data.shape[0]
#    for m in range(0,12):
    if True:
        m = 6
        sp1.plot(x[m::12], sic_arctic[m::12], 'r')
        sp3.plot(x[m::12], sic_antarctic[m::12], 'b')
        sp1.plot(x2[m::12], had_sic_arctic[m::12], 'k')
        sp3.plot(x2[m::12], had_sic_antarctic[m::12], '#808080')
        
        sp0.plot(x[m::12], sst_arctic[m::12], 'r', lw=2.0)
        sp2.plot(x[m::12], sst_antarctic[m::12], 'b', lw=2.0)
        sp0.plot(x2[m-1::12], had_sst_arctic[m::12], 'k', lw=2.0)
        sp2.plot(x2[m-1::12], had_sst_antarctic[m::12], '#808080', lw=2.0)
        
    plt.show()
    
    sst_fh.close()
    sic_fh.close()

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
    deg = 1
    opts, args = getopt.getopt(sys.argv[1:], 'r:s:e:n:f:a:i:v:d:m',
                               ['run_type=', 'ref_start=', 'ref_end=', 'neofs=', 'eof_year=', 'sample=', 'intvarmode=',
                                'varneofs=', 'monthly', 'degree'])

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
        if opt in ['--degree', '-d']:
            deg = int(val)
            
    # get the SST and SIC names
    sst_fname = get_syn_sst_filename(run_type, ref_start, ref_end, neofs, eof_year, sample, intvarmode, monthly)
    print sst_fname
    sic_fname = sst_fname.replace("ssts", "sic")[:-3] + "_" + str(deg) + ".nc"
    histo_sy, histo_ey, rcp_sy, rcp_ey = get_start_end_periods()
    hadisst_ey = 2010
#    hadisst_fname = get_HadISST_month_smooth_filename(histo_sy, hadisst_ey, 400)
    hadisst_fname = get_HadISST_input_filename(400)
    print hadisst_fname
    
    # plot the SST and SIC per month
    plot_SST_SIC_extent(sst_fname, sic_fname, hadisst_fname)