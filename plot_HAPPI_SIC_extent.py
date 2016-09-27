#! /usr/bin/env python  
#############################################################################
#
# Program : plot_syn_SIC_extent.py
# Author  : Neil Massey
# Purpose : Plot a timeseries of the ice extent of the synthetic SSTs
# Inputs  : run_type  : rcp4.5 | rc8.5 | histo
#           ref_start : year to start reference period, 1850->2005
#           ref_end   : year to end reference period, 1850->2005
#           year      : year to take warming at as an alternative to warm
# Notes   : all reference values are calculated from the historical run_type
#           CMIP5 ensemble members are only included if their historical run 
#           includes the reference period
#           requires Andrew Dawsons eofs python libraries:
#            http://ajdawson.github.io/eofs/
# Output  : 
# Date    : 15/09/15
#
#############################################################################

import sys,os,getopt
from calc_sea_ice_extent import *
sys.path.append("../CREDIBLE_SST")
from cmip5_functions import get_output_directory
from scipy.io.netcdf import *
import numpy
import matplotlib.pyplot as plt

#############################################################################

def get_sic_fname(sy, ey):
    out_dir = "/Users/Neil/Coding/HAPPI_output/HAPPI_SST_RCP26_"+str(sy)+"_"+str(ey)
    out_path = out_dir + "/HAPPI15_SIC_RCP26_"+str(sy)+"01_"+str(ey)+"12_added_to_OSTIA_200601_201512_1DEG_MMM.nc"
    return out_path

#############################################################################

def get_sic_extent_fname(run_type, ref_sy, ref_ey):
    in_fname = get_sic_fname(run_type, ref_start, ref_end, 0)
    out_fname = in_fname[:-21] + "sic_extent_yrmn.nc"
    return out_fname

#############################################################################

def calc_sic_extent_data(sy, ey):
    # 1. load the data
    # 2. calculate the sea ice extent for arctic and antarctic
    # 3. store it in a numpy array
    a = 0
    fname = get_sic_fname(sy,ey)
    fh = netcdf_file(fname)
    sic = fh.variables["sic"][:]
    lats = fh.variables["latitude"][:]
    lonv = fh.variables["longitude"]
    mv = fh.variables["sic"]._attributes["_FillValue"]
    lon_d = lonv[1] - lonv[0]
    grid_areas = calc_grid_areas(lats, lon_d)
    grid_areas = grid_areas.reshape([1, grid_areas.shape[0], 1])

    arctic_sic, antarctic_sic = calc_sea_ice_extent(sic, lats, lon_d, mv, 1e-3, grid_areas=grid_areas)

    sic_all_extent = numpy.zeros([2, arctic_sic.shape[0]], 'f')
    sic_all_extent[0] = arctic_sic
    sic_all_extent[1] = antarctic_sic
    a += 1
    fh.close()
    return sic_all_extent

#############################################################################

def plot_sic_extents(sp0, m, time_data, sic_extent_data):
    # load the 99 percentiles of the sic_extent data
    c = ['r', 'g', 'b', 'm', 'k']
    ls = []
    ldata = sic_extent_data[m::12]
    l = sp0.plot(time_data, ldata, c[0]+'-', lw=2, alpha=1.0, zorder=0)
    
#############################################################################

def plot_hadisst_sic_extent(sp0, h, m, yr, time_data):
    # load the hadisst_sic
    fname = "/Users/Neil/ClimateData/HadISST2/HadISST.2.1.0.0_realisation_dec2010_400.nc"
    ncfh = netcdf_file(fname)
    sic_var = ncfh.variables["sic"]
    lat_var = ncfh.variables["latitude"]
    mv = sic_var._attributes["_FillValue"]
    idx = m+(yr-1850)*12
    sic_data = sic_var[idx]
    print sic_data
    lat_data = lat_var[:]
    sic_arctic, sic_antarctic = calc_sea_ice_extent(sic_data, lat_data, 1.0, mv, S=1e-3)

    sic_extent = numpy.zeros([time_data.shape[0]], 'f')
    if h==0:
        sic_extent[:] = sic_arctic[:]
    else:
        sic_extent[:] = sic_antarctic[:]
    
    sp0.plot(time_data, sic_extent, 'k-', alpha=1.0, lw=3, zorder=2)
    ncfh.close()

#############################################################################


if __name__ == "__main__":
    year_start = -1
    year_end = -1
    opts, args = getopt.getopt(sys.argv[1:], 'y:z:',
                               ['year_start=', 'year_end='])

    for opt, val in opts:
        if opt in ['--year_start', '-y']:
            year_start = int(val)
        if opt in ['--year_end', '-z']:
            year_end = int(val)
            
    sic_extents = calc_sic_extent_data(year_start, year_end)
        
    sp0 = plt.subplot(211)
    sp1 = plt.subplot(212)
    time_data = numpy.arange(year_start,year_end+1)
    sp0.set_xlim([time_data[0],time_data[-1]])
    sp1.set_xlim([time_data[0],time_data[-1]])
 
    nh_mon = 8 # september
    sh_mon = 2 # february
    plot_sic_extents(sp0, nh_mon, time_data, sic_extents[0])
    plot_sic_extents(sp1, sh_mon, time_data, sic_extents[1])
    plot_hadisst_sic_extent(sp0, 0, nh_mon, 2007, time_data)
    plot_hadisst_sic_extent(sp1, 1, sh_mon, 2007, time_data)
     
    sp0.set_xlabel("Year")
    sp0.set_ylabel("Sea ice total area, $km^2$")
    sp0.set_title("Total sea ice area NH")
     
    sp1.set_xlabel("Year")
    sp1.set_ylabel("Sea ice total area, $km^2$")
    sp1.set_title("Total sea ice area SH")
    plt.gcf().set_size_inches(12, 8, forward=True)
     
    plt.tight_layout()
     
    outname = "HAPPI_sic_ext_ts_NH_SH.pdf"
    plt.savefig(outname)
    print outname
