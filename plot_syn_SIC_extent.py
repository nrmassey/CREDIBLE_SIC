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
from netcdf_file import *
import numpy
import matplotlib.pyplot as plt

#############################################################################

def get_sic_fname(run_type, ref_sy, ref_ey, sample):
    out_dir_base = "/Users/Neil/Coding/CREDIBLE_output/output/"
    out_dir = out_dir_base + "HadISST_1899_2005_"+run_type+"_2006_2100_r"+str(ref_sy)+"_"+str(ref_ey)+"_y2050/yrmns/sic/"
    out_path = out_dir + "/HadISST_1899_2005_"+run_type+"_2006_2100_r"+str(ref_sy)+"_"+str(ref_ey)+"_y2050_f2050_n6_a"+str(sample)+"_varmon_sic_yrmn.nc"
    return out_path

#############################################################################

def get_sic_extent_fname(run_type, ref_sy, ref_ey):
    in_fname = get_sic_fname(run_type, ref_start, ref_end, 0)
    out_fname = in_fname[:-21] + "sic_extent_yrmn.nc"
    return out_fname

#############################################################################

def calc_sic_extent_data(run_type, ref_sy, ref_ey):
    # 1. load the data
    # 2. calculate the sea ice extent for arctic and antarctic
    # 3. store it in a numpy array
    for a in range(1,100):
        fname = get_sic_fname(run_type, ref_start, ref_end, a)
        fh = netcdf_file(fname)
        sic = fh.variables["sic"][:]
        lats = fh.variables["latitude"][:]
        lonv = fh.variables["longitude"]
        mv = fh.variables["sic"]._attributes["_FillValue"]
        lon_d = lonv[1] - lonv[0]
        if a == 1:
            grid_areas = calc_grid_areas(lats, lon_d)
            grid_areas = grid_areas.reshape([1, grid_areas.shape[0], 1])

        arctic_sic, antarctic_sic = calc_sea_ice_extent(sic, lats, lon_d, mv, 1e-3, grid_areas=grid_areas)

        if a == 1:
            sic_all_extent = numpy.zeros([2, 99, arctic_sic.shape[0]], 'f')
        sic_all_extent[0, a-1] = arctic_sic
        sic_all_extent[1, a-1] = antarctic_sic
        fh.close()
        
    # save the file out
    
    out_fname = get_sic_extent_fname(run_type, ref_sy, ref_ey)
    out_fh = netcdf_file(out_fname, "w")
    regi_out_dim = out_fh.createDimension("region", sic_all_extent.shape[0])
    samp_out_dim = out_fh.createDimension("sample", sic_all_extent.shape[1])
    time_out_dim = out_fh.createDimension("time", sic_all_extent.shape[2])
    
    regi_out_var = out_fh.createVariable("region", sic_all_extent.dtype, ("region",))
    samp_out_var = out_fh.createVariable("sample", sic_all_extent.dtype, ("sample",))
    time_out_var = out_fh.createVariable("time", sic_all_extent.dtype, ("time",))
    data_out_var = out_fh.createVariable("sic_extent", sic_all_extent.dtype, ("region", "sample", "time",))
    
    # get the time data
    fname = get_sic_fname(run_type, ref_start, ref_end, 1)
    fh = netcdf_file(fname)
    time_out_var._attributes = fh.variables["time"]._attributes
    time_out_var[:] = fh.variables["time"][:]
    samp_out_var[:] = numpy.arange(1,100)
    regi_out_var[:] = numpy.arange(0,2)
    data_out_var[:] = sic_all_extent[:]

    fh.close()
    out_fh.close()

#############################################################################

def plot_sic_extents(sp0, time_data, sic_extent_fname, h=0):
    # load the 99 percentiles of the sic_extent data
    fh = netcdf_file(sic_extent_fname)
    sic_data = fh.variables["sic_extent"][:]
    sp0.fill_between(time_data, sic_data[h,0], sic_data[h,-1], facecolor='r', alpha=0.2, zorder=0)
    for a in range(0, sic_data.shape[1]):
        ldata = sic_data[h,a]
        sp0.plot(time_data, ldata, 'r-', lw=1, alpha=0.5, zorder=0)
    
    fh.close()

#############################################################################

def plot_hadisst_sic_extent(sp0, hemi):
    # load the hadisst_sic
    fname = "/Users/Neil/ClimateData/HadISST2/HadISST.2.1.0.0_realisation_dec2010_400_yrmn.nc"
    ncfh = netcdf_file(fname)
    sic_var = ncfh.variables["sic"]
    lat_var = ncfh.variables["latitude"]
    mv = sic_var._attributes["_FillValue"]
    sic_data = sic_var[:]
    print numpy.max(sic_data)
    lat_data = lat_var[:]
    sic_arctic, sic_antarctic = calc_sea_ice_extent(sic_data, lat_data, 1.0, mv, S=1e-3)
    
    hadisst_x = numpy.arange(1850,2011)
    if hemi == 0:
        sic_extent = sic_arctic
    else:
        sic_extent = sic_antarctic
    sp0.plot(hadisst_x, sic_extent, 'b-', alpha=1.0, lw=2, zorder=2)
    sp0.text(hadisst_x[-1]+2, sic_extent[-1], "HadISST2", color='b')
    ncfh.close()

#############################################################################

def plot_cmip5_sic_extent(sp0, rcp, hemi):
    dirc = "/Users/Neil/Coding/CREDIBLE_output/output/"+rcp+"_2006_2100/sic/"
    if hemi == 0:
        fname = "atlas_sic_OImon_arctic_"+rcp+"_ens_mean_200601-210012_1x1_yrmns.nc"
    else:
        fname = "atlas_sic_OImon_antarctic_"+rcp+"_ens_mean_200601-210012_1x1_yrmns.nc"

    ncfh = netcdf_file(dirc+fname)
    sic_var = ncfh.variables["sic"]
    lat_var = ncfh.variables["latitude"]
    mv = sic_var._attributes["_FillValue"]
    sic_data = numpy.array(sic_var[:])
    sic_data[sic_data < 0] = 0
    lat_data = lat_var[:]
    sic_arctic, sic_antarctic = calc_sea_ice_extent(sic_data, lat_data, 1.0, mv, S=1e-3)
    
    cmip_x = numpy.arange(2006,2101)
    if hemi == 0:
        sic_extent = sic_arctic
    else:
        sic_extent = sic_antarctic
    sp0.plot(cmip_x, sic_extent, 'k-', alpha=1.0, lw=2, zorder=2)
    sp0.text(cmip_x[-1]-5, sic_extent[-1], "CMIP5 MM", color='k')
    ncfh.close()

#############################################################################


if __name__ == "__main__":
    ref_start = -1
    ref_end = -1
    run_type = ""
    eof_year = -1
    monthly = False
    hemi = 0
    opts, args = getopt.getopt(sys.argv[1:], 'r:s:e:f:mh:',
                               ['run_type=', 'ref_start=', 'ref_end=',
                                'monthly', 'hemisphere='])

    for opt, val in opts:
        if opt in ['--run_type', '-r']:
            run_type = val
        if opt in ['--ref_start', '-s']:
            ref_start = int(val)
        if opt in ['--ref_end', '-e']:
            ref_end = int(val)
        if opt in ['--monthly', '-m']:
            monthly = True
        if opt in ['--hemisphere', '-h']:
            hemi = int(val)
            
    sic_extent_fname = get_sic_extent_fname(run_type, ref_start, ref_end)
    if not os.path.exists(sic_extent_fname):
        calc_sic_extent_data(run_type, ref_start, ref_end)
        
    sp0 = plt.subplot(111)
    time_data = numpy.arange(1899,2101)
    sp0.set_xlim([1900,2099])
    if hemi==1:
        sp0.set_ylim([3.5e8,9e8])
        hemistring="Antarctic"
    elif hemi==0:
        sp0.set_ylim([2.5e8,8.5e8])
        hemistring="Arctic"

    plot_sic_extents(sp0, time_data, sic_extent_fname, hemi)
    plot_hadisst_sic_extent(sp0, hemi)
    plot_cmip5_sic_extent(sp0, run_type, hemi)
    
    sp0.set_xlim([1900,2099])
    sp0.set_xlabel("Year")
    sp0.set_ylabel("Sea ice total area, $km^2$")
    sp0.set_title("Total sea ice area, " + hemistring + ", " + run_type.upper())
    outname = "synth_sic_ext_ts_"+hemistring+"_"+run_type+".pdf"
#    plt.gcf().set_size_inches(12, 3, forward=True)
    plt.savefig(outname)
    print outname
