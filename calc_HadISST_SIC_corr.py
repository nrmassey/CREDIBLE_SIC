#! /usr/bin/env python  
#############################################################################
#
# Program : calc_HadISST_SIC_corr.py
# Author  : Neil Massey
# Purpose : Calculate the residuals if a linear regression is performed between
#           the sea-ice coverage (SIC) and the sea-surface temperature (SST)
#           in the HadISST SST dataset
# Inputs  : 
# Notes   : 
# Output  : 
# Date    : 14/04/15
#
#############################################################################

import os, sys, getopt
sys.path.append("../CREDIBLE_SST")
from create_HadISST_sst_anoms import get_HadISST_input_filename, get_HadISST_year_mean_filename
from netcdf_file import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import scipy.stats
import numpy
import cartopy.crs as ccrs
import matplotlib.colors as col

#############################################################################

def create_color_map():
    # this is the color scale from Neu et al 2013 (BAMS)
    levels = numpy.array([0,0.0001,0.001,0.01,0.1,1.0,2.0,3.0,4.0,5.0,6.0], 'f')
    cmap = ["#ffffff", "#f8e73d", "#e9a833", "#009a4a",
            "#00aeca", "#0088b7", "#295393", "#2e3065", 
            "#973683", "#ff0000", "#999999"]
    ccmap, norm = col.from_levels_and_colors(levels, cmap, 'max')
    return ccmap, norm, levels

#############################################################################

def calc_HadISST_SIC_corr():
    hadisst_name = get_HadISST_input_filename(400)
    nc_fh = netcdf_file(hadisst_name)
    sst_var = nc_fh.variables["sst"]
    sic_var = nc_fh.variables["sic"]
    
    start = 1850
    d1 = 1850
    d2 = 2010
    nm=12
    s = (d1-start)*nm
    e = (d2-start)*nm

    stderr = numpy.zeros([sst_var.shape[2], sst_var.shape[1]], 'f')

    sst_data = sst_var[s:e:nm,155,177] - 273.15
    sic_data = sic_var[s:e:nm,155,177]
    print s,e,nm

#    plt.plot(sic_data,'k-')
#    plt.plot([x for x in range(0,sic_data.shape[0],nm)],[(0.0,1.0) for x in range(0,sic_data.shape[0],nm)])
#    plt.show()
    plt.plot(sst_data,sic_data,'ro')
    plt.show()
    sys.exit()

    for j in range(stderr.shape[1]):
        for i in range(0, stderr.shape[0]):
            sst_data = sst_var[s:e,j,i]
            sic_data = sic_var[s:e,j,i]
            if abs(sst_data[0]) > 1e3:
                stderr[i,j] = 0.0
            else:
                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(sst_data,sic_data)
                if numpy.isfinite(std_err):
                    stderr[i,j] = std_err
                else:
                    std_err = 0.0

    X = nc_fh.variables["longitude"][:]
    Y = nc_fh.variables["latitude"][:]
    
    cmap, norm, levels = create_color_map()
    
    projection = ccrs.PlateCarree()
    sp = plt.subplot(111, projection=projection)
    X = numpy.append(X, X[-1]+X[-1]-X[-2])
    stderr = numpy.append(stderr, stderr[0,:])
    stderr = stderr.reshape((X.shape[0], Y.shape[0]))
    
    sp.pcolormesh(X, Y, stderr.T, cmap=cmap, norm=norm,
                  transform=ccrs.PlateCarree())

    sp.set_aspect(1.0)
    sp.coastlines(lw=1.0)
    sp.gridlines()

    fig = plt.gcf()
    cax = fig.add_axes([0.90, 0.2, 0.02, 0.6])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
                                   spacing='uniform', extend='both')
    
    plt.show()
        
    nc_fh.close()

#############################################################################

if __name__ == "__main__":
    calc_HadISST_SIC_corr()