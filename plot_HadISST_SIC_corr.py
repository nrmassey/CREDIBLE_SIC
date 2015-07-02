#! /usr/bin/env python  
#############################################################################
#
# Program : plot_HadISST_SIC_corr.py
# Author  : Neil Massey
# Purpose : Plot the results of fitting a varying degree polynomial to a grid
#           box in the HadISST dataset
# Date    : 10/06/15
#
#############################################################################

import os, sys, getopt
sys.path.append("../CREDIBLE_SST")
from create_HadISST_sst_anoms import get_HadISST_input_filename, get_HadISST_year_mean_filename
from create_HadISST_SST_SIC_mapping import *
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

def calc_HadISST_SIC_corr(sp0, rn):
    hadisst_name = get_HadISST_input_filename(rn)
    print hadisst_name
    nc_fh = netcdf_file(hadisst_name)
    sst_var = nc_fh.variables["sst"]
    sic_var = nc_fh.variables["sic"]
    
    start = 1850
    d1 = 1978
    d2 = 2010
    nm=12
    s = (d1-start)*nm
    e = (d2-start)*nm

    stderr = numpy.zeros([sst_var.shape[2], sst_var.shape[1]], 'f')

    c = ['ko','bo','go','ro','yo','co','mo','ks','bs','gs','rs','ys']
    lines = []
    mons=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    for m in range(0, 12):
        sst_data = numpy.array(sst_var[s+m:e+m:nm,14,14] - 273.15) # wonky with deg 4
        sic_data = numpy.array(sic_var[s+m:e+m:nm,14,14])          # wonky
#        sst_data = numpy.array(sst_var[s+m:e+m:nm,33,121] - 273.15)
#        sic_data = numpy.array(sic_var[s+m:e+m:nm,33,121])
        sic_idx = numpy.where(sic_data > 0.0)
        sst_data = sst_data[sic_idx]
        sic_data = sic_data[sic_idx]
        sst_idx = numpy.where(sst_data < -1.75)
        sic_data[sst_idx] = 1.0
        if sst_data.shape[0] > 5:
            deg = find_polyfit(sst_data, sic_data)
#            deg = 4
            pf = numpy.polyfit(sst_data, sic_data, deg)
            R = numpy.max(sst_data)+0.1 - numpy.min(sst_data)-0.1
            S = R/100
            if S < 0.001:
                S = 0.001
            ip = numpy.arange(numpy.min(sst_data)-0.1, numpy.max(sst_data)+0.1, S)
            pp = calc_polynomial(pf, ip, deg)
            pp[pp>1.0] = 1.0
            pp[pp<0.0] = 0.0
            sp0.plot(ip,pp,c[m][0]+"-")
        l = sp0.plot(sst_data,sic_data, c[m])
        lines.append(l[0])
    sp0.legend(lines, mons)
    
#############################################################################

if __name__ == "__main__":
    sp0 = plt.subplot(111)
    hadisst_rns = [1059, 115, 1169, 1194, 1346, 137, 1466, 396, 400, 69]
    calc_HadISST_SIC_corr(sp0, hadisst_rns[0])
    sp0.set_xlabel("SST")
    sp0.set_ylabel("SIC")
    plt.show()