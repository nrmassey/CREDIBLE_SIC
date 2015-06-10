#! /usr/bin/env python  
#############################################################################
#
# Program : plot_CMIP5_GMSST_siex_corr.py
# Author  : Neil Massey
# Purpose : Plot the scatter / correlation between the GMSST anomaly and the
#           sea-ice extent anomaly in the CMIP5 ensemble
# Inputs  : run_type  : rcp4.5 | rc8.5 | histo
#           ref_start : year to start reference period, 1850->2005
#           ref_end   : year to end reference period, 1850->2005
#           run_type  : historical | rcp45 | rcp85
# Notes   : all reference values are calculated from the historical run_type
#           anomalies are expressed in terms of this reference value
#           CMIP5 ensemble members are only included if their historical run 
#           includes the reference period
# Output  : in the output/ directory filename is:
#            
# Date    : 15/04/15
#
#############################################################################

import os, sys, getopt
sys.path.append("../CREDIBLE_SST")
from filter_cmip5_members import read_cmip5_index_file
from create_CMIP5_GMT_GMSST_anom_ts import get_gmt_gmsst_anom_ts_fname
from create_CMIP5_siex_ts import get_siex_anom_ts_fname
from netcdf_file import *
import matplotlib.pyplot as plt
import numpy
from scipy.stats import linregress

#############################################################################

def plot_CMIP5_GMSST_siex_corr(run_type, ref_start, ref_end, monthly=False):
    # get the GMT/GMSST anomaly timeseries filename
    gmt_gmsst_fname = get_gmt_gmsst_anom_ts_fname(run_type, ref_start, ref_end, monthly)
    siex_fname = get_siex_anom_ts_fname(run_type, ref_start, ref_end, monthly)

    # load the gmsst data
    fh_gmsst = netcdf_file(gmt_gmsst_fname)
    gmsst = fh_gmsst.variables["tos"][:]
    fh_gmsst.close()
    
    # load the sea ice extent data
    fh_siex = netcdf_file(siex_fname)
    nh_siex = fh_siex.variables["nh_siex"][:]
    sh_siex = fh_siex.variables["sh_siex"][:]
    fh_siex.close()
    
    y = 2050-1899
    
    # calculate the percentile values for the 1st and 99th for the sea ice extent
    p1  = numpy.percentile(nh_siex.flatten(), 1)
    p99 = numpy.percentile(nh_siex.flatten(), 99)
    
    y_gmsst = gmsst[:,y].flatten()
    y_nh_siex = nh_siex[:,y].flatten()
    y_nh_siex = y_nh_siex[y_gmsst < 1e10].flatten()
    y_gmsst = y_gmsst[y_gmsst < 1e10].flatten()
    
    plt.plot(y_gmsst, y_nh_siex, 'k.')
    s, i, r, p, err = linregress(y_gmsst, y_nh_siex)
    p = numpy.array([numpy.min(y_gmsst), numpy.max(y_gmsst)], 'f')
    plt.plot(p, p*s+i, '-', lw=2.0)
    plt.show()
    
    for e in range(0, nh_siex.shape[0]):
        if gmsst[e,0] > 1e10:
            continue
        x = numpy.where((nh_siex[e] < p1))
        if x[0].shape[0] == 0:
#            plt.plot(gmsst[e], nh_siex[e], '.')
            s, i, r, p, err = linregress(gmsst[e], nh_siex[e])
            p = numpy.array([numpy.min(gmsst[e]), numpy.max(gmsst[e])], 'f')
            plt.plot(p, p*s+i, '-', lw=2.0)
    plt.show()

    # calculate the percentile values for the 1st and 99th for the sea ice extent
    p1  = numpy.percentile(sh_siex.flatten(), 1)
    p99 = numpy.percentile(sh_siex.flatten(), 99)
    
    for e in range(0, sh_siex.shape[0]):
        if gmsst[e,0] > 1e10:
            continue
        x = numpy.where((sh_siex[e] < p1))
        if x[0].shape[0] == 0:
#            plt.plot(gmsst[e], sh_siex[e], '.')
            s, i, r, p, err = linregress(gmsst[e], sh_siex[e])
            p = numpy.array([numpy.min(gmsst[e]), numpy.max(gmsst[e])], 'f')
            plt.plot(p, p*s+i, '-', lw=2.0)
    plt.show()

#############################################################################

if __name__ == "__main__":
    ref_start = -1
    ref_end = -1
    run_type = ""
    monthly = False
    opts, args = getopt.getopt(sys.argv[1:], 'r:s:e:m',
                               ['run_type=', 'ref_start=', 'ref_end=',
                                'monthly'])

    for opt, val in opts:
        if opt in ['--run_type', '-r']:
            run_type = val
        if opt in ['--ref_start', '-s']:
            ref_start = int(val)
        if opt in ['--ref_end', '-e']:
            ref_end = int(val)
        if opt in ['--monthly', '-m']:
            monthly = True
            
    plot_CMIP5_GMSST_siex_corr(run_type, ref_start, ref_end, monthly)