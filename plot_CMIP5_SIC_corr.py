#! /usr/bin/env python  
#############################################################################
#
# Program : plot_CMIP5_SIC_corr.py
# Author  : Neil Massey
# Purpose : Plot the results of fitting a varying degree polynomial to a grid
#           box in the CMIP5 dataset
# Date    : 10/06/15
#
#############################################################################

import os, sys, getopt, glob
sys.path.append("../CREDIBLE_SST")
from cmip5_functions import get_output_directory
from create_HadISST_sst_anoms import get_HadISST_input_filename
from scipy.io.netcdf import *
import matplotlib.pyplot as plt
import numpy
import scipy.stats as stats

#############################################################################

def get_sic_ens_mean_fname(run_type, model, anoms=False):
    out_dir = get_output_directory(run_type, 1986, 2005)
    fng = out_dir+"/concat_sic_anoms/atlas_sic_OImon_"+model+"_"+run_type+"_ens_mean_??????-??????_1x1"
    fng += "_anoms.nc"
    print fng
    fname = glob.glob(fng)[0]
    return fname

#############################################################################

def get_sst_ens_mean_fname(run_type, model, anoms=False):
    out_dir = get_output_directory(run_type, 1986, 2005)
    fng = out_dir+"/concat_sst_anoms/atlas_tos_Omon_"+model+"_"+run_type+"_ens_mean_??????-??????_1x1"
    fng += "_anoms.nc"
    fname = glob.glob(fng)[0]
    return fname

#############################################################################

def plot_CMIP5_SST_SIC_corr(sp0, run_type, model, yr):
    nm = 12
    m = 9
    d = int((yr-2015)/10)
    
    c = ['ko','bo','go','ro','yo','co','mo','ks','bs','gs','rs','ys']
    lines = []
    mons=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    gx=37
    gy=5

    fname_sic = get_sic_ens_mean_fname(run_type, model, True)
    fname_sst = get_sst_ens_mean_fname(run_type, model, True)
    fh_sic = netcdf_file(fname_sic, 'r')
    fh_sst = netcdf_file(fname_sst, 'r')
    #
    var_sic = fh_sic.variables["sic"]
    var_sst = fh_sst.variables["tos"]
    #
    lsst_data = var_sst[m::12, gy, gx]
    lsic_data = var_sic[m::12, gy, gx]
    sd = d*10 - 15
    if sd < 0:
        sd = 0
    # end decade
    ed = d*10 + 15
    if ed >= lsst_data.shape[0]:
        ed = lsst_data.shape[0]
    # get the data to determine the (linear) mapping
    sst_map = lsst_data[sd:ed].squeeze()
    sic_map = lsic_data[sd:ed].squeeze()
    print sst_map.shape[0]

    D,res,rnk,sng,rcd = numpy.polyfit(sst_map, sic_map, 1,full=True)
    print D,res, rnk
    sst_r = numpy.arange(numpy.min(sst_map),numpy.max(sst_map)+0.01, 0.001)
    sic_r = sst_r*D[0] + D[1]
    sp0.plot(sst_r, sic_r, 'k', lw=2.0)
    sp0.plot(sst_map, sic_map, 'r.', ms=15)

#############################################################################

if __name__ == "__main__":
    
    run_type = "rcp45"
    sp0 = plt.subplot(111)
    yr = 2055
    plot_CMIP5_SST_SIC_corr(sp0, run_type, "arctic", yr)
#    sp0.set_xlim([-1.8,1.8])
    sp0.set_xlabel("SST anomaly $^\circ C$")
    sp0.set_ylabel("SIC anomaly $^\circ C$")
#    sp0.set_ylim([0.0, 1.0])
    plt.tight_layout()
    plt.savefig("SST_SIC_corr.pdf")