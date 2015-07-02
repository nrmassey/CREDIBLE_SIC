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
from create_HadISST_sst_anoms import get_HadISST_input_filename, get_HadISST_year_mean_filename
from create_HadISST_SST_SIC_mapping import find_polyfit, calc_polynomial
from netcdf_file import *
import matplotlib.pyplot as plt
import numpy
import scipy.stats as stats

#############################################################################

def get_sic_ens_mean_fname(run_type, model, anoms=False):
    out_dir = get_output_directory(run_type, 2006, 2100)
    fng = out_dir+"/sic/atlas_sic_OImon_"+model+"_"+run_type+"_ens_mean_??????-??????_1x1"
    if anoms:
        fng += "_anoms.nc"
    else:
        fng += ".nc"
    fname = glob.glob(fng)[0]
    return fname

#############################################################################

def get_sst_ens_mean_fname(run_type, model, anoms=False):
    out_dir = get_output_directory(run_type, 2006, 2100)
    fng = out_dir+"/tos/atlas_tos_Omon_"+model+"_"+run_type+"_ens_mean_??????-??????_1x1"
    if anoms:
        fng += "_anoms.nc"
    else:
        fng += ".nc"
    fname = glob.glob(fng)[0]
    return fname

#############################################################################

def plot_CMIP5_SST_SIC_corr(sp0, run_type, model_list, sy, ey, anoms=False):
    nm = 12
    
    si = (sy-2006)*nm
    ei = (ey-sy)*nm + si
    
    c = ['ko','bo','go','ro','yo','co','mo','ks','bs','gs','rs','ys']
    lines = []
    mons=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    mod = 0
    for m in range(0, 12):
#    if True:
#        m = 5
        all_sst_data = numpy.zeros([0], 'f')
        all_sic_data = numpy.zeros([0], 'f')
        for model in model_list:
            #
            fname_sic = get_sic_ens_mean_fname(run_type, model, anoms)
            fname_sst = get_sst_ens_mean_fname(run_type, model, anoms)
            #
            fh_sic = netcdf_file(fname_sic, 'r')
            fh_sst = netcdf_file(fname_sst, 'r')
            #
            var_sic = fh_sic.variables["sic"]
            var_sst = fh_sst.variables["tos"]
            #
            if anoms:
                sst_data = numpy.array(var_sst[si+m:ei+m:nm,14,14])
            else:
                sst_data = numpy.array(var_sst[si+m:ei+m:nm,23,15] - 273.15)
            sic_data = numpy.array(var_sic[si+m:ei+m:nm,23,15])
            if numpy.max(sic_data) > 90.0:
                sic_data *= 0.01
            
            if anoms:
                sic_idx = numpy.where((sic_data > -1000.0) & (sst_data < 1000))
            else:
                sst_idx = numpy.where(sst_data < -1.75)
                sic_data[sst_idx] = 1.0
                sic_idx = numpy.where((sic_data > 0.0) & (sst_data < 1000))
            sst_data = sst_data[sic_idx]
            sic_data = sic_data[sic_idx]
            if sst_data.shape[0] != 0:
#                p5 = numpy.percentile(sst_data, 5)
#                p95 = numpy.percentile(sst_data, 95)
#                sst_idx = numpy.where((sst_data > p5) & (sst_data < p95))
#                sst_data = sst_data[sst_idx]
#                sic_data = sic_data[sst_idx]
                # add to all data
                all_sst_data = numpy.append(all_sst_data, sst_data)
                all_sic_data = numpy.append(all_sic_data, sic_data)
            #
            fh_sic.close()
            fh_sst.close()
            mod += 1
            
        # got all the data now fit a polynomial to it
        if all_sst_data.shape[0] > 5:
            deg = find_polyfit(all_sst_data, all_sic_data)
            pf = numpy.polyfit(all_sst_data, all_sic_data, deg)
            R = numpy.max(all_sst_data)+0.1 - numpy.min(all_sst_data)-0.1
            S = R/100
            if S < 0.001:
                S = 0.001
            ip = numpy.arange(numpy.min(all_sst_data)-0.1, numpy.max(all_sst_data)+0.1, S)
            pp = calc_polynomial(pf, ip, deg)
            if not anoms:
                pp[pp>1.0] = 1.0
                pp[pp<0.0] = 0.0
            else:
                pp[pp>1.0] = 1.0
                pp[pp<-1.0] = -1.0
            sp0.plot(ip,pp,c[m][0]+"-")

        l = sp0.plot(all_sst_data,all_sic_data, c[m])


#############################################################################

if __name__ == "__main__":
    
    # top 10 Arctic and Antarctic sea-ice models in CMIP5
    top_10_arctic = ["CESM1-BGC", "GFDL-CM3", "CCSM4", "CESM1-CAM5", 
                     "EC-EARTH", "MIROC5", "ACCESS1-3", "CMCC-CMS",
                     "NorESM1-M", "ACCESS1-0"]
    top_10_arctic = ["arctic"]
    top_10_antarctic = ["CMCC-CM", "NorESM1-M", "MIROC-ESM", 
                        "GISS-E2-H-CC", "ACCESS1-0", "MRI-CGCM3", 
                        "EC-EARTH", "CMCC-CMS", "MIROC-ESM-CHEM", 
                        "bcc-csm1-1"]
    run_type = "rcp45"
    sy = 2010
    ey = 2040
    sp0 = plt.subplot(111)
    plot_CMIP5_SST_SIC_corr(sp0, run_type, top_10_arctic, sy, ey, True)
    sp0.set_xlim([-1.8,1.8])
    sp0.set_xlabel("SST")
    sp0.set_ylabel("SIC")
#    sp0.set_ylim([0.0, 1.0])
    plt.show()