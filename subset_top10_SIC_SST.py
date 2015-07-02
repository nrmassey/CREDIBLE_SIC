#! /usr/bin/env python  
#############################################################################
#
# Program : subset_top10_SIC_SST.py
# Author  : Neil Massey
# Purpose : Create a subset of the top ten (by SIE RMS score from Q.Shu et 
#           al 2015) members of the CMIP5 ensemble
# Inputs  : 
# Notes   : 
# Output  : 
#            
# Date    : 15/06/15
#
#############################################################################

import os, sys, getopt
from cmip5_functions import get_cmip5_path, get_cmip5_tos_fname, get_cmip5_sic_fname, get_output_directory
from create_CMIP5_sst_anoms import create_remapped_field
from cdo import *

#############################################################################

def get_ens_mean_regrid_string(ens_fnames, var):
    for ens_fn in ens_fnames:
        cdo_string = create_remapped_field(ens_fn, 2006, 2100, var, True)
        print cdo_string

#############################################################################

def get_ensemble_fnames_sic(run_type, model):
    # get all the ensemble members for a particular model for sic
    path = get_cmip5_path()
    path += model + "/" + run_type + "/mon/OImon/"
    ens_dirs = os.listdir(path)
    ens_mems = []
    for dir in ens_dirs:
        ens_mem = os.listdir(path+dir)[0]
        ens_mems.append(path+dir+"/"+ens_mem)
    return ens_mems

#############################################################################

def get_ensemble_fnames_tos(run_type, model):
    # get all the ensemble members for a particular model for sic
    path = get_cmip5_path()
    path += model + "/" + run_type + "/mon/Omon/"
    ens_dirs = os.listdir(path)
    ens_mems = []
    for dir in ens_dirs:
        ens_mem = os.listdir(path+dir)[0]
        ens_mems.append(path+dir+"/"+ens_mem)
    return ens_mems

#############################################################################

def get_ensemble_mean_string(ens_fnames, start_year, end_year):
    strg = ""
    for ens_fn in ens_fnames:
        if str(start_year) in ens_fn and (str(2300) or str(end_year) in ens_fn):
            strg += "-selyear,"+str(start_year)+"/"+str(end_year) + " " +ens_fn + " "
    return strg

#############################################################################

def get_remap_string(var):
    lsm_grid = "/soge-home/staff/coml0118/EXCLUDEFROMBACKUP/LSM/hadisst_grid"
    cdo_string = " -fillmiss "
    if var == "tos":
        cdo_string += " -setctomiss,0 "         # fix for MIROC model
        cdo_string += " -setrtomiss,1e3,1e20 "  # fix for IPSL model
    cdo_string += " -selname,\"" + var + "\""
    return cdo_string, lsm_grid

#############################################################################

def remap(model_name, var):
    sy = 2006
    ey = 2100

    if var == "tos":
        ens_fnames = get_ensemble_fnames_tos(run_type, model)
    elif var == "sic":
        ens_fnames = get_ensemble_fnames_sic(run_type, model)
    
    ens_mean_string = get_ensemble_mean_string(ens_fnames, sy, ey)
    out_dir = get_output_directory(run_type, 1986, 2005)
    run_n = ens_fnames[0].split("/")[9]
    out_ens = out_dir + "/" + ens_fnames[0].split("/")[10].replace(run_n, "ens_mean")
    cdo.ensmean(input=ens_mean_string, output=out_ens)
    remap_string, lsm_grid = get_remap_string(var)
    out_remap = out_ens[:-3] + "_1x1.nc"
    cdo.remapbil(lsm_grid, input=out_ens, output=out_remap)
    print out_remap
    os.remove(out_ens)

#############################################################################

if __name__ == "__main__":
    
    # top 10 Arctic and Antarctic sea-ice models in CMIP5
    top_10_arctic = ["CESM1-BGC", "GFDL-CM3", "CCSM4", "CESM1-CAM5", 
                     "EC-EARTH", "MIROC5", "ACCESS1-3", "CMCC-CMS",
                     "NorESM1-M", "ACCESS1-0"]
    top_10_antarctic = ["CMCC-CM", "NorESM1-M", "MIROC-ESM", 
                        "GISS-E2-H-CC", "ACCESS1-0", "MRI-CGCM3", 
                        "EC-EARTH", "CMCC-CMS", "MIROC-ESM-CHEM", 
                        "bcc-csm1-1"]
                        
    run_type = "rcp45"
    cdo = Cdo()
    
#    for model in top_10_arctic:
#        remap(model, "tos")
#        remap(model, "sic")        
        
    for model in top_10_antarctic:
        if model in top_10_arctic:  # already processed!
            print model
            continue
        remap(model, "tos")
        remap(model, "sic")        
