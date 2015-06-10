#! /usr/bin/env python  
#############################################################################
#
# Program : create_CMIP5_sic_anoms.py
# Author  : Neil Massey
# Purpose : Create concatenated CMIP5 sic anomalies over the two CMIP5 periods:
#           historical (1899->2005) and future RCP (2006->2100)
# Inputs  : run_type  : rcp4.5 | rc8.5 | histo
#           ref_start : year to start reference period, 1850->2005
#           ref_end   : year to end reference period, 1850->2005
#           run_type  : historical | rcp45 | rcp85
# Notes   : all reference values are calculated from the historical run_type
#           CMIP5 ensemble members are only included if their historical run 
#           includes the reference period
# Output  : in the output/ directory filename is:
#            
# Date    : 15/04/15
#
#############################################################################

import os, sys, getopt
sys.path.append("../CREDIBLE_SST")
from cmip5_functions import get_output_directory, get_cmip5_sic_fname
from filter_cmip5_members import read_cmip5_index_file, read_cmip5_model_mean_index_file
from create_CMIP5_sst_anoms import get_start_end_periods, create_remapped_field

from cdo import *

#############################################################################

def get_concat_sic_output_path(run_type, ref_start, ref_end):
    out_path = get_output_directory(run_type, ref_start, ref_end)+"/concat_sic_anoms/"
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    return out_path

#############################################################################

def get_concat_anom_sic_output_fname(idx0, idx1, run_type, ref_start, ref_end):
    histo_sy, histo_ey, rcp_sy, rcp_ey = get_start_end_periods()
    out_name = idx0 + "_" + idx1 + "_" +\
               "histo" + "_" + str(histo_sy) + "_" + str(histo_ey) + "_" +\
               run_type + "_" + str(rcp_sy) + "_" + str(rcp_ey) + ".nc"
    out_path = get_concat_sic_output_path(run_type, ref_start, ref_end) + out_name
    return out_path

#############################################################################

def create_concat_sic_anoms(run_type, ref_start, ref_end, start_idx, end_idx, monthly):
    # Build a time series of concatenated sic anomalies (wrt 1986->2005)
    # from 1899->2100
    # get the filtered set of cmip5 models / runs
    cmip5_rcp_idx = read_cmip5_index_file(run_type, ref_start, ref_end)
    n_ens = len(cmip5_rcp_idx)
    if (end_idx == 0):
        end_idx = n_ens

    # create the cdo object
    cdo = Cdo()
    
    # variable name
    sic_var_name = "sic"
    histo_sy, histo_ey, rcp_sy, rcp_ey = get_start_end_periods()
    
    for idx in range(start_idx, end_idx):
        print cmip5_rcp_idx[idx][0]
        # get the rcp and historical sic filenames
        sic_rcp_fname = get_cmip5_sic_fname(run_type, cmip5_rcp_idx[idx][0], cmip5_rcp_idx[idx][1])
        sic_histo_fname = get_cmip5_sic_fname("historical", cmip5_rcp_idx[idx][0], cmip5_rcp_idx[idx][1])
        if sic_rcp_fname == "" or sic_histo_fname == "":
            continue
        # calculate the reference values - these are the 20 year means of each month
        # i.e. January 1986->2005 mean, February 1986->2005 mean
        ref_string = " -ymonmean "
        ref_string += " -selvar," +sic_var_name+\
                      " -selyear,"+str(ref_start)+"/"+str(ref_end)+" "+sic_histo_fname
        # subtract from the historical
        hist_string = " -ymonsub -selvar,"+sic_var_name+\
                      " -selyear,"+str(histo_sy)+"/"+str(histo_ey)+\
                      " "+sic_histo_fname +\
                      " "+ref_string
        # subtract from the projection
        rcp_string = " -ymonsub -selvar,"+sic_var_name+\
                     " -selyear,"+str(rcp_sy)+"/"+str(rcp_ey)+\
                     " "+sic_rcp_fname +\
                     " "+ref_string
        # concatenate the files together and save
        cdo.cat(input=hist_string + " " + rcp_string, output="tmp_sic.nc")
        
        # remap and save to correct filename
        out_path = get_concat_anom_sic_output_fname(cmip5_rcp_idx[idx][0], 
                                                    cmip5_rcp_idx[idx][1], 
                                                    run_type, ref_start, ref_end)

        sic_remap_string = create_remapped_field("tmp_sic.nc", histo_sy, rcp_ey, sic_var_name, False)
        cdo.addc(0,input=sic_remap_string, output="tmp_sic2.nc")
        # set lsm
        lsm_path = "/soge-home/staff/coml0118/LSM/HadISST2_lsm.nc"
        cdo.add(input=" -smooth9 tmp_sic2.nc " + lsm_path, output=out_path)

        os.remove("tmp_sic.nc")
        os.remove("tmp_sic2.nc")
        
#############################################################################

if __name__ == "__main__":
    ref_start = -1
    ref_end = -1
    run_type = ""
    start_idx = 0
    end_idx = 0
    monthly = False
    opts, args = getopt.getopt(sys.argv[1:], 'r:s:e:i:j:m',
                               ['run_type=', 'ref_start=', 'ref_end=', 
                                'st_idx=', 'ed_idx=', 'monthly'])

    for opt, val in opts:
        if opt in ['--run_type', '-r']:
            run_type = val
        if opt in ['--ref_start', '-s']:
            ref_start = int(val)
        if opt in ['--ref_end', '-e']:
            ref_end = int(val)
        if opt in ['--st_idx', '-i']:
            start_idx = int(val)
        if opt in ['--ed_idx', '-j']:
            end_idx = int(val)
        if opt in ['--monthly', '-m']:
            monthly = True
            
    create_concat_sic_anoms(run_type, ref_start, ref_end, start_idx, end_idx, monthly)