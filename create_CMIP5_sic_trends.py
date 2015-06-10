#! /usr/bin/env python  
#############################################################################
#
# Program : create_CMIP5_sic_trends.py
# Author  : Neil Massey
# Purpose : Create trends from the concatenated CMIP5 sic anomalies over the 
#           two CMIP5 periods: historical (1899->2005) and future RCP (2006->2100)
#           The trend for each timestep (month) is the gradient of the months
#           over 40 years, centred (20th year) on the timestep (month)
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

#############################################################################

if __name__ == "__main__":
    ref_start = -1
    ref_end = -1
    run_type = ""
    start_idx = 0
    end_idx = 0
    opts, args = getopt.getopt(sys.argv[1:], 'r:s:e:i:j:',
                               ['run_type=', 'ref_start=', 'ref_end=', 'st_idx=', 'ed_idx='])

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
            
    create_concat_sic_anoms(run_type, ref_start, ref_end, start_idx, end_idx)