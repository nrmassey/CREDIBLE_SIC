#! /usr/bin/env python
########################################################################################

import os, sys, getopt
sys.path.append("../CREDIBLE_SST")
from cmip5_functions import load_data, get_missing_value
from create_HadISST_SST_SIC_mapping import get_HadISST_monthly_anomaly_filenames, get_HadISST_monthly_ref_filenames
from create_HadISST_sst_anoms import get_HadISST_input_filename, get_HadISST_output_directory 
from calc_sea_ice_extent import *
from scipy import interpolate
import matplotlib.pyplot as plt
from netcdf_file import *

import pyximport
pyximport.install(setup_args={'include_dirs':[numpy.get_include()]})
sys.path.append("/Users/Neil/python_lib")
from window_smooth import window_smooth_3D
from fill_ice import *
from remove_isolated_ice import *
from sst_sic_mapping import *

mon = 11

########################################################################################

def load_hadisst_data():
    rn = 400
    start = 1899
    end = 2010
    sst_fname, sic_fname = get_HadISST_monthly_anomaly_filenames(start, end, rn)
    sst_data = load_data(sst_fname, "sst")
    sic_data = load_data(sic_fname, "sic")
    mv = get_missing_value(sst_fname, "sst")
    
    sst_ref_fname, sic_ref_fname = get_HadISST_monthly_ref_filenames(start, end, rn)
    sic_ref_data = load_data(sic_ref_fname, "sic")
    sst_ref_data = load_data(sst_ref_fname, "sst")
    
    sic_hadisst_fname = get_HadISST_input_filename(rn)
    sic_hadisst = load_data(sic_hadisst_fname, "sic")
    
    return sst_data, sic_data, sic_ref_data, sst_ref_data, sic_hadisst, mv

########################################################################################

def get_HadISST_lon_lat_time_vars():
    rn = 400
    start = 1899
    end = 2010
    sst_fname, sic_fname = get_HadISST_monthly_anomaly_filenames(start, end, rn)

    fh = netcdf_file(sst_fname)
    lon_var = fh.variables["longitude"]
    lat_var = fh.variables["latitude"]
    time_var = fh.variables["time"]
    return lon_var, lat_var, time_var

########################################################################################

if __name__ == "__main__":
    # these are anomalies
    sst_data, sic_data, sic_ref, sst_ref, sic_hadisst, mv = load_hadisst_data()
    
    # this is the full hadisst data
    hadisst_sst_fname = get_HadISST_input_filename(400)
    sst_hadisst = load_data(hadisst_sst_fname, "sst")
    
    lon_var, lat_var, time_var = get_HadISST_lon_lat_time_vars()
    
    # subset to 1850->2010
    had_sic_data = sic_hadisst[:]
    had_sst_data = sst_hadisst[:]
    S = 0
    E = sst_hadisst.shape[0]
    
    years = [[y-5, y+5] for y in range(1850,2010,10)]
            
#    mapping = calc_sic_mapping(sst_data, sic_data, mv)
#    save_mapping("test_map.nc", mapping, lat_var, lon_var, years, mv)

    # calculate the anomaly
    n_rpts = had_sst_data.shape[0] / sst_ref.shape[0]
    had_sst_anom = had_sst_data - numpy.tile(sst_ref, [n_rpts,1,1])

    in_map, mv = load_mapping("test_map.nc")
    
    out_sic = calc_sic_from_sst(had_sst_anom, in_map, mv)
    # add the sic reference back on
    out_sic = out_sic + numpy.tile(sic_ref, [n_rpts,1,1])
    out_sic[out_sic < mv] = mv
    
    print "Filling ice holes"
    sic_filled = fill_ice(out_sic, mv)
    
    print "Removing isolated ice"
    sic_removed = remove_isolated_ice(sic_filled, mv)

    # smooth the ice with a 3x1 smoothing window
    weights = numpy.ones([1,3,1], 'f')
    print "Smoothing"
    sic_smooth = window_smooth_3D(sic_removed, weights, mv, smooth_zero=True)
    # fix range
    sic_smooth[(sic_smooth < 0.0) & (sic_smooth != mv)] = 0.0
    sic_smooth[sic_smooth > 1.0] = 1.0
    
    # remove the block of sea ice in the Baltic sea caused by using the
    # 1986->2005 mean in months 04 to 11
    for m in range(4,12):
        sic_smooth[m::12,28:31,210:213] = 0.0
    
    save_sic("test_sic.nc", sic_smooth, lon_var, lat_var, time_var, mv)

    sic_jan = had_sic_data[mon::12] 
    syn_sic_jan = sic_smooth[mon::12]

    # calculate sea ice extent
    lats = numpy.arange(90,-90,-1)
    d_lon = 1.0
    print numpy.max(syn_sic_jan), numpy.min(syn_sic_jan)
    sie_jan_arctic, sie_jan_antarctic = calc_sea_ice_extent(sic_jan, lats, d_lon, mv)
    syn_jan_arctic, syn_jan_antarctic = calc_sea_ice_extent(syn_sic_jan, lats, d_lon, mv)

    y_s = 1850 + numpy.arange(S/12,E/12)
    
    sp0 = plt.subplot(111)
    print y_s.shape, sie_jan_arctic.shape, syn_jan_arctic.shape
    sp0.plot(y_s, sie_jan_arctic*(1e-3)**2, 'b-')
    sp0.plot(y_s, syn_jan_arctic*(1e-3)**2, 'r-')
    sp0.set_ylim([2.5e8,10.5e8])
    hemistring="Arctic"

#    V = sie_jan_arctic - syn_jan_arctic
#    sp0.plot(y_s, V, 'k--')
#    sp0.set_xlim([y_s[0], y_s[-20]])
#    sp0.set_ylim([sie_jan_arctic[0]*0.9, sie_jan_arctic[-20]*1.1])
    plt.show()