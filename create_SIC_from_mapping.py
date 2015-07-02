#! /usr/bin/env python  
#############################################################################
#
# Program : create_SIC_from_mapping.py
# Author  : Neil Massey
# Purpose : Create a sea-ice concentration file (SIC) from an SST file and
#           a mapping file between SST and SIC.
#           To create the mapping run ./create_HadISST_SST_SIC_mapping.py
#                                  or ./create_CMIP5_SST_SIC_mapping.py
# Inputs  : input   : file of sst values
#           sic_fit : file of fitted sst to sic values
#           output  : output file name
# Output  : in the output directory:
#           output_name
#            
# Date    : 11/06/15
#
#############################################################################

import os, sys, getopt

import numpy
import pyximport
pyximport.install(setup_args={'include_dirs':[numpy.get_include()]})

sys.path.append("../CREDIBLE_SST")
sys.path.append("/Users/Neil/python_lib")
from create_HadISST_SST_SIC_mapping import *
from window_smooth import window_smooth_3D

from netcdf_file import *
import datetime

#############################################################################

def create_SIC_from_mapping(sst_data, time, fit_data, deg_data, mv, anoms=False):
    # create the output for the sic - same shape as the sst input
    sic_out = numpy.zeros(sst_data.shape, 'f')

    bd = datetime.date(1,1,1)

    # loop over every latitude, longitude and month
    for m in range(0, 12):
        sd = bd + datetime.timedelta(time[m])
        for lat in range(0, sst_data.shape[1]):
            for lon in range(0, sst_data.shape[2]):
                # get the sst values - a slice through all the months for the current month
                sst_V = sst_data[m::12,lat,lon]                
                # get the number of degrees for the polyfit
                n_degs = deg_data[m, lat, lon]
                # get the polyfit coefficients
                p_fit = fit_data[:, m, lat, lon]
                if n_degs > 0:
                    # calculate the sic concentration
                    sic_V = calc_polynomial(p_fit, sst_V, n_degs)
                    sic_V[numpy.isnan(sic_V)] = mv
                    if not anoms:
                        sic_V[sic_V > 1.0] = 1.0
                        sic_V[sic_V < 0.0] = 0.0            # 0.0 cutoff to mirror that of 
                                                            # calculating CMIP5 ensemble mean
                    sic_out[m::12, lat, lon] = sic_V
                else:
                    sic_out[m::12, lat, lon] = 0.0

    # remove isolated pixels
#    for t in range(0, sic_out.shape[0]):
#        for lat in range(1, sic_out.shape[1]-1):
#            for lon in range(0, sic_out.shape[2]):
#                if sic_out[t,lat,lon] != 0.0 and sic_out[t,lat-1,lon] == 0.0 and sic_out[t,lat+1,lon] == 0.0:
#                    sic_out[t,lat,lon] = 0.0
    
    # smooth the data with a longitudinal smoother
#    weights = numpy.ones([1,3,1],'f')      # longitudinal smoother
#    smoothed_sic = window_smooth_3D(sic_out, weights, mv, smooth_zero=False)
#    sic_out = smoothed_sic
    return sic_out

#############################################################################

def write_SIC(output, sic_out, lon_var, lat_var, time_var, mv):
    # write out the file
    fh_out = netcdf_file(output, 'w')
    
    # create the dimensions - longitude, latitude, month and polynomial degree
    lon_out_dim = fh_out.createDimension("longitude", lon_var.shape[0])
    lat_out_dim = fh_out.createDimension("latitude", lat_var.shape[0])
    time_out_dim = fh_out.createDimension("time", time_var.shape[0])
    
    # create the variables to go with it
    lon_out_var = fh_out.createVariable("longitude", lon_var[:].dtype, ("longitude",))
    lat_out_var = fh_out.createVariable("latitude", lat_var[:].dtype, ("latitude",))
    time_out_var = fh_out.createVariable("time", time_var[:].dtype, ("time",))
    
    # write out the data - copy the lat / lon data, month and degree are just integer sequences
    lon_out_var[:] = numpy.array(lon_var[:])
    lat_out_var[:] = numpy.array(lat_var[:])
    time_out_var[:] = numpy.array(time_var[:])
    
    # copy the time, lat and lon attributes
    lon_out_var._attributes = lon_var._attributes
    lat_out_var._attributes = lat_var._attributes
    time_out_var._attributes = time_var._attributes
    
    # create the data output
    data_out_var = fh_out.createVariable("sic", 'f', ("time", "latitude", "longitude"))
    data_out_var._attributes["_FillValue"] = mv
    data_out_var._attributes["missing_value"] = mv
    data_out_var._attributes["standard_name"] = "sea_ice_area_fraction"
    data_out_var._attributes["long_name"] = "sea_ice_area_concentration"
    data_out_var._attributes["units"] = "1"
    
    # write the data
    data_out_var[:] = sic_out[:]
    fh_out.close()
    
#############################################################################

def create_SIC_from_mapping_file(input, output, sic_fit):
    # load the sst data, the latitude and the longitude
    fh_sst = netcdf_file(input)
    sst_var = fh_sst.variables["sst"]
    lon_var = fh_sst.variables["longitude"]
    lat_var = fh_sst.variables["latitude"]
    time_var = fh_sst.variables["time"]
    mv = sst_var._attributes["_FillValue"]
    
    # load the mapping data
    fh_sic_map = netcdf_file(sic_fit)
    fit_var = fh_sic_map.variables["polyfit"]
    deg_var = fh_sic_map.variables["degrees"]
    
    sst_data = numpy.array(sst_var[:])
    fit_data = fit_var[:]
    deg_data = deg_var[:]
    
    sic_out = create_SIC_from_mapping(sst_data, fit_data, deg_data, mv)
    write_SIC(output, sic_out, lon_var, lat_var, time_var, mv)
    
    fh_sic_map.close()
    fh_sst.close()

#############################################################################

if __name__ == "__main__":
    in_file = ""
    out_file = ""
    fit_file = ""
    opts, args = getopt.getopt(sys.argv[1:], 'i:f:o:',
                               ['input=', 'sic_fit=', 'output='])

    for opt, val in opts:
        if opt in ['--input', '-i']:
            in_file = val
        if opt in ['--output', '-o']:
            out_file = val
        if opt in ['--sic_fit', '-f']:
            fit_file = val

    create_SIC_from_mapping_file(in_file, out_file, fit_file)