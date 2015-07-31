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
from scipy import interpolate

from netcdf_file import *
import datetime

#############################################################################

def create_coefficient_splines(sst_sic_map, m, lat, lon, mv):
    # create the coefficient values that has been interpolated using a spline
    # over the time period between the start and end year
    
    # start / end of different datasets
    hadisst_sy = 1899
    cmip5_sy = 2006
    cmip5_ey = 2100

    # years and periods
    years = get_year_intervals()
    
    # number of current coefficients / intervals
    N = len(years)
    # array to hold the years
    y_s = numpy.zeros([N,], 'f')
    # numpy arrays of coefficients - two as we're just doing a linear fit
    coeffs = numpy.zeros([2, N], 'f')

    # counter into the fit coefficients
    y1 = 0
    y2 = 0
    for y in years:
        # calculate the mid year of the period
        y_s[y1] = (years[y1][0] + years[y1][1])*0.5
        # assign the coefficients
        if y[0] < 2006:
            coeffs[:,y1] = sst_sic_map[0][:,y1,m,lat,lon]
        else:
            coeffs[:,y1] = sst_sic_map[2][:,y2,m,lat,lon]
            y2 +=1
        y1+=1

    # if it's the missing value than reset to zero
    coeffs[coeffs==mv] = 0.0
    # now interpolate using a smoothed spline
    y_new = numpy.arange(hadisst_sy, cmip5_ey+1)
    # new coefficients
    new_coeffs = numpy.zeros([2, y_new.shape[0]], 'f')
    for c in range(0, 2):
        s = interpolate.UnivariateSpline(y_s, coeffs[c], s=0)
        new_coeffs[c,:] = s(y_new)
        # overwrite the end and beginning coefficients
        new_coeffs[c,0:6] = coeffs[c,0]
        new_coeffs[c,-6:] = coeffs[c,-1]
    
    return new_coeffs

#############################################################################

def create_SIC_from_mapping(sst_data, time, sst_sic_map, mv, anoms=True):
    # create the output for the sic - same shape as the sst input
    sic_out = numpy.zeros(sst_data.shape, 'f')

    bd = datetime.date(1,1,1)

    # loop over every latitude, longitude and month
    for m in range(0, 12):
        sd = bd + datetime.timedelta(time[m])
        for lat in range(0, sst_data.shape[1]):
            for lon in range(0, sst_data.shape[2]):
                # create the coefficient values for this grid box and month
                cf_vals = create_coefficient_splines(sst_sic_map, m, lat, lon, mv)
                # get the sst values - a slice through all the months for the current month
                sst_V = sst_data[m::12,lat,lon]
                if(sst_V[0] == mv):
                    sic_out[m::12, lat, lon] = mv
                else:
                    sic_out[m::12, lat, lon] = sst_V * cf_vals[0] + cf_vals[1]
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