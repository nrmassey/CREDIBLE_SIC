#! /usr/bin/env python  
#############################################################################
#
# Program : sst_sic_mapping.py
# Author  : Neil Massey
# Purpose : Create a mapping between the sea-surface temperature (SST) and
#           sea-ice concentration (SIC).  The algorithm requires SST and SIC
#           anomalies.
#           The mapping is created with the following algorithm:
#
#           For every grid box
#               For every month
#                   Remove sample pairs (SST,SIC) where SIC == 0.0
#                   For sample pairs where all SST == 0.0, SIC is the mean 
#                   For the remaining sample pairs
#                       Fit a polynomial of degree d to the sample pairs, where
#                           x = SST and y = SIC
#                           y = p[0]*pow(x,d) + p[1]*pow(x,d-1) + p[d]
#
# Date    : 10/06/15
#
##############################################################################

import numpy
from scipy import interpolate
from scipy.io.netcdf import *

##############################################################################

def calc_sic_mapping(sst_data, sic_data, mv, n_degs=2):
    # calculate the mapping between sst and sic for each grid point
    n_map = sst_data.shape[0] / (10*12) # one mapping per decade
    print sst_data.shape[0]
    n_lat = sst_data.shape[1]
    n_lon = sst_data.shape[2]
    mapping = numpy.zeros([n_map, n_degs+1, 12, n_lat, n_lon], 'f')
    # calculate the mappings per month and grid box
    for m in range(0, 12):
        for lat in range(0, n_lat):
            for lon in range(0, n_lon):
                # get monthly data for this latitude--longitude
                lsst_data = sst_data[m::12,lat,lon]
                lsic_data = sic_data[m::12,lat,lon]
                if lsic_data[0] != mv:
                    # check whether grid box is ice free - don't waste cycles determining the mapping here
                    if numpy.mean(numpy.abs(lsic_data)) < 0.02:
                        mapping[0,:,m,lat,lon] = 0
                        continue
                    # loop over each decade in the month data
                    for d in range(0, n_map):
                        # start decade
                        sd = d*10 - 15
                        if sd < 0:
                            sd = 0
                        # end decade
                        ed = d*10 + 15
                        if ed >= lsst_data.shape[0]:
                            ed = lsst_data.shape[0]
                        # get the data to determine the (linear) mapping
                        sst_map = lsst_data[sd:ed]
                        sic_map = lsic_data[sd:ed]
                        # determine the linear fit between the sst and sic
                        if numpy.sum(numpy.abs(sst_map)) == 0:
                            mapping[d,0:n_degs,m,lat,lon] = 0.0
                            mapping[d,n_degs,m,lat,lon] = numpy.mean(sic_map)
                        else:
                            D,res,rnk,sng,rcd = numpy.polyfit(sst_map, sic_map, n_degs,full=True)
                            if res < 1.0:
                                mapping[d,:,m,lat,lon] = D
                            else:
                                mapping[d,0,m,lat,lon] = 0
                                mapping[d,n_degs,m,lat,lon] = numpy.mean(sic_map)
                else:
                    mapping[:,:,m,lat,lon] = mv
    return mapping

########################################################################################

def calc_sic_from_sst(sst_data, mapping, mv, n_degs=2):
    n_map = mapping.shape[0]
    n_lat = sst_data.shape[1]
    n_lon = sst_data.shape[2]
    out_sic = numpy.zeros(sst_data.shape, 'f')
    X = numpy.zeros([n_map+2])
    X[-1] = n_map * 10
    X[1:-1] = numpy.arange(0,n_map)*10+5
    Sc = float(X[-1]-X[0]+1) / (sst_data.shape[0]/12)
    
    X_s = numpy.arange(X[0], X[-1], Sc)
    # calculate the mappings per month and grid box
    for m in range(0, 12):
        for lat in range(0, n_lat):
            for lon in range(0, n_lon):
                # check for mapping error
                if mapping[0,0,m,lat,lon] == mv:
                    out_sic[m::12,lat,lon] = mv
                    continue
                # get the local sst data
                lsst_data = sst_data[m::12,lat,lon]
                # create the interpolated constants
                new_coeffs = numpy.zeros([n_degs+1, X_s.shape[0]], 'f')
                for c in range(0, n_degs+1):
                    M = mapping[:,c,m,lat,lon]
                    C = numpy.zeros(M.shape[0]+2)
                    C[0] = M[0]
                    C[1:-1] = M
                    C[-1] = M[-1]
                    s = interpolate.interp1d(X, C, kind="linear")
                    new_coeffs[c] = s(X_s)
                # copy the coefficients to last point
                new_coeffs[:,-1] = new_coeffs[:,-2]
                # reconstruct
                V = 0
                for p in range(0, n_degs+1):
                    E = new_coeffs.shape[1]
                    V += (lsst_data[:E])**(n_degs-p) * new_coeffs[p,:]
                if numpy.isnan(V).any():
                    V[:] = mv
                S = V.shape[0]*12
                out_sic[m:S:12,lat,lon] = V
    # fill in the last year as the same as the penultimate year
    for m in range(0,12):
        out_sic[-12+m] = out_sic[-24+m]
    return out_sic
    
########################################################################################

def save_mapping(out_fname, mapping, in_lat_var, in_lon_var, years, mv, n_degs=2):
    # open the file
    out_fh = netcdf_file(out_fname, "w")
    
    # create the output years
    out_years = []
    for y in years:
        out_years.append(int((y[0]+y[1])*0.5))
    
    # create the dimensions - longitude, latitude, month and polynomial degree
    lon_out_dim = out_fh.createDimension("longitude", in_lon_var.shape[0])
    lat_out_dim = out_fh.createDimension("latitude", in_lat_var.shape[0])
    mon_out_dim = out_fh.createDimension("month", 12)
    coeff_out_dim = out_fh.createDimension("coeff", mapping.shape[1])
    year_out_dim = out_fh.createDimension("year", len(out_years))
    
    # create the variables to go with it
    lon_out_var = out_fh.createVariable("longitude", in_lon_var[:].dtype, ("longitude",))
    lat_out_var = out_fh.createVariable("latitude", in_lat_var[:].dtype, ("latitude",))
    mon_out_var = out_fh.createVariable("month", 'i', ("month",))
    coeff_out_var = out_fh.createVariable("coeff", 'i', ("coeff",))
    year_out_var = out_fh.createVariable("year", "i", ("year",))
    
    # write out the data - copy the lat / lon data, month and degree are just integer sequences
    lon_out_var[:] = in_lon_var[:]
    lat_out_var[:] = in_lat_var[:]
    mon_out_var[:] = [x for x in range(0, 12)]
    coeff_out_var[:] = [x for x in range(0, mapping.shape[1])]
    year_out_var[:] = out_years
    
    # copy the lat and lon attributes
    lon_out_var._attributes = in_lon_var._attributes
    lat_out_var._attributes = in_lat_var._attributes
    
    # create the output data variables
    data_out_var = out_fh.createVariable("polyfit", 'f', ("year", "coeff", "month", "latitude", "longitude"))
    # degrees of polyfit
    degrees_out_var = out_fh.createVariable("degrees", 'i', ("year", "month", "latitude", "longitude"))
    
    # attributes
    data_out_var._attributes["_FillValue"] = mv
    
    # data
    data_out_var[:] = mapping[:]
    degrees_out_var[:] = n_degs
    out_fh.close() 
    
#############################################################################

def load_mapping(in_fname):
    in_fh = netcdf_file(in_fname, "r")
    in_var = in_fh.variables["polyfit"]
    in_data = in_var[:]
    mv = in_var._attributes["_FillValue"]
    
    in_fh.close()
    
    return in_data, mv

#############################################################################

def save_sic(output, sic_out, lon_var, lat_var, time_var, mv):
    # write out the file
    fh_out = netcdf_file(output, 'w')
    # create the dimensions - longitude, latitude, month and polynomial degree
    lon_out_dim = fh_out.createDimension("longitude", lon_var.shape[0])
    lat_out_dim = fh_out.createDimension("latitude", lat_var.shape[0])
    time_out_dim = fh_out.createDimension("time", sic_out.shape[0])
    
    # create the variables to go with it
    lon_out_var = fh_out.createVariable("longitude", lon_var[:].dtype, ("longitude",))
    lat_out_var = fh_out.createVariable("latitude", lat_var[:].dtype, ("latitude",))
    time_out_var = fh_out.createVariable("time", time_var[:].dtype, ("time",))
    
    # write out the data - copy the lat / lon data, month and degree are just integer sequences
    lon_out_var[:] = numpy.array(lon_var[:])
    lat_out_var[:] = numpy.array(lat_var[:])
    time_out_var[:] = numpy.array(time_var[:sic_out.shape[0]])
    
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