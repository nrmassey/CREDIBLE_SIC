#! /usr/bin/env python  
#############################################################################
#
# Program : create_HadISST_SST_SIC_mapping.py
# Author  : Neil Massey
# Purpose : Create a mapping between the sea-surface temperature (SST) and
#           sea-ice concentration (SIC) for the HadISST2 dataset
#           The mapping is created with the following algorithm:
#
#           For every grid box
#               For every month
#                   For sample pairs (SST,SIC) where SST < -1.79, SIC = 1.0
#                   Remove sample pairs (SST,SIC) where SIC == 0.0
#                   If there are more than 5 sample pairs
#                       Fit a polynomial of degree d to the sample pairs, where
#                           x = SST and y = SIC
#                           y = p[0]*pow(x,d) + p[1]*pow(x,d-1) + p[d]
#
# Inputs  : start : start year to fit the polynomial to
#           end   : end year to fit the polynomial to
#           deg   : degree of the polynomial to fit (recommend 4)
#           rn    : run number
# Output  : in the output directory:
#           hadisst_polyfit_<start>_<end>_<deg>_<rn>.nc
#           a netcdf file containing the polynomial fit for every month and every grid point
#           array format is: [d,m,lat,lon], where d is position in the polynomial above
#               m is month, lat and lon are the grid box
#            
# Date    : 10/06/15
#
#############################################################################

import os, sys, getopt
sys.path.append("../CREDIBLE_SST")
from cmip5_functions import get_output_directory
from create_HadISST_sst_anoms import get_HadISST_input_filename, get_HadISST_output_directory
from netcdf_file import *
import numpy

#############################################################################

def find_polyfit(sst_data, sic_data, max_deg):
    # find the highest polynomial degree which fits the relationship between
    # sst_data and sic_data but doesn't have any inflection points
    fit=False
    deg = max_deg     # maximum degree 4 function
    while not fit:
        pf = numpy.polyfit(sst_data, sic_data, deg)
        R = numpy.max(sst_data)+0.1 - numpy.min(sst_data)-0.1
        S = R/100
        if S < 0.001:
            S = 0.001
        ip = numpy.arange(numpy.min(sst_data)-0.1, numpy.max(sst_data)+0.1, S)
        d2x = calc_polynomial_d2x(pf, ip, deg)
        dx = calc_polynomial_dx(pf, ip, deg)
        if numpy.where(dx < 0.0)[0] != [] and numpy.where(dx > 0.0)[0] != []:
            deg -=1
            if deg == 0:
                fit=True
        else:
            fit=True
    return deg

#############################################################################

def calc_polynomial(pf, ip, deg):
    pp = numpy.zeros(ip.shape[0], 'f')
    for i in range(0, deg+1):
        pp += pf[i] * numpy.power(ip,deg-i)
    return pp

#############################################################################

def calc_polynomial_dx(pf, ip, deg):
    # calculate the first derivative of the polynomial
    pp = numpy.zeros(ip.shape, 'f')
    for i in range(0, deg):
        pp += (deg-i) * pf[i] * numpy.power(ip,deg-i-1)
    return pp

#############################################################################

def calc_polynomial_d2x(pf, ip, deg):
    # calculate the second derivative of the polynomial
    pp = numpy.zeros(ip.shape, 'f')
    for i in range(0, deg-1):
        pp += (deg-i)*(deg-i-1)*pf[i]*numpy.power(ip,deg-i-2)
    return pp

#############################################################################

def save_SST_SIC_mapping(polyfits, degrees, out_fname, in_lat_var, in_lon_var, mv):
    # open the file
    out_fh = netcdf_file(out_fname, "w")
    
    # create the dimensions - longitude, latitude, month and polynomial degree
    lon_out_dim = out_fh.createDimension("longitude", in_lon_var.shape[0])
    lat_out_dim = out_fh.createDimension("latitude", in_lat_var.shape[0])
    mon_out_dim = out_fh.createDimension("month", 12)
    coeff_out_dim = out_fh.createDimension("coeff", polyfits.shape[0])
    
    # create the variables to go with it
    lon_out_var = out_fh.createVariable("longitude", in_lon_var[:].dtype, ("longitude",))
    lat_out_var = out_fh.createVariable("latitude", in_lat_var[:].dtype, ("latitude",))
    mon_out_var = out_fh.createVariable("month", 'i', ("month",))
    coeff_out_var = out_fh.createVariable("coeff", 'i', ("coeff",))
    
    # write out the data - copy the lat / lon data, month and degree are just integer sequences
    lon_out_var[:] = in_lon_var[:]
    lat_out_var[:] = in_lat_var[:]
    mon_out_var[:] = [x for x in range(0, 12)]
    coeff_out_var[:] = [x for x in range(0, polyfits.shape[0])]
    
    # copy the lat and lon attributes
    lon_out_var._attributes = in_lon_var._attributes
    lat_out_var._attributes = in_lat_var._attributes
    
    # create the output data variables
    data_out_var = out_fh.createVariable("polyfit", 'f', ("coeff", "month", "latitude", "longitude"))
    # degrees of polyfit
    degrees_out_var = out_fh.createVariable("degrees", 'i', ("month", "latitude", "longitude"))
    
    # attributes
    data_out_var._attributes["_FillValue"] = mv
    
    # data
    data_out_var[:] = polyfits[:]
    degrees_out_var[:] = degrees[:]
    out_fh.close()  

#############################################################################

def get_HadISST_SST_SIC_mapping_fname(start, end, rn, deg):
    out_dir = get_output_directory("HadISST", start, end)
    out_name = "hadisst_polyfit_"+str(start)+"_"+str(end)+"_"+str(deg)+"_"+str(rn)+".nc"
    return out_dir + "/" + out_name

#############################################################################

def calc_polyfits(sst_var, sic_var, max_deg, start_idx, end_idx, start_lat=0, end_lat=-1, anoms=False):

    # create the storage - [polyfit coefficients, month number, 
    mv = sst_var._attributes["_FillValue"]
    polyfits = numpy.ones([max_deg+1, 12, sst_var.shape[1], sst_var.shape[2]], numpy.float32) * mv
    degrees = numpy.zeros([12, sst_var.shape[1], sst_var.shape[2]], 'i')

    s = start_idx
    e = end_idx
    nm=12

    # reassign end_lat if necessary
    if end_lat == -1:
        end_lat = polyfits.shape[2]

    # loop over the data and perform the polyfit
    for m in range(0, 12):
        for lat in range(start_lat, end_lat):
            for lon in range(0, polyfits.shape[3]):
                # get the data to fit the polynomial to
                sst_data = sst_var[s+m:e+m:nm, lat, lon]
                sic_data = sic_var[s+m:e+m:nm, lat, lon]
                if (sst_data[0] == mv or sic_data[0] == mv):
                    polyfits[0,m,lat,lon] = mv
                    continue
                # adjust the data - first filter for where sic_data is greater than 0.1
                # 0.1 is used as when using the ensemble mean from CMIP5 there is some
                # error around the LSM due to the different LSMs used to construct the
                # ensemble mean
                if not anoms:
                    sic_idx = numpy.where(sic_data > 0.1)
                    sst_data = sst_data[sic_idx]
                    sic_data = sic_data[sic_idx]
                    # now set all sic to 1.0 where sst < -1.79 C
                    sst_idx = numpy.where(sst_data < -1.79 + 273.15)    # need to convert to Kelvins
                    sic_data[sst_idx] = 1.0
                else:
                    # continue if all in the sea ice zero
                    if numpy.mean(numpy.abs(sic_data)) < 0.02:
                        continue

                # only do the polyfit if there is more than 5 years of data
                if sst_data.shape[0] > 5:
                    deg = find_polyfit(sst_data, sic_data, max_deg)
                    polyfits[0:deg+1,m,lat,lon] = numpy.polyfit(sst_data, sic_data, deg)
                    degrees[m,lat,lon] = deg
    return polyfits, degrees, mv

#############################################################################

def create_HadISST_SST_SIC_mapping(start, end, rn, deg):
    # load the HadISST file - get the name from the run number
    hadisst_name = get_HadISST_input_filename(rn)
    nc_fh = netcdf_file(hadisst_name)
    
    # load the sst and sic
    sst_var = nc_fh.variables["sst"]
    sic_var = nc_fh.variables["sic"]
    # load the lat and lon
    lon_var = nc_fh.variables["longitude"]
    lat_var = nc_fh.variables["latitude"]
        
    # get the start and end points in the time axis
    base = 1850    # start of HadISST data
    nm=12
    s = (start-base)*nm
    e = (end-base)*nm

    # get the output filename
    out_name = get_HadISST_SST_SIC_mapping_fname(start, end, rn, deg)

    # do the calculation
    polyfits, degrees, mv = calc_polyfits(sst_var, sic_var, deg, s, e)

    # save the polyfits
    save_SST_SIC_mapping(polyfits, degrees, out_name, lat_var, lon_var, mv)
    nc_fh.close()

#############################################################################

if __name__ == "__main__":
    start = 1978
    end   = 2010
    run_n = 400
    deg   = 4
    opts, args = getopt.getopt(sys.argv[1:], 's:e:n:d:',
                               ['start=', 'end=', 'runn=', 'deg='])

    for opt, val in opts:
        if opt in ['--start', '-s']:
            start = int(val)
        if opt in ['--end', '-e']:
            end = int(val)
        if opt in ['--runn', '-n']:
            run_n = val
        if opt in ['--deg', '-d']:
            deg = int(val)
            
    create_HadISST_SST_SIC_mapping(start, end, run_n, deg)