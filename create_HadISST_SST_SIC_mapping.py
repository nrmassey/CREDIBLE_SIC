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
from cmip5_functions import load_data
from create_CMIP5_sst_anoms import get_start_end_periods, save_3d_file
from create_HadISST_sst_anoms import get_HadISST_input_filename, get_HadISST_output_directory 
from create_HadISST_sst_anoms import get_HadISST_reference_fname, get_HadISST_annual_cycle_residuals_fname
from create_HadISST_sst_anoms import save_3d_file

import numpy
import pyximport
pyximport.install(setup_args={'include_dirs':[numpy.get_include()]})
sys.path.append("/Users/Neil/python_lib")
from running_gradient_filter import *

from netcdf_file import *

#############################################################################

def get_year_intervals():
#    years = [[1899,1940], [1929,1970], [1959,2000], [1989,2010],
#             [2006,2040], [2039,2070], [2059,2100]]
    years = [[1899+x, 1910+x] for x in range(0, 100,10)]
    years.append([2006,2010])
    for y in range(2009, 2090, 10):
        years.append([y, y+11])
    return years

#############################################################################

def find_polyfit(sst_data, sic_data, max_deg):
    # find the highest polynomial degree which fits the relationship between
    # sst_data and sic_data but doesn't have any inflection points
    fit=False
    deg = max_deg     # maximum degree 4 function
    while not fit:
        try:
            pf = numpy.polyfit(sst_data, sic_data, deg)
        except:
            deg -= 1
            continue
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

def save_SST_SIC_mapping(polyfits, degrees, out_fname, years, start_year, end_year, in_lat_var, in_lon_var, mv):
    # open the file
    out_fh = netcdf_file(out_fname, "w")
    
    # create the output years
    out_years = []
    for y in years:
        if y[0] >= start_year and y[1] <= end_year:
            out_years.append(int((y[0]+y[1])*0.5))
    
    # create the dimensions - longitude, latitude, month and polynomial degree
    lon_out_dim = out_fh.createDimension("longitude", in_lon_var.shape[0])
    lat_out_dim = out_fh.createDimension("latitude", in_lat_var.shape[0])
    mon_out_dim = out_fh.createDimension("month", 12)
    coeff_out_dim = out_fh.createDimension("coeff", polyfits.shape[0])
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
    coeff_out_var[:] = [x for x in range(0, polyfits.shape[0])]
    year_out_var[:] = out_years
    
    # copy the lat and lon attributes
    lon_out_var._attributes = in_lon_var._attributes
    lat_out_var._attributes = in_lat_var._attributes
    
    # create the output data variables
    data_out_var = out_fh.createVariable("polyfit", 'f', ("coeff", "year", "month", "latitude", "longitude"))
    # degrees of polyfit
    degrees_out_var = out_fh.createVariable("degrees", 'i', ("year", "month", "latitude", "longitude"))
    
    # attributes
    data_out_var._attributes["_FillValue"] = mv
    
    # data
    data_out_var[:] = polyfits[:]
    degrees_out_var[:] = degrees[:]
    out_fh.close() 
    print out_fname

#############################################################################

def get_HadISST_SST_SIC_mapping_fname(start, end, rn, deg, anoms=True):
    out_dir = get_HadISST_output_directory(start, end, rn)
    out_name = "hadisst_polyfit_"+str(start)+"_"+str(end)+"_"+str(deg)+"_"+str(rn)
    if anoms:
        out_name += "_anoms.nc"
    else:
        out_name += ".nc"
    return out_dir + "/" + out_name

#############################################################################

def calc_polyfits(sst_var, sic_var, max_deg, mv, start_lat=0, end_lat=-1):

    # create the storage - [polyfit coefficients, month number, 
    polyfits = numpy.ones([max_deg+1, 12, sst_var.shape[1], sst_var.shape[2]], numpy.float32) * mv
    degrees = numpy.zeros([12, sst_var.shape[1], sst_var.shape[2]], 'i')

    nm = 12
    s = 0
    e = sst_var.shape[0] - nm

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
                # continue if all in the sea ice zero
                if numpy.mean(numpy.abs(sic_data)) < 0.02 or\
                   numpy.mean(numpy.abs(sst_data)) < 0.02:
                    continue

                # only do the polyfit if there is more than 5 years of data
                if sst_data.shape[0] > 5:
                    deg = find_polyfit(sst_data, sic_data, max_deg)
                    polyfits[0:deg+1,m,lat,lon] = numpy.polyfit(sst_data, sic_data, deg)
                    degrees[m,lat,lon] = deg
    return polyfits, degrees

#############################################################################

def get_HadISST_monthly_anomalies(rn):
    # load the HadISST file expressed as anomalies from monthly means of the
    # period 1986->2005 (ymonmean)
    # load the HadISST file - get the name from the run number
    hadisst_name = get_HadISST_input_filename(rn)
    nc_fh = netcdf_file(hadisst_name)
    
    # load the sst and sic
    sst_var = nc_fh.variables["sst"]
    sic_var = nc_fh.variables["sic"]
    # load the lat and lon
    lon_var = nc_fh.variables["longitude"]
    lat_var = nc_fh.variables["latitude"]
    mv = nc_fh.variables["sst"]._attributes["_FillValue"]

    # get the sst and sic reference
    histo_sy = 1899
    histo_ey = 2010
    ref_start = 1986
    ref_end = 2005
    hadisst_sst_ref_fname = get_HadISST_reference_fname(histo_sy, histo_ey, ref_start, ref_end, rn)
    hadisst_sic_ref_fname = hadisst_sst_ref_fname[:-3] + "_sic.nc"
    
    # get the sst annual cycle
    hadisst_ac_ref_fname = get_HadISST_annual_cycle_residuals_fname(histo_sy, histo_ey, ref_start, ref_end, rn)
    
    # load the data
    hadisst_sst_ref = load_data(hadisst_sst_ref_fname, "sst")
    hadisst_sic_ref = load_data(hadisst_sic_ref_fname, "sic")
    hadisst_ac_ref  = load_data(hadisst_ac_ref_fname, "sst")
    hadisst_sst     = load_data(hadisst_name, "sst")
    hadisst_sic     = load_data(hadisst_name, "sic")

    # subset the hadisst sst and sic data
    hadisst_sst = hadisst_sst[:].byteswap().newbyteorder()
    hadisst_sic = hadisst_sic[:].byteswap().newbyteorder()
        
    # tile the sst ac residuals and subtract from the sst. then subtract the reference pattern
    n_rpts = hadisst_sst.shape[0] / 12
    hadisst_ac_tile = numpy.tile(hadisst_ac_ref, [n_rpts,1,1])
    hadisst_sst_anoms = hadisst_sst - hadisst_ac_tile - hadisst_sst_ref

    # tile the sic reference and subtract from the sic
    hadisst_sic_ref_tile = numpy.tile(hadisst_sic_ref, [n_rpts,1,1])
    hadisst_sic_anoms = hadisst_sic - hadisst_sic_ref_tile
    
    # restore the lsm
    hadisst_sst_anoms[hadisst_sst==mv] = mv
    hadisst_sic_anoms[hadisst_sst==mv] = mv
    
    # smooth the sea-ice and sst
    smooth = False
    if smooth:
        P = 40
        smoothed_sst_anoms = running_gradient_3D(hadisst_sst_anoms, P, mv)
        smoothed_sic_anoms = running_gradient_3D(hadisst_sic_anoms, P, mv)
    else:
        smoothed_sst_anoms = hadisst_sst_anoms
        smoothed_sic_anoms = hadisst_sic_anoms
    
    return smoothed_sst_anoms, smoothed_sic_anoms

#############################################################################

def create_HadISST_SST_SIC_mapping(rn, deg):
    # get the latitude, longitude and missing value
    hadisst_name = get_HadISST_input_filename(rn)
    nc_fh = netcdf_file(hadisst_name)
    # load the lat and lon and variable definititions 
    lon_var = nc_fh.variables["longitude"]
    lat_var = nc_fh.variables["latitude"]
    # get the mv
    mv = nc_fh.variables["sst"]._attributes["_FillValue"]

    # get the sic and sst anomalies
    hadisst_sst_anoms, hadisst_sic_anoms = get_HadISST_monthly_anomalies(rn)

    # test - save the anomaly files out
#    sst_var = nc_fh.variables["sst"]
#    sic_var = nc_fh.variables["sic"]
#    time_var = nc_fh.variables["time"]
#    save_3d_file("test_sst.nc", hadisst_sst_anoms, lon_var, lat_var, sst_var._attributes, time_var)
#    save_3d_file("test_sic.nc", hadisst_sic_anoms, lon_var, lat_var, sic_var._attributes, time_var)

    # get the start and end periods
    years = get_year_intervals()
    hadisst_by = 1850
    hadisst_ey = 2010

    n_years = 0
    for y in years:
        if y[1] <= hadisst_ey:
            n_years += 1
    
    # create the output
    output_degrees  = numpy.zeros([n_years, 12, lat_var.shape[0], lon_var.shape[0]], 'f')
    output_polyfits = numpy.zeros([deg+1, n_years, 12, lat_var.shape[0], lon_var.shape[0]], 'f')

    yi = 0
    # loop over each year
    for year in years:
        if year[1] > hadisst_ey:
            break
        # calculate the start and end indices
        start = year[0]
        end = year[1]
        si = (start-hadisst_by)*12
        ei = (end-hadisst_by)*12
        
        # extrapolate beyond the 10 years to ensure a good continuity
        if si > 120:
            si -= 120
        if ei < hadisst_sst_anoms.shape[0] - 120:
            ei += 120
        
        # subset the data
        hadisst_sst_anoms_local = hadisst_sst_anoms[si:ei]
        hadisst_sic_anoms_local = hadisst_sic_anoms[si:ei]

        # do the calculation
        polyfits, degrees = calc_polyfits(hadisst_sst_anoms_local, hadisst_sic_anoms_local, deg, mv)

        # assign to outputs
        output_polyfits[:,yi,:,:,:] = polyfits
        output_degrees[yi] = degrees
        yi += 1

    # get the output filename
    out_name = get_HadISST_SST_SIC_mapping_fname(years[0][0], hadisst_ey, rn, deg)
    # save the polyfits
    save_SST_SIC_mapping(output_polyfits, output_degrees, out_name, years, years[0][0], hadisst_ey, lat_var, lon_var, mv)
        
    nc_fh.close()

#############################################################################

if __name__ == "__main__":
    run_n = 400
    deg   = 4
    opts, args = getopt.getopt(sys.argv[1:], 'n:d:',
                               ['runn=', 'deg='])

    for opt, val in opts:
        if opt in ['--runn', '-n']:
            run_n = val
        if opt in ['--deg', '-d']:
            deg = int(val)
            
    create_HadISST_SST_SIC_mapping(run_n, deg)