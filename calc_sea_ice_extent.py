#! /usr/bin/env python  
#############################################################################
#
# Program : calc_sea_ice_extent.py
# Author  : Neil Massey
# Purpose : Calculate the total sea ice extent (in km^2) per timestep and
#           per hemisphere (i.e. Arctic and Antarctic) 
# Inputs  : same as ./create_HadISST_CMIP5_syn_SSTs.py
# Output  : in the output directory:
#           output_name
#            
# Date    : 29/06/15
#
#############################################################################

import numpy

earth_r = 6372797 # in metres, mean radius

#############################################################################

def calc_grid_areas(lats, d_lon, r = earth_r):
    """Calculate the surace area of a grid box on a sphere defined by two points, given as two latitude longitude pairs.
Parameters : lat1, lon1 : latitude and longitude of left hand corner
           : lat2, lon2 : latitude and longitude of right hand corner
Returns    : area in metres squared"""

    d_lat = lats[0] - lats[1]
    lats2_rad = numpy.radians(lats) + d_lat * 0.5
    lats1_rad = numpy.radians(lats) - d_lat * 0.5
    d_lon_rad = numpy.radians(d_lon)

    areas = r**2*d_lon_rad * (numpy.sin(lats2_rad) - numpy.sin(lats1_rad))
    return areas

#############################################################################

def calc_sea_ice_extent(sic_data, lats, d_lon, mv, S=1.0):
    # calculate the sea ice extent by multiplying the fractional sea-ice
    # coverage by the grid box area and the summing over the arctic and
    # antarctic separately
    
    # half the number of latitudes
    n_lats2 = lats.shape[0] * 0.5
    # calculate grid areas and reshape to the sic_data shape
    grid_areas = calc_grid_areas(lats, d_lon)
    grid_areas = grid_areas.reshape([1, grid_areas.shape[0], 1])
    
    # mask the sic_data
    sic_data_mv = numpy.ma.masked_equal(sic_data, mv)
    
    # calculate the sic_extent for each grid box
    sic_extent = sic_data_mv * grid_areas
    # calculate the arctic and antarctic SIC by summing over the lat and lon axis
    # after subsetting to the hemisphere
    arctic_extent = numpy.sum(numpy.sum(sic_extent[:,0:n_lats2,:], axis=2), axis=1)
    antarctic_extent = numpy.sum(numpy.sum(sic_extent[:,n_lats2:,:], axis=2), axis=1)
    return arctic_extent*(S**2), antarctic_extent*(S**2)

#############################################################################

def calc_sst_arctic_means(sst_data, lats, d_lon, mv):
    # calculate the "arctic" and "antarctic" temperature means
    # these are actually the upper and lower third of the grid
    S = 0.5
    arcti = lats.shape[0] * S
    anti =  lats.shape[0] - lats.shape[0] * S
    
    grid_areas = calc_grid_areas(lats, d_lon)
    grid_areas_mv = numpy.tile(grid_areas, (sst_data.shape[2],1)).T
    grid_areas_mv[sst_data[0] == mv] = 0.0
    sum_gav_arcti = numpy.sum(numpy.sum(grid_areas_mv[0:arcti], axis=1), axis=0)
    sum_gav_anti = numpy.sum(numpy.sum(grid_areas_mv[anti:], axis=1), axis=0)
    grid_areas = grid_areas.reshape(1, grid_areas.shape[0], 1)

    # mask the sst data
    sst_data_mv = numpy.ma.masked_equal(sst_data, mv)
    print sst_data_mv.shape
    # calculate the sst multiplied by the area for each grid box
    sst_area = sst_data_mv * grid_areas
    # now calc the arctic and antarctic mean SST
    grid_areas_mv = numpy.array(grid_areas)
    arctic_sst = numpy.ma.sum(numpy.ma.sum(sst_area[:,0:arcti,:], axis=2), axis=1) / sum_gav_arcti
    anti_sst = numpy.ma.sum(numpy.ma.sum(sst_area[:,anti:,:], axis=2), axis=1) / sum_gav_anti
    
    return arctic_sst, anti_sst
    
#############################################################################

if __name__ == "__main__":
    lats = [89.5 - x for x in range(0, 180)]
    print lats
    print calc_grid_areas(lats, 1.0)