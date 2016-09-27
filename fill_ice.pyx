#cython: boundscheck=False
#cython: wraparound=True

#############################################################################
#
# Program : fill_ice.pyx
# Author  : Neil Massey
# Purpose : Find zero points in the ice data and fill in with the average of 
#           the surrounding 8 (non-zero, non-mv) grid boxes.
#           Input parameters are:
#               data  - field to be smoothed (numpy array)
#               mv    - missing data value (filled data does not take into
#                        account missing data)
# Date    : 26/06/15
#
#############################################################################

cimport numpy
import numpy

cpdef fill_ice(numpy.ndarray[float, ndim=3] data, float mv=2e20):

    # create the storage for the smoothed output
    cdef nt = data.shape[0]         # number of t points
    cdef ny = data.shape[1]         # number of y points
    cdef nx = data.shape[2]         # number of x points
    
    cdef numpy.ndarray[float, ndim=3] out_data = numpy.zeros([nt, ny, nx], 'f')
    
    cdef int t, x, y, i, j          # indices into the array
    cdef int nnz                    # counts

    for t in range(0, nt):
        for y in range(0, ny):
            for x in range(0, nx):
                V = data[t,y,x]
                nnz = 0             # number of surrounding grid points that are not zero
                if V == mv:
                    out_data[t,y,x] = mv
                    continue
                if V >= 0.001:      # fill holes of zero only (rounding error to 3sf here)
                    out_data[t,y,x] = V
                    continue
                    
                sum_V = 0.0
                sum_W = 0.0
                for i in range(-1,2):
                    for j in range(-1,2):
                        if i == 0 and j == 0:
                            continue        # do not compare yrslf!
                        V0 = data[t,y+j,x+i]
                        # get the number of surrounding grid points that are either
                        # less than or greater than the current grid point
                        if V0 >= 0.001 or V0 == mv:
                            nnz += 1
                            if V0 != mv:
                                sum_W += 1
                                sum_V += V0
                # do the averaging if nx >=5 or ny >=5 - i.e. this grid box is
                # surrounded on all but one sides by non-zero or mv grid boxes
                if nnz >= 5 and sum_W != 0:
                    out_data[t,y,x] = sum_V / sum_W
                else:
                    out_data[t,y,x] = 0
    return out_data
