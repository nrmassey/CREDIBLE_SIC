#cython: boundscheck=False
#cython: wraparound=False

#############################################################################
#
# Program : remove_isolated_ice.pyx
# Author  : Neil Massey
# Purpose : Remove ice that has a concentration of zero on the latitude above
#           and below - i.e. isolated ice pockets
#           Input parameters are:
#               data  - field to be smoothed (numpy array)
#               mv    - missing data value (filled data does not take into
#                        account missing data)
# Date    : 26/06/15
#
#############################################################################

cimport numpy
import numpy

cpdef remove_isolated_ice(numpy.ndarray[float, ndim=3] data, float mv=2e20):

    # create the storage for the smoothed output
    cdef nt = data.shape[0]         # number of t points
    cdef ny = data.shape[1]         # number of y points
    cdef nx = data.shape[2]         # number of x points
    
    cdef numpy.ndarray[float, ndim=3] out_data = numpy.zeros([nt, ny, nx], 'f')
    
    cdef int t, x, y, i, j          # indices into the array
    cdef int nnz                    # counts

    out_data[:,0,:] = data[:,0,:]
    out_data[:,ny-1,:] = data[:,ny-1,:]

    for t in range(0, nt):
        for y in range(1, ny-1):
            for x in range(0, nx):
                V = data[t,y,x]
                if V == mv:
                    out_data[t,y,x] = V
                    continue
                V0 = data[t,y-1,x]
                V1 = data[t,y+1,x]
                if V != 0.0 and V0 == 0.0 and V1 == 0.0:
                    out_data[t,y,x] = 0.0
                else:
                    out_data[t,y,x] = data[t,y,x]
    return out_data