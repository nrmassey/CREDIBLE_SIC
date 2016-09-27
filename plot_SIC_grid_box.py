import sys
from netcdf_file import *
import matplotlib.pyplot as plt
import numpy

fh = netcdf_file("../CREDIBLE_output/output/HadISST_1899_2005_rcp45_2006_2100_r1986_2005_y2050/varmon/sic/HadISST_1899_2005_rcp45_2006_2100_r1986_2005_y2050_f2050_n6_a50_varmon_sic_mon.nc")
fh_sst = netcdf_file("../CREDIBLE_output/output/HadISST_1899_2005_rcp45_2006_2100_r1986_2005_y2050/varmon/sst/HadISST_1899_2005_rcp45_2006_2100_r1986_2005_y2050_f2050_n6_a50_varmon_ssts_mon.nc")
sic_var = fh.variables["sic"]
sst_var = fh_sst.variables["sst"]
m=3
sic_ts_data = sic_var[m::12,0:60,:]
sst_ts_data = sst_var[m::12,0:60,:]
dates = numpy.arange(1899,2100)
sic_ms_data = numpy.ma.masked_less(sic_ts_data, -1)
sst_ms_data = numpy.ma.masked_less(sst_ts_data, -1)
sp0 = plt.subplot(211)
sp1 = plt.subplot(212)
sp0.plot(dates, sic_ms_data[:,90-68,4:15])
sp1.plot(dates, sst_ms_data[:-1,90-68,4:15])

#for x in range(0, 360):
#    v1 = sic_ms_data[2004-1899,90-68,x:x+1]
#    v2 = sic_ms_data[2010-1899,90-68,x:x+1]
#    if not v1.mask:
#        print x, v2-v1
    
plt.show()
fh.close()
fh_sst.close()