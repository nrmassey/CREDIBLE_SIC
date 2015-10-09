
import sys
sys.path.append("../CREDIBLE_SST")
from cmip5_functions import load_data

from create_CMIP5_SST_SIC_mapping import load_CMIP5_anom_data
from create_HadISST_CMIP5_syn_SSTs import get_syn_sst_filename
from create_HadISST_SST_SIC_mapping import load_HadISST_anom_data, load_HadISST_ref_data
from calc_sea_ice_extent import *
import matplotlib.pyplot as plt

def calc_GMSST(fld, lats, mv):
    field = numpy.ma.masked_equal(fld, mv)
    wts = numpy.cos(numpy.deg2rad(lats))
    GMSST = numpy.ma.average(numpy.ma.mean(field, axis=2),axis=1,weights=wts)
    return GMSST

def plot_sic_sst():
    rtype = "rcp45"
    cmip5_sst_anoms, cmip5_sic_anoms, cmip5_mv = load_CMIP5_anom_data(rtype)
    hadisst_sst_anoms, hadisst_sic_anoms, hadisst_mv = load_HadISST_anom_data(400)
    hadisst_sst_ref, hadisst_sic_ref = load_HadISST_ref_data(400)
    
    a = 50
    syn_SSTs_fname = get_syn_sst_filename(rtype, 1986, 2005, 6, 2050, a, 2, True)
    syn_SSTs = load_data(syn_SSTs_fname, "sst")
    syn_sst_mv = numpy.min(syn_SSTs)
    
    syn_SIC_fname = syn_SSTs_fname.replace("ssts", "sic")
    syn_SIC_fname = syn_SIC_fname.replace("sst", "sic")
    syn_SIC_all = load_data(syn_SIC_fname, "sic")
    syn_sic_mv = numpy.min(syn_SIC_all)
    
    m = 0 # month
    cmip5_sub_sic_anoms = cmip5_sic_anoms[m::12]
    hadisst_sub_sic_anoms = hadisst_sic_anoms[m::12]
    syn_sub_sic = syn_SIC_all[m::12]
    lats = numpy.array([90-x for x in range(0,180)])

    hadisst_X = numpy.arange(1850,2011)
    cmip5_X = numpy.arange(2006,2101)
    syn_X = numpy.arange(1899,2100)
    
    # reconstruct from the anomalies
    cmip5_sic = cmip5_sub_sic_anoms + hadisst_sic_ref[m]
    hadisst_sic = hadisst_sub_sic_anoms + hadisst_sic_ref[m]
    cmip5_sic[cmip5_sic < cmip5_mv] = cmip5_mv
    hadisst_sic[hadisst_sic < hadisst_mv] = hadisst_mv

    # calculate sea-ice extent    
    cmip5_sic_arc_extent_anom, cmip5_sic_ant_extent_anom = calc_sea_ice_extent(cmip5_sic, lats, 1.0e-9, cmip5_mv)
    hadisst_sic_arc_extent_anom, hadisst_sic_ant_extent_anom = calc_sea_ice_extent(hadisst_sic, lats, 1.0e-9, hadisst_mv)
    syn_sic_arc_extent, syn_sic_ant_extent = calc_sea_ice_extent(syn_sub_sic, lats, 1.0e-9, syn_sic_mv)

    # plot sea-ice extent
    sp0 = plt.subplot(211)
    sp0.plot(hadisst_X, hadisst_sic_arc_extent_anom, 'k-', lw=2.0)
    sp0.plot(cmip5_X, cmip5_sic_arc_extent_anom, 'r-', lw=1.5)
    sp0.plot(syn_X, syn_sic_arc_extent, "b-", lw=2.0)
    
    hadisst_sst = hadisst_sst_anoms[m::12] + hadisst_sst_ref[m]
    cmip5_sst = cmip5_sst_anoms[m::12] + hadisst_sst_ref[m]
    hadisst_sst[hadisst_sst < hadisst_mv] = hadisst_mv
    cmip5_sst[cmip5_sst < cmip5_mv] = cmip5_mv
    # calc NH GMSST and plot
    cmip5_sst_NH_anoms = calc_GMSST(cmip5_sst[:,:90,:], lats[:90], cmip5_mv)
    hadisst_sst_NH_anoms = calc_GMSST(hadisst_sst[:,:90,:], lats[:90], hadisst_mv)
    syn_sst_anoms = calc_GMSST(syn_SSTs[m::12,:90,:], lats[:90], syn_sst_mv)
    
    sp1 = plt.subplot(212)
    syn_X = numpy.arange(1899,2101)
    sp1.plot(hadisst_X, hadisst_sst_NH_anoms, 'k-', lw=2.0)
    sp1.plot(cmip5_X, cmip5_sst_NH_anoms, 'r-', lw=1.5)
    sp1.plot(syn_X, syn_sst_anoms, 'b-', lw=2.0)
    
    plt.show()

def plot_sic_sst_anom_relationship():
    rtype = "rcp45"
    cmip5_sst_anoms, cmip5_sic_anoms, cmip5_mv = load_CMIP5_anom_data(rtype)
    m = 0
    cmip5_mon_sst = cmip5_sst_anoms[m::12]
    cmip5_mon_sic = cmip5_sic_anoms[m::12]
    
    dec_st = 0
    dec_end = 100
    x = 99
    y = 33
    n_degs = 1
    
    cmip5_sub_sst = cmip5_mon_sst[dec_st:dec_end,y,x]
    cmip5_sub_sic = cmip5_mon_sic[dec_st:dec_end,y,x]
    
    D,res,rnk,sng,rcd = numpy.polyfit(cmip5_sub_sst, cmip5_sub_sic, n_degs,full=True)
    S = D[0] * cmip5_sub_sst + D[1]
    
    sp0 = plt.subplot(111)
    sp0.plot(cmip5_sub_sst, cmip5_sub_sic, 'ro')
    sp0.plot(cmip5_sub_sst, S, 'bo')
    sp0.set_xlabel("sst anom")
    sp0.set_ylabel("sic anom")
    plt.show()

if __name__ == "__main__":
#    plot_sic_sst()
    plot_sic_sst_anom_relationship()