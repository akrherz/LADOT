import iemplot
import netCDF4
import numpy
import mesonet
import sys
import mx.DateTime

ts = mx.DateTime.DateTime(2000,int(sys.argv[1]), 1)

cyears = [1971, 1976, 1981, 1986, 1991, 1996]
syears = [2041, 2046, 2051, 2056, 2061, 2066]

def compute_monthly_avg(varname, years, month):
    total = None
    for year in years:
        nc = netCDF4.Dataset("data/mon_%s_HRM3_hadcm3_%s010100.nc" % (varname, 
             year))
        if total is None:
            total = numpy.sum( nc.variables[varname][month-1:60:12,:,:], 0)
        else:
            total +=  numpy.sum( nc.variables[varname][month-1:60:12,:,:], 0)
        nc.close()

    return total / 30.0


cmonth = compute_monthly_avg('pr', cyears, ts.month)
smonth = compute_monthly_avg('pr', syears, ts.month)
#cmonth = compute_monthly_avg('tas', cyears, ts.month)
#smonth = compute_monthly_avg('tas', syears, ts.month)


#delta = smonth - cmonth
#delta = mesonet.k2f(smonth) - mesonet.k2f(cmonth)

delta = (smonth - cmonth)  / cmonth * 100.0
delta = delta * 10800.0 / 25.4 # inch

#maxv = numpy.ceil( numpy.max(numpy.abs(delta)) )

nc = netCDF4.Dataset("data/mon_pr_HRM3_hadcm3_1968010100.nc")
lats = nc.variables['lat'][:]
lons = nc.variables['lon'][:]
nc.close()

# Write to summary file
"""
nc = netCDF4.Dataset('monthly_deltas.nc', 'a')
nc.variables['deltat'][ts.month-1,:,:] = delta
#if ts.month == 1:
#	nc.variables['lat'][:] = lats
#	nc.variables['lon'][:] = lons

nc.close()
"""

cfg = {
 'cnLevelSelectionMode': 'ManualLevels',
 'nglSpreadColorStart': -1,
 'nglSpreadColorEnd': 3,
 'cnMinLevelValF': -100,
 'cnMaxLevelValF': 100,
 'cnLevelSpacingF': 5,
 #'cnMinLevelValF': -5,
 #'cnMaxLevelValF': 5,
 #'cnLevelSpacingF': 0.5,
 'wkColorMap': 'posneg_1',
 #'_title': 'HRM3-HADCM3 Temperature Delta (Scenario - Contemp)',
 '_title': 'HRM3-HADCM3 Precipitation Delta (Scenario - Contemp)',
 '_valid': ts.strftime("%B"),
 '_conus': True,
 #'_louisiana': True,
 'lbTitleString': '%',
 #'lbTitleString': 'C',
 '_showvalues': False,
 '_format': '%.0f',
 '_watermark': False,
}

tmpfp = iemplot.simple_grid_fill(lons, lats, delta, cfg)

iemplot.postprocess(tmpfp, None, fname='conus_precip_delta_%s.png' % (ts.strftime("%b"),))
#iemplot.postprocess(tmpfp, None, fname='bogus.png')
