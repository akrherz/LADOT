import netCDF4
import numpy as np
import sys
from pyiem.datatypes import temperature
import datetime

ts = datetime.datetime(2000, int(sys.argv[1]), 1)

cyears = [1971, 1976, 1981, 1986, 1991, 1996]
syears = [2041, 2046, 2051, 2056, 2061, 2066]


def compute_monthly_avg(varname, years, month):
    total = None
    for year in years:
        nc = netCDF4.Dataset(
            "../data/mon_%s_HRM3_hadcm3_%s010100.nc" % (varname, year))
        if total is None:
            total = np.sum(nc.variables[varname][month-1:60:12, :, :], 0)
        else:
            total += np.sum(nc.variables[varname][month-1:60:12, :, :], 0)
        nc.close()

    return total / 30.0


# cmonth = compute_monthly_avg('pr', cyears, ts.month)
# smonth = compute_monthly_avg('pr', syears, ts.month)
cmonth = compute_monthly_avg('tas', cyears, ts.month)
smonth = compute_monthly_avg('tas', syears, ts.month)


delta = temperature(smonth, 'K').value('F') - temperature(cmonth,
                                                          'K').value('F')
# delta = mesonet.k2f(smonth) - mesonet.k2f(cmonth)


# delta = (smonth - cmonth) / cmonth * 100.0
# delta = delta * 10800.0 / 25.4  # inch


nc = netCDF4.Dataset("../data/mon_pr_HRM3_hadcm3_1968010100.nc")
lats = nc.variables['lat'][:]
lons = nc.variables['lon'][:]
nc.close()

nc = netCDF4.Dataset('../data/monthly_deltas.nc', 'a')
nc.variables['deltat'][ts.month - 1, :, :] = delta
if ts.month == 1:
    nc.variables['lat'][:] = lats
    nc.variables['lon'][:] = lons

nc.close()
