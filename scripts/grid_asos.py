# Grid hourly ASOS data please
# Temperature (K),
# Wind speed (mps),
# % Sun shine,
# Precipitation,
# Relative humidity.

import netCDF4
import numpy as np
import datetime
import psycopg2.extras
from scipy.interpolate import griddata
from pyiem.datatypes import temperature
from pyiem import meteorology
from pyiem.network import Table as NetworkTable

asos = psycopg2.connect(database='asos', host='iemdb', user='nobody')
acursor = asos.cursor(cursor_factory=psycopg2.extras.DictCursor)

badtimes = [
 ]
BASE = datetime.datetime(1970, 1, 1)
STOP = datetime.datetime(2012, 1, 1)
sts = datetime.datetime(2010, 1, 1)
ets = datetime.datetime(2012, 1, 1)

# -94.043341 28.926454,-88.816516 33.019544
WEST = -94.25
SOUTH = 28.75
EAST = -88.75
NORTH = 33.25
NX = 22
NY = 22
DELTAX = (EAST - WEST) / float(NX)
DELTAY = (NORTH - SOUTH) / float(NY)
XAXIS = WEST + DELTAX * np.arange(0, NX)
YAXIS = SOUTH + DELTAY * np.arange(0, NY)
XI, YI = np.meshgrid(XAXIS, YAXIS)


def create_netcdf():
    nc = netCDF4.Dataset("data/asosgrid.nc", 'w')
    # Dimensions
    nc.createDimension('time', int((STOP - BASE).hours))
    nc.createDimension('lat', int((NORTH-SOUTH)/DELTAY))
    nc.createDimension('lon', int((EAST-WEST)/DELTAX))
    # Variables
    tm = nc.createVariable('time', 'd', ('time',))
    tm.units = 'hours since 1970-01-01'
    tm[:] = range(int((STOP - BASE).total_seconds() / 3600.))

    lat = nc.createVariable('lat', 'd', ('lat',))
    lat.units = 'degrees north'
    lat.long_name = 'Latitude'
    lat.axis = 'Y'
    lat[:] = np.arange(SOUTH, NORTH, DELTAY)

    lon = nc.createVariable('lon', 'd', ('lon',))
    lon.units = 'degrees east'
    lon.long_name = 'Longitude'
    lon.axis = 'X'
    lon[:] = np.arange(WEST, EAST, DELTAX)

    tmpk = nc.createVariable('tmpk', 'f', ('time', 'lat', 'lon'))
    tmpk.units = 'K'
    tmpk._FillValue = 1.e20
    tmpk.missing_value = 1.e20
    tmpk.long_name = 'Surface Air Temperature'

    smps = nc.createVariable('smps', 'f', ('time', 'lat', 'lon'))
    smps.units = 'm s-1'
    smps._FillValue = 1.e20
    smps.missing_value = 1.e20
    smps.long_name = '10m Wind Speed'

    skyc = nc.createVariable('skyc', 'f', ('time', 'lat', 'lon'))
    skyc.units = '%'
    skyc._FillValue = 1.e20
    skyc.missing_value = 1.e20
    skyc.long_name = 'Sky Coverage'

    p01m = nc.createVariable('p01m', 'f', ('time', 'lat', 'lon'))
    p01m.units = 'mm'
    p01m._FillValue = 1.e20
    p01m.missing_value = 1.e20
    p01m.long_name = 'Precipitation'

    relh = nc.createVariable('relh', 'f', ('time', 'lat', 'lon'))
    relh.units = '%'
    relh._FillValue = 1.e20
    relh.missing_value = 1.e20
    relh.long_name = 'Relative Humidity'

    nc.close()


def grid_wind(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        if rs[i]['max_sknt'] is not None:
            lats.append(nt.sts[rs[i]['station']]['lat'])
            lons.append(nt.sts[rs[i]['station']]['lon'])
            vals.append(rs[i]['max_sknt'] * 0.514)  # knots to mps
    if len(vals) < 4:
        print "No WIND data at all for time: %s" % (ts,)
        return
    grid = griddata((lons, lats), np.array(vals), (XI, YI), method='nearest')
    offset = int((ts - BASE).total_seconds() / 3600.)
    if grid is not None:
        nc.variables['smps'][offset, :, :] = np.where(grid < 0., 0., grid)
    else:
        print "WIND gridding failed, len vals %s" % (len(vals),)


def grid_relh(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        if rs[i]['max_tmpf'] is not None and rs[i]['max_dwpf'] is not None:
            lats.append(nt.sts[rs[i]['station']]['lat'])
            lons.append(nt.sts[rs[i]['station']]['lon'])
            vals.append((meteorology.relh(temperature(rs[i]['max_tmpf'], 'F'),
                                          temperature(rs[i]['max_dwpf'], 'F'))
                         ).value('%'))
    if len(vals) < 4:
        print "No RELH data at all for time: %s" % (ts,)
        return
    grid = griddata((lons, lats), np.array(vals), (XI, YI), method='nearest')
    offset = int((ts - BASE).total_seconds() / 3600.)
    if grid is not None:
        nc.variables['relh'][offset, :, :] = np.where(grid < 12., 12., grid)
    else:
        print "RELH gridding failed, len vals %s" % (len(vals),)


def grid_skyc(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        v = max(rs[i]['max_skyc1'], rs[i]['max_skyc2'], rs[i]['max_skyc3'])
        if v is not None:
            lats.append(nt.sts[rs[i]['station']]['lat'])
            lons.append(nt.sts[rs[i]['station']]['lon'])
            vals.append(float(v))
    if len(vals) < 4:
        print "No SKYC data at all for time: %s" % (ts,)
        return
    grid = griddata((lons, lats), np.array(vals), (XI, YI), method='nearest')
    offset = int((ts - BASE).total_seconds() / 3600.)
    if grid is not None:
        gt = np.where(grid > 0., grid, 0.0)
        nc.variables['skyc'][offset, :, :] = np.where(gt > 100., 100., gt)
    else:
        print "SKYC gridding failed, len vals %s" % (len(vals),)
        print vals


def grid_tmpf(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        if rs[i]['max_tmpf'] is not None:
            lats.append(nt.sts[rs[i]['station']]['lat'])
            lons.append(nt.sts[rs[i]['station']]['lon'])
            vals.append(temperature(rs[i]['max_tmpf'], 'F').value('K'))
    if len(vals) < 4:
        print "No TMPF data at all for time: %s" % (ts,)
        return
    grid = griddata((lons, lats), np.array(vals), (XI, YI))
    offset = int((ts - BASE).total_seconds() / 3600.)
    if grid is not None:
        nc.variables['tmpk'][offset, :, :] = grid
    else:
        print "TMPK gridding failed, len vals %s" % (len(vals),)


def grid_p01m(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        lats.append(nt.sts[rs[i]['station']]['lat'])
        lons.append(nt.sts[rs[i]['station']]['lon'])
        vals.append(rs[i]['max_p01m'])

    grid = griddata((lons, lats), np.array(vals), (XI, YI), method='nearest')
    offset = int((ts - BASE).total_seconds() / 3600.)
    if grid is not None:
        nc.variables['p01m'][offset, :, :] = np.where(grid > 0., grid, 0.)
    else:
        print "P01M gridding failed, len vals %s" % (len(vals),)


def grid_hour(nc, ts, ids):

    sql = """SELECT station,
         max(case when tmpf > -60 and tmpf < 130
         THEN tmpf else null end) as max_tmpf,
         max(case when sknt > 0 and sknt < 100
         then sknt else 0 end) as max_sknt,
         max(getskyc(skyc1)) as max_skyc1,
         max(getskyc(skyc2)) as max_skyc2,
         max(getskyc(skyc3)) as max_skyc3,
         max(case when p01i > 0 and p01i < 1000
         then p01i * 25.4 else 0 end) as max_p01m,
         max(case when dwpf > -60 and dwpf < 100
         THEN dwpf else null end) as max_dwpf from alldata
         WHERE station in %s and
         valid >= '%s' and valid < '%s' GROUP by station
         """ % (ids, ts.strftime("%Y-%m-%d %H:%M"),
                (ts + datetime.timedelta(hours=1)).strftime("%Y-%m-%d %H:%M"))
    acursor.execute(sql)
    if acursor.rowcount > 4:
        rows = []
        for row in acursor:
            rows.append(row)
        grid_tmpf(nc, ts, rows)
        grid_relh(nc, ts, rows)
        grid_wind(nc, ts, rows)
        grid_skyc(nc, ts, rows)
        grid_p01m(nc, ts, rows)
    else:
        print(("%s has %02i entries, FAIL"
               ) % (ts.strftime("%Y-%m-%d %H:%M"), acursor.rowcount))

#
# create_netcdf()
# sys.exit()
# Include Florida data so that we get some analysis over the gulf
nt = NetworkTable(('MS_ASOS', 'LA_ASOS', 'AR_ASOS',
                   'TX_ASOS', 'OK_ASOS', 'FL_ASOS'))
ids = str(nt.sts.keys())
ids = "(%s)" % (ids[1:-1],)
nc = netCDF4.Dataset("../data/asosgrid.nc", 'a')
# """
# Redo 31 Dec each year
for yr in range(1970, 2011):
    for hr in range(24):
        now = datetime.datetime(yr, 12, 31, hr)
        print now
        grid_hour(nc, now, ids)
"""
now = sts
while now < ets:
    print now
    grid_hour(nc, now)
    now += datetime.timedelta(hours=1)
"""
nc.close()
