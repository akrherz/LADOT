# Grid daily COOP data please
# high
# low
# precip

import sys
import netCDF3
import numpy
import mx.DateTime
from pyIEM import iemdb, mesonet
import Ngl
i = iemdb.iemdb()
mesosite = i['mesosite']
coop = i['coop']
locs = {}

badtimes = []
BASE = mx.DateTime.DateTime(1970,1,1)
STOP = mx.DateTime.DateTime(2012,1,1)
sts = mx.DateTime.DateTime(1970,1,1)
ets = mx.DateTime.DateTime(2012,1,1)

# -94.043341 28.926454,-88.816516 33.019544
WEST  = -94.25
SOUTH = 28.75
EAST  = -88.75
NORTH =  33.25
NX = 22
NY = 22
DELTAX = (EAST - WEST) / float(NX) 
DELTAY = (NORTH - SOUTH) / float(NY)
XAXIS = WEST + DELTAX * numpy.arange(0, NX)
YAXIS = SOUTH + DELTAY * numpy.arange(0, NY)


def load_stationtable():
    rs = mesosite.query("""SELECT id, x(geom) as lon, y(geom) as lat from
         stations where network = 'LACLIMATE'""").dictresult()
    for i in range(len(rs)):
        locs[ rs[i]['id'] ] = rs[i]
        locs[rs[i]['id'] ]['lat'] += (0.00001 * i) # Jitter

def reset_precip(yr):
    # Check precip
    for stationid in locs.keys():
        rs2 = coop.query("""SELECT sum(precip) from alldata_la WHERE
          station = '%s' and year = %s""" % (stationid, yr)).dictresult()
        if rs2[0]['sum'] > 20.:
            locs[ stationid ]['iprecip'] = True
        else:
            locs[ stationid ]['iprecip'] = False

def create_netcdf():
    nc = netCDF3.Dataset("data/coopgrid.nc", 'w')
    # Dimensions
    nc.createDimension('time', 0)
    nc.createDimension('lat', int((NORTH-SOUTH)/DELTAY) )
    nc.createDimension('lon', int((EAST-WEST)/DELTAX)   )
    # Variables
    tm = nc.createVariable('time', 'd', ('time',))
    tm.units = 'days since 1970-01-01'
    tm[:] = range( int((STOP - BASE).days) )

    lat = nc.createVariable('lat', 'd', ('lat',))
    lat.units = 'degrees north'
    lat.long_name = 'Latitude'
    lat.axis = 'Y'
    lat[:] = numpy.arange(SOUTH, NORTH, DELTAY)

    lon = nc.createVariable('lon', 'd', ('lon',))
    lon.units = 'degrees east'
    lon.long_name = 'Longitude'
    lon.axis = 'X'
    lon[:] = numpy.arange(WEST, EAST, DELTAX)

    high = nc.createVariable('high', 'f', ('time','lat','lon'))
    high.units = 'K'
    high._FillValue = 1.e20
    high.missing_value = 1.e20
    high.long_name = 'Daily High Temperature'

    low = nc.createVariable('low', 'f', ('time','lat','lon'))
    low.units = 'K'
    low._FillValue = 1.e20
    low.missing_value = 1.e20
    low.long_name = 'Daily Low Temperature'

    pday = nc.createVariable('p01d', 'f', ('time','lat','lon'))
    pday.units = 'mm'
    pday._FillValue = 1.e20
    pday.missing_value = 1.e20
    pday.long_name = 'Daily Precipitation'

    nc.close()

def grid_high(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        if rs[i]['high'] is not None:
            lats.append(  locs[rs[i]['station']]['lat'] )
            lons.append(  locs[rs[i]['station']]['lon'] )
            vals.append( mesonet.f2k( rs[i]['high'] ) )
    grid = Ngl.natgrid(lons, lats, vals, XAXIS, YAXIS)
    offset = int((ts - BASE).days ) 
    if grid is not None:
        nc.variables['high'][offset,:,:] = grid.transpose()
    else:
        print "HIGH gridding failed, len vals %s" % (len(vals),)

def grid_low(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        if rs[i]['low'] is not None:
            lats.append(  locs[rs[i]['station']]['lat'] )
            lons.append(  locs[rs[i]['station']]['lon'] )
            vals.append( mesonet.f2k( rs[i]['low'] ) )
    grid = Ngl.natgrid(lons, lats, vals, XAXIS, YAXIS)
    offset = int((ts - BASE).days ) 
    if grid is not None:
        nc.variables['low'][offset,:,:] = grid.transpose()
    else:
        print "LOW gridding failed, len vals %s" % (len(vals),)


def grid_p01d(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        if not locs[rs[i]['station']]['iprecip'] or rs[i]['p01d'] is None:
            continue
        lats.append(  locs[rs[i]['station']]['lat'] )
        lons.append(  locs[rs[i]['station']]['lon'] )
        vals.append( rs[i]['p01d'] )
    grid = Ngl.natgrid(lons, lats, vals, XAXIS, YAXIS)
    offset = int((ts - BASE).days )
    if grid is not None:
        gt = grid.transpose()
        nc.variables['p01d'][offset,:,:] = numpy.where(gt > 0., gt, 0.)
    else:
        print "P01D gridding failed, len vals %s" % (len(vals),)

def grid_hour(nc, ts):
    ids = `locs.keys()`
    ids = "(%s)" % (ids[1:-1],)
    sql = """SELECT station,
         high, low, precip * 25.4 as p01d
         from alldata_la
         WHERE station in %s and 
         day = '%s'""" % ( ids, ts.strftime("%Y-%m-%d") )
    print sql
    rs = coop.query( sql ).dictresult()
    if len(rs) > 4:
        grid_high(nc, ts, rs)
        grid_low(nc, ts, rs)
        grid_p01d(nc, ts, rs)
    else:
        print "%s has %02i entries, FAIL" % (ts.strftime("%Y-%m-%d %H:%M"), 
            len(rs))

#
#create_netcdf()
#sys.exit()
load_stationtable()
nc = netCDF3.Dataset("data/coopgrid.nc", 'a')
now = sts
#now = mx.DateTime.DateTime(1980,1,1)
while now < ets:
  if now.day == 1 and now.month == 1:
    reset_precip(now.year)
  if now not in badtimes:
    grid_hour(nc , now)
  now += mx.DateTime.RelativeDateTime(days=1)

nc.close()
