# Grid hourly ASOS data please
# Temperature (K),
# Wind speed (mps),
# % Sun shine, 
# Precipitation, 
# Relative humidity.

import sys
import netCDF3
import numpy
import mx.DateTime
from pyIEM import iemdb, mesonet
import Ngl
i = iemdb.iemdb()
mesosite = i['mesosite']
asos = i['asos']
locs = {}

badtimes = [
 ]
BASE = mx.DateTime.DateTime(1970,1,1)
STOP = mx.DateTime.DateTime(2012,1,1)
sts = mx.DateTime.DateTime(2010,1,1)
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
    #sql = """SELECT id, x(geom) as lon, y(geom) as lat from
    #     stations where network IN ('MS_ASOS','LA_ASOS','AR_ASOS',
    #     'TN_ASOS','AL_ASOS', 'FL_ASOS') and id in ('EGI', 'ELD', 'ESF', 
    #     'HOT', 'HSV', 'LCH', 'LFT', 'MEM', 'MGM', 'MLU',
    #     'MOB', 'MSY', 'NBG', 'NEW', 'NMM', 'NPA', 'NSE', 'PBF', 'PNS', 'POE',
    #     'SHV', 'TCL', 'TXK')"""
    sql = """SELECT id, x(geom) as lon, y(geom) as lat from
         stations where network IN ('MS_ASOS','LA_ASOS','AR_ASOS',
         'TX_ASOS','OK_ASOS')"""
    rs = mesosite.query( sql ).dictresult()
    for i in range(len(rs)):
        locs[ rs[i]['id'] ] = rs[i]
    ids = `locs.keys()`
    ids = "(%s)" % (ids[1:-1],)

def create_netcdf():
    nc = netCDF3.Dataset("data/asosgrid.nc", 'w')
    # Dimensions
    nc.createDimension('time', int((STOP - BASE).hours) )
    nc.createDimension('lat', int((NORTH-SOUTH)/DELTAY) )
    nc.createDimension('lon', int((EAST-WEST)/DELTAX)   )
    # Variables
    tm = nc.createVariable('time', 'd', ('time',))
    tm.units = 'hours since 1970-01-01'
    tm[:] = range( int((STOP - BASE).hours) )

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

    tmpk = nc.createVariable('tmpk', 'f', ('time','lat','lon'))
    tmpk.units = 'K'
    tmpk._FillValue = 1.e20
    tmpk.missing_value = 1.e20
    tmpk.long_name = 'Surface Air Temperature'

    smps = nc.createVariable('smps', 'f', ('time','lat','lon'))
    smps.units = 'm s-1'
    smps._FillValue = 1.e20
    smps.missing_value = 1.e20
    smps.long_name = '10m Wind Speed'

    skyc = nc.createVariable('skyc', 'f', ('time','lat','lon'))
    skyc.units = '%'
    skyc._FillValue = 1.e20
    skyc.missing_value = 1.e20
    skyc.long_name = 'Sky Coverage'

    p01m = nc.createVariable('p01m', 'f', ('time','lat','lon'))
    p01m.units = 'mm'
    p01m._FillValue = 1.e20
    p01m.missing_value = 1.e20
    p01m.long_name = 'Precipitation'

    relh = nc.createVariable('relh', 'f', ('time','lat','lon'))
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
            lats.append(  locs[rs[i]['station']]['lat'] )
            lons.append(  locs[rs[i]['station']]['lon'] )
            vals.append( rs[i]['max_sknt'] * 0.514 ) # knots to mps
    if len(vals) < 4:
        print "No WIND data at all for time: %s" % (ts,)   
        return
    grid = Ngl.natgrid(lons, lats, vals, XAXIS, YAXIS)
    grid = Ngl.natgrid(lons, lats, vals, XAXIS, YAXIS)
    offset = int((ts - BASE).hours )
    if grid is not None:
        gt = grid.transpose()
        nc.variables['smps'][offset,:,:] = numpy.where(gt < 0., 0., gt)
    else:
        print "WIND gridding failed, len vals %s" % (len(vals),)


def grid_relh(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        if rs[i]['max_tmpf'] is not None and rs[i]['max_dwpf'] is not None:
            lats.append(  locs[rs[i]['station']]['lat'] )
            lons.append(  locs[rs[i]['station']]['lon'] )
            vals.append( mesonet.relh( rs[i]['max_tmpf'], rs[i]['max_dwpf'] ) )
    if len(vals) < 4:
        print "No RELH data at all for time: %s" % (ts,)   
        return
    grid = Ngl.natgrid(lons, lats, vals, XAXIS, YAXIS)
    offset = int((ts - BASE).hours )
    if grid is not None:
        gt = grid.transpose()
        gt = numpy.where(gt < 100.1, gt, 100.)
        nc.variables['relh'][offset,:,:] = numpy.where(gt < 12., 12., gt)
    else:
        print "RELH gridding failed, len vals %s" % (len(vals),)

def grid_skyc(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        v =  max(rs[i]['max_skyc1'], rs[i]['max_skyc2'], rs[i]['max_skyc3'])
        if v is not None:
            lats.append(  locs[rs[i]['station']]['lat'] )
            lons.append(  locs[rs[i]['station']]['lon'] )
            vals.append( float(v) )
    if len(vals) < 4:
        print "No SKYC data at all for time: %s" % (ts,)   
        return
    grid = Ngl.natgrid(lons, lats, vals, XAXIS, YAXIS)
    offset = int((ts - BASE).hours )
    if grid is not None:
        gt = grid.transpose()
        gt = numpy.where(gt > 0., gt, 0.0)
        nc.variables['skyc'][offset,:,:] = numpy.where(gt > 100., 100., gt)
    else:
        print "SKYC gridding failed, len vals %s" % (len(vals),)
        print vals

def grid_tmpf(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        if rs[i]['max_tmpf'] is not None:
            lats.append(  locs[rs[i]['station']]['lat'] )
            lons.append(  locs[rs[i]['station']]['lon'] )
            vals.append( mesonet.f2k( rs[i]['max_tmpf'] ) )
    if len(vals) < 4:
        print "No TMPF data at all for time: %s" % (ts,)   
        return
    grid = Ngl.natgrid(lons, lats, vals, XAXIS, YAXIS)
    offset = int((ts - BASE).hours ) 
    if grid is not None:
        nc.variables['tmpk'][offset,:,:] = grid.transpose()
    else:
        print "TMPK gridding failed, len vals %s" % (len(vals),)

def grid_p01m(nc, ts, rs):
    lats = []
    lons = []
    vals = []
    for i in range(len(rs)):
        lats.append(  locs[rs[i]['station']]['lat'] )
        lons.append(  locs[rs[i]['station']]['lon'] )
        vals.append( rs[i]['max_p01m'] )
    grid = Ngl.natgrid(lons, lats, vals, XAXIS, YAXIS)
    offset = int((ts - BASE).hours )
    if grid is not None:
        gt = grid.transpose()
        nc.variables['p01m'][offset,:,:] = numpy.where(gt > 0., gt, 0.)
    else:
        print "P01M gridding failed, len vals %s" % (len(vals),)

def grid_hour(nc, ts):
    ids = `locs.keys()`
    ids = "(%s)" % (ids[1:-1],)
    sql = """SELECT station,
         max(case when tmpf > -60 and tmpf < 130 THEN tmpf else null end) as max_tmpf,
         max(case when sknt > 0 and sknt < 100 then sknt else 0 end) as max_sknt,
         max(getskyc(skyc1)) as max_skyc1,
         max(getskyc(skyc2)) as max_skyc2,
         max(getskyc(skyc3)) as max_skyc3,
         max(case when p01i > 0 and p01i < 1000 then p01i * 25.4 else 0 end) as max_p01m,
         max(case when dwpf > -60 and dwpf < 100 THEN dwpf else null end) as max_dwpf from t%s  
         WHERE station in %s and 
         valid >= '%s' and valid < '%s' GROUP by station""" % (
         ts.gmtime().year, ids, 
         ts.strftime("%Y-%m-%d %H:%M"),
     (ts + mx.DateTime.RelativeDateTime(hours=1)).strftime("%Y-%m-%d %H:%M") )
    rs = asos.query( sql ).dictresult()
    if len(rs) > 4:
        grid_tmpf(nc, ts, rs)
        grid_relh(nc, ts, rs)
        grid_wind(nc, ts, rs)
        grid_skyc(nc, ts, rs)
        grid_p01m(nc, ts, rs)
    else:
        print "%s has %02i entries, FAIL" % (ts.strftime("%Y-%m-%d %H:%M"), 
            len(rs))

#
#create_netcdf()
#sys.exit()
load_stationtable()
nc = netCDF3.Dataset("data/asosgrid.nc", 'a')
now = sts
#now = mx.DateTime.DateTime(1980,1,1)
while now < ets:
  #print now
  #if now not in badtimes:
  grid_hour(nc , now)
  now += mx.DateTime.RelativeDateTime(hours=1)

nc.close()
