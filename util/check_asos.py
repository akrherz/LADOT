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

nc = netCDF3.Dataset("data/asosgrid.nc")
now = sts
while now < ets:
    offset = int((now - BASE).hours )
    print now.year, numpy.average( nc.variables['smps'][offset,:,:] )
    now += mx.DateTime.RelativeDateTime(years=1)

nc.close()
