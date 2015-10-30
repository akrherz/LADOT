import netCDF4
from pyiem.plot import MapPlot
import datetime

nc = netCDF4.Dataset('../data/asosgrid.nc', 'r')
tmpk = nc.variables['tmpk']
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
BASE = datetime.datetime(1970, 1, 1)

sts = datetime.datetime(1971, 12, 30, 12)
for i in range(72):
    ts = sts + datetime.timedelta(hours=i)
    offset = int((ts - BASE).total_seconds() / 3600)

    data = tmpk[offset, :, :]

    m = MapPlot(sector='state', state='LA', title='%s' % (ts, ))
    p = m.pcolormesh(lon, lat, data, range(270, 310, 2))
    m.postprocess(filename='test%03i.png' % (i,))
    m.close()

nc.close()
