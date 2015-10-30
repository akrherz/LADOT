import netCDF4
from pyiem.plot import MapPlot
import calendar
import matplotlib.cm as cm

nc = netCDF4.Dataset('../data/monthly_deltas.nc', 'r')
deltat = nc.variables['deltat']
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
for i in range(12):

    data = deltat[i, :, :]

    m = MapPlot(sector='state', state='LA', nologo=True,
                caption='LADOT Project, map by daryl herzmann',
                title="Temperature Deltas $^\circ$C for %s" % (
                    calendar.month_name[i+1], ))
    p = m.pcolormesh(lon, lat, data, range(-10, 11, 2),
                     cmap=cm.get_cmap('RdBu_r'), units='C')
    m.postprocess(filename=('monthly_tdeltas_%s.png'
                            ) % (calendar.month_name[i+1], ))
    m.close()

nc.close()
