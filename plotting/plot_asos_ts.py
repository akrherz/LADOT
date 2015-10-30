import netCDF4
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import datetime
from pyiem.datatypes import temperature

nc = netCDF4.Dataset('../data/asosgrid.nc', 'r')
tmpk = nc.variables['tmpk']
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
BASE = datetime.datetime(1970, 1, 1)

sts = datetime.datetime(1975, 12, 30, 7)
offset = int((sts - BASE).total_seconds() / 3600)

data = temperature(tmpk[offset:(offset+72), 8, 7], 'K').value('F')

(fig, ax) = plt.subplots(1, 1)

ax.bar(np.arange(72), data, facecolor='b', edgecolor='b')
ax.set_ylim(min(data)-3, max(data)+3)
ax.grid(True)
fig.savefig('test.png')
plt.close()

nc.close()
 