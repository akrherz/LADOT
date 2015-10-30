from matplotlib import pyplot as plt
import numpy
import netCDF4
from pyiem.datatypes import distance
import numpy as np
import psycopg2
import calendar
pgconn = psycopg2.connect(database='coop', host='iemdb', user='nobody')
cursor = pgconn.cursor()

nc = netCDF4.Dataset('../data/monthly_deltas.nc')
lats = nc.variables['lat'][:]
lons = nc.variables['lon'][:]
deltap = nc.variables['deltap'][:]


def find_delta_yx(lon, lat):
    mindist = 100.0
    (ys, xs) = numpy.shape(lats)
    for y in range(ys):
        for x in range(xs):
            dist = ((lats[y, x] - lat) ** 2 + (lons[y, x] - lon) ** 2)**0.5
            if dist < mindist:
                mindist = dist
                minx = x
                miny = y

    return miny, minx


def get_month(month):
    vals = []
    for line in open('/mesonet/share/pickup/ladot/virtual_station.dat'):
        tokens = line.split(",")
        lat = float(tokens[3])
        lon = float(tokens[4])

        y, x = find_delta_yx(lon, lat)

        vals.append(deltap[month - 1, y, x])

    return np.average(vals)

cursor.execute("""with data as (
  SELECT station, year, month, sum(precip) from alldata_la
  GROUP by station, month, year)

    SELECT month, avg(sum) from data GROUP by month ORDER by month ASC""")
avgs = []
for row in cursor:
    avgs.append(float(row[1]))

deltas = [get_month(x) for x in range(1, 13)]
deltas = np.array(deltas)
avgs = np.array(avgs)
avgsMM = distance(avgs, 'IN').value('MM')
scenMM = avgsMM * (100. + deltas) / 100.

fig, ax = plt.subplots(2, 1)

width = 0.40
bar1 = ax[0].bar(np.arange(12), avgsMM, width, color='b', label='Contemporary')
bar2 = ax[0].bar(np.arange(12) + width, scenMM, width, color='r',
                 label='Future')

ax[0].set_xticks(np.arange(12) + width)
ax[0].set_xticklabels(calendar.month_abbr[1:])
ax[0].set_ylabel("Precipitation [mm]")
ax[0].set_title("Louisiana Average Monthly Precipitation")
ax[0].grid(True)
ax[0].set_ylim(0, 200)
ax[0].legend(ncol=2, fontsize=12)

ax[1].bar(np.arange(12), deltas, color='k', width=0.8)
ax[1].grid(True)
ax[1].set_ylabel("Precipitation Delta [%]")
ax[1].set_xticks(np.arange(12) + width)
ax[1].set_xticklabels(calendar.month_abbr[1:])

fig.savefig("test.png")
