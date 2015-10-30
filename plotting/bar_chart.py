import netCDF4
import matplotlib.pyplot as plt
import numpy
import psycopg2

def find_yx(nc, lon, lat):
    lats = nc.variables['lat'][:]
    lons = nc.variables['lon'][:]
    nc.close()

    mindist = 100.0
    (ys,xs) = numpy.shape(lats)
    for y in range(ys):
        for x in range(xs):
            dist = ((lats[y,x] - lat) ** 2 + (lons[y,x] - lon) ** 2)**0.5
            if dist < mindist:
                mindist = dist
                minx = x
                miny = y

    return miny, minx


def doit(lon, lat, name):
  nc = netCDF4.Dataset('../data/monthly_deltas.nc')
  (y,x) = find_yx(nc, lon, lat) 

  T = nc.variables['deltat'][:,y,x]
  print T
  P = nc.variables['deltap'][:,y,x]
  nc.close()

  (fig, ax) = plt.subplots(2,1, sharex=True)

  ax[0].set_title("HRM3-HADCM3 Monthly Deltas for Louisiana Climate Zone %s" % (name,))
  ax[0].bar(numpy.arange(1,13)-0.4, T, fc='r')
  ax[0].grid(True)
  ax[0].set_ylabel("Temperature Change $^{\circ}\mathrm{F}$")
  ax[0].set_xlim(0,13)
  ax[0].set_ylim(0, max(T) + 1)
  for i in range(12):
    ax[0].text(i+1, T[i] + 0.1, "%.1f" % (T[i],), ha='center')


  bars = ax[1].bar(numpy.arange(1,13)-0.4, P, fc='b')
  for bar in bars:
    if bar.get_y() < 0:
        bar.set_facecolor('r')
  for i in range(12):
    va = 'bottom'
    delta = 2
    if P[i] < 0:
        va = 'top'
        delta = -2
    ax[1].text(i+1, P[i] + delta, "%.0f" % (P[i],), ha='center', va=va)
  ax[1].set_xticklabels( ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec') )
  ax[1].set_xticks( numpy.arange(1,13) )
  ax[1].grid(True)
  ax[1].set_ylabel("Precip Change %")
  ax[1].set_ylim( min(P) - 10, max(P) + 10)

  fig.savefig('mon_%s_chart.png' % (name,))
  del fig

POSTGIS = psycopg2.connect(database='postgis', host='iemdb', user='nobody')
pcursor = POSTGIS.cursor()

pcursor.execute("""select ST_x(ST_Centroid(the_geom)), ST_y(ST_Centroid(the_geom)),
name from climate_div WHERE st = 'LA'""")

for row in pcursor:
    doit(row[0], row[1], row[2])
