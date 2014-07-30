import netCDF3
from matplotlib import pyplot as plt
import numpy
from scipy import stats

def k2f(ar):
  return (ar - 273.5 ) * 9.0/5.0 + 32.0

nc = netCDF3.Dataset("data/coopgrid.nc", 'r')
high = nc.variables["low"]

fig = plt.figure()
ax = fig.add_subplot(111)
width = 0.35

n, bins, patches = ax.hist(numpy.ravel( k2f(high[:,10,10])), 100, normed=1)

#ax.set_xticklabels( range(1970,2010,5) )
ax.set_xlabel("Temperature [F]")
ax.set_ylabel("Frequency")
ax.set_title("Daily Low Temperature Distribution for Mississippi (1970-2009)")

#ax.set_ylim( 40, 90 )
ax.grid( True )

fig.savefig("test.png")

