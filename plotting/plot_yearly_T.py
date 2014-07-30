import netCDF4 as netCDF3
from matplotlib import pyplot as plt
import numpy
from scipy import stats

def k2f(ar):
  return (ar - 273.5 ) * 9.0/5.0 + 32.0

nc = netCDF3.Dataset("data/yr_coopgrid.nc", 'r')
high = nc.variables["high"]
low = nc.variables["low"]

hdata = numpy.zeros((40,), 'f')
ldata = numpy.zeros((40,), 'f')
for i in range(40):
  hdata[i] = numpy.average(high[i,:,:])
  ldata[i] = numpy.average(low[i,:,:])

fig = plt.figure()
ax = fig.add_subplot(111)
width = 0.35
bar1 = ax.bar( numpy.arange(40), k2f(hdata), width, color='r' )
bar2 = ax.bar( numpy.arange(40) + width, k2f(ldata), width, color='b' )

ax.set_xticklabels( range(1970,2010,5) )
ax.set_ylabel("Temperature [F]")
ax.set_xlabel("Year")
ax.set_title("Mean Yearly Temperature for Louisiana (1970-2009)")

h_slope, intercept, h_r_value, p_value, std_err = stats.linregress(numpy.arange(40),k2f(hdata))
ax.plot([0,40], [intercept, (40 * h_slope) + intercept], color='r')
print h_r_value, p_value

l_slope, intercept, l_r_value, p_value, std_err = stats.linregress(numpy.arange(40),k2f(ldata))
ax.plot([0,40], [intercept, (40 * l_slope) + intercept], color='b')
print l_r_value, p_value

ax.legend( (bar1[0],bar2[0]), (r"$High: \frac{dT}{dyear} = %.3f , R^2 = %.2f$" % (h_slope, h_r_value ** 2), r"$Low: \frac{dT}{dyear} = %.3f , R^2 = %.2f$" % (l_slope, l_r_value ** 2)) )

ax.set_ylim( 40, 90 )

fig.savefig("test.png")

