import netCDF3
from matplotlib import pyplot as plt
import numpy
from scipy import stats


nc = netCDF3.Dataset("data/yr_coopgrid.nc", 'r')
precip = nc.variables["p01d"]

pdata = numpy.zeros((40,), 'f')
for i in range(40):
  pdata[i] = numpy.average(precip[i,:,:]) * 365.0 / 25.4 # in/yr

fig = plt.figure()
ax = fig.add_subplot(111)
width = 0.70
bar1 = ax.bar( numpy.arange(40), pdata, width, color='b' )

ax.set_xticklabels( range(1970,2010,5) )
ax.set_ylabel("Precipitation [inch]")
ax.set_xlabel("Year")
ax.set_title("Mean Yearly Precipitation for Mississippi (1970-2009)")

slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.arange(40),pdata)
ax.plot([0,40], [intercept, (40 * slope) + intercept], color='b')
print r_value, p_value

ax.legend( (bar1[0],), (r"$Precip: \frac{dP}{dyear} = %.3f , R^2 = %.2f$" % (slope, r_value ** 2),))

ax.set_ylim( 40, 90 )

fig.savefig("test.png")

