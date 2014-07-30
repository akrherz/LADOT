import netCDF4
import numpy

nc = netCDF4.Dataset('monthly_deltas.nc', 'w')
nc.createDimension('month', 12)
nc.createDimension('x', 155)
nc.createDimension('y', 130)

lat = nc.createVariable('lat', 'd', ('y','x'))
lon = nc.createVariable('lon', 'd', ('y','x'))
deltat = nc.createVariable('deltat', 'f', ('month', 'y','x'))
deltap = nc.createVariable('deltap', 'f', ('month', 'y','x'))

nc.close()
