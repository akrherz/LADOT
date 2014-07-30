# Generate .hcd file for Pavement Design Guide

#YYYYMMDDHH,Temperature (F),Wind speed (mph),% Sun shine, Precipitation, Relative humidity.

import os
import pdb
import math
import numpy
import numpy.ma
import netCDF4
import iemdb
POSTGIS = iemdb.connect('postgis')
pcursor = POSTGIS.cursor()
import mx.DateTime, sys
from pyIEM import mesonet

anc = netCDF4.Dataset("data/asosgrid.nc", 'r')
atmpk = anc.variables['tmpk']
asmps = anc.variables['smps']
askyc = anc.variables['skyc']
arelh = anc.variables['relh']
ap01m = anc.variables['p01m']

cnc = netCDF4.Dataset("data/coopgrid.nc", 'r')
ahigh = cnc.variables['high']
alow = cnc.variables['low']
ap01d = cnc.variables['p01d']

def hourly_fitter_temp(asos, base, trange):
    """
    Use the hourly fit of asos data to do something with the COOP data
    """
    weights = (asos - numpy.min(asos)) / (numpy.max(asos) - numpy.min(asos))
    #if (base + trange) > 100:
    #   print
    #   print trange, base
    #   print weights
    #   print base + ( trange * weights )
      
    return base + ( trange * weights )


def hourly_fitter_precip(asos, coop):
    """
    Use the hourly fit of asos data to do something with the COOP data
    """
    if coop  == 0:
        return numpy.ma.array([0.]*len(asos))
    if numpy.sum(asos) == 0 and coop > 0:
        asos = numpy.ma.zeros( (len(asos)) )
        asos[1:5] = [1,2,1,1] # Fake storm
    weights = asos / numpy.sum( asos )
    return coop * weights


def boundschk(vname, val, pval, lower, upper):
  #print vname, val 
  v = numpy.ma.where( val.filled() >= lower, val.filled(), pval[:len(val)])
  v = numpy.ma.where( v <= upper, v, pval[:len(val)])
  return v, v

def k2f(thisk):
  return (9.00/5.00 * (thisk - 273.15) ) + 32.00


def computeIJ(lon, lat):
  lats = anc.variables['lat'][:]
  lons = anc.variables['lon'][:]
  mindist = 100000
  for j in range(len(lats)):
    for i in range(len(lons)):
      dist = math.sqrt( (lons[i] - lon)**2 + (lats[j] - lat)**2 )
      if dist < mindist:
        mindist = dist
        mini = i
        minj = j

  return mini, minj

DELTAT = [1.3, 1.5, 3.2, 2.4, 2.7, 1.9, 2.8, 2., 1.7, 1.9, 2.0, 1.4]
DELTAP = [-0.208, -0.133,  0.407,  0.144,  0.388, 0.132, -0.174,  
           0.173,  0.682, -0.197,  0.374, -0.114]

def find_delta_yx(lon, lat):
    nc = netCDF4.Dataset('monthly_deltas.nc')
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

def get_deltas(lon, lat):
    """ 
    Compute what our deltas are for this location, moved to a climate
    sector first!
    """
    pcursor.execute("""SELECT ST_x(ST_Centroid(the_geom)), ST_y(ST_Centroid(the_geom))
    from climate_div WHERE st = 'LA' and ST_Contains(the_geom, ST_GeomFromText('SRID=4326;POINT(%s %s)'))""" % (lon, lat))
    row = pcursor.fetchone()
    lon = row[0]
    lat = row[1]
    y,x = find_delta_yx(lon, lat)
    nc = netCDF4.Dataset('monthly_deltas.nc')

    T = nc.variables['deltat'][:,y,x]
    P = nc.variables['deltap'][:,y,x] / 25.4
    nc.close()
    return T, P

def runner():
  for line in open('/mesonet/share/pickup/ladot/virtual_station2.dat'):
    tokens = line.split(",")
    stid = tokens[0]
    #if stid != '00070':
    #  continue
    lat = float(tokens[3])
    lon = float(tokens[4])

    DELTAT, DELTAP = get_deltas( lon, lat )
    gridx, gridy = computeIJ( lon, lat )

    print stid, tokens[1], gridx, gridy

    s_atmpk = atmpk[:,gridy,gridx]
    s_asmps = asmps[:,gridy,gridx]
    s_askyc = askyc[:,gridy,gridx]
    s_ap01m = ap01m[:,gridy,gridx]
    s_arelh = arelh[:,gridy,gridx]

    s_high = ahigh[:,gridy,gridx]
    s_low = alow[:,gridy,gridx]
    s_p01d = ap01d[:,gridy,gridx]

    out = open("%s.hcd" % (stid,), 'w')
    fmt = "%s,%.1f,%.1f,%.1f,%.2f,%.1f\n"
    sts = mx.DateTime.DateTime(2010,1,1,7)
    ets = mx.DateTime.DateTime(2050,1,1,7)
    BASE = mx.DateTime.DateTime(2010,1,1,0)
    END = mx.DateTime.DateTime(2050,1,1,0)
    MAXSZ = int((END - BASE).hours )
    MAXDSZ = int((END - BASE).days )
    interval = mx.DateTime.RelativeDateTime(days=1)

    now = sts
    p_tmpf = [0]*24
    p_mph = [0]*24
    p_psun = [0]*24
    p_phour = [0]*24
    p_relh = [0]*24

    # We need to bootstrap the first 7 hours
    tmpf = k2f( s_atmpk[:7] ) # Stored in K
    mph = s_asmps[:7] * 2.0
    psun = 100. - s_askyc[:7]
    phour = s_ap01m[:7] / 25.4 # Convert to inches
    relh = s_arelh[:7]
    for i in range(len(tmpf)):
      ts = now - mx.DateTime.RelativeDateTime(hours=(7-i))
      out.write(fmt % (ts.strftime("%Y%m%d%H"), tmpf[i], mph[i], 
                     psun[i], phour[i], relh[i]) )

 
    while now < ets:
      high_off = DELTAT[ now.month - 1]
      low_off = DELTAT[ now.month - 1]

      aoffset1 = int((now - BASE).hours )
      aoffset2 = int(((now + mx.DateTime.RelativeDateTime(days=1)) - BASE).hours )
      if aoffset1 < 0:
        aoffset1 = 0
      if aoffset2 >= MAXSZ:
        aoffset2 = MAXSZ
      coffset = int((now - BASE).days ) + 1
      if coffset >= MAXDSZ:
        coffset = MAXDSZ - 1

      tmpf = k2f( s_atmpk[aoffset1:aoffset2] ) # Stored in K
      mph = s_asmps[aoffset1:aoffset2] * 2.0
      psun = 100. - s_askyc[aoffset1:aoffset2]
      phour = s_ap01m[aoffset1:aoffset2] / 25.4 # Convert to inches
      relh = s_arelh[aoffset1:aoffset2]
      high = k2f( s_high[coffset] ) + high_off
      low = k2f( s_low[coffset] ) + low_off
      p01d = s_p01d[coffset] / 25.4 # Convert to inches
      # Apply the delta
      if p01d > 0.01:
        p01d += DELTAP[ now.month - 1 ] * p01d

      # we smear the temperature data
      tmpf = hourly_fitter_temp(tmpf, low, high - low) 
      tmpf, p_tmpf = boundschk('t',tmpf, p_tmpf, -20., 120.)

      # we smear the precipitation data
      phour = hourly_fitter_precip(phour, p01d)
      phour, p_phour = boundschk('p',phour, p_phour, 0.0, 10.) # 0,10 inches
      #if p01d > 4:
      #  print phour, p01d

      # can't touch these
      mph, p_mph = boundschk('m',mph, p_mph, 0.0, 100.)
      psun, p_psun = boundschk('s',psun, p_psun, 0.0, 100.)
      relh, p_relh = boundschk('r',relh, p_relh, 0.0, 100.) # 0,100 %

      for i in range(len(tmpf)):
        ts = now + mx.DateTime.RelativeDateTime(hours=i)
        s = fmt % (ts.strftime("%Y%m%d%H"), tmpf[i], mph[i], 
                       psun[i], phour[i], relh[i]) 
        if s.find('nan') > -1:
            pdb.set_trace()
        out.write(s)
      now += interval
    out.close()
    os.system("python qc_hcd.py %s.hcd" % (stid,))

runner()
#print boundschk([0,12,3,4], [10,-4,13,14], 0, 10)
