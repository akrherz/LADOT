# Need something to generate the station table for MPEG

from pyIEM import iemdb
import urllib2, re
i = iemdb.iemdb()
postgis = i['postgis']

rs = postgis.query("SELECT x(centroid(geom)) as lon, y(centroid(geom)) as lat, name from nws_ugc WHERE ugc ~* 'LAC' ORDER by name ASC").dictresult()

for i in range(len(rs)):
  lat = rs[i]['lat']
  lon = rs[i]['lon']
  r = urllib2.urlopen('http://www.earthtools.org/height/%s/%s' % (lat,lon)).read()
  newelev =  float(re.findall("<meters>(.*)</meters>", r)[0]) * 3.28
  # Weather station number, weather station abbreviation, location (city|state), latitude, longitude, elevation, first date in file (YYYMMDD)
  print "%05i,%s| LA,%s,%.3f,%.3f,%.1f,19700101" % (i, rs[i]['name'].upper(), rs[i]['name'], rs[i]['lat'], rs[i]['lon'], newelev)
