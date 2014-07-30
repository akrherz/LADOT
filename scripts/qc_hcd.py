# Need something to account for my incompetence 

import sys
import numpy
import mx.DateTime

tmpf = numpy.zeros( (350640), numpy.float )
relh = numpy.zeros( (350640), numpy.float )
p01i = numpy.zeros( (350640), numpy.float )
smps = numpy.zeros( (350640), numpy.float )
psun = numpy.zeros( (350640), numpy.float )


cnt = 0
for line in open(sys.argv[1]):
  # ts, tmpf[i], mph[i], psun[i], phour[i], relh[i]
  tokens = line.split(",")
  tmpf[cnt] = float(tokens[1])
  smps[cnt] = float(tokens[2])
  psun[cnt] = float(tokens[3])
  p01i[cnt] = float(tokens[4])
  relh[cnt] = float(tokens[5])
  if float(tokens[5]) < 5:
    print cnt, float(tokens[5])
  cnt += 1

if cnt != 350640:
  print "ERROR, invalid filelength: %s" % (cnt,)

print "Air Temp  [ F ] MIN: %6.2f MAX: %6.2f   AVG: %6.2f" % (numpy.min(tmpf), 
  numpy.max(tmpf), numpy.average(tmpf))
print "Wnd Speed [mph] MIN: %6.2f MAX: %6.2f   AVG: %6.2f" % (numpy.min(smps), 
  numpy.max(smps), numpy.average(smps))
print "Perct Sun [ %% ] MIN: %6.2f MAX: %6.2f   AVG: %6.2f" % (numpy.min(psun), 
  numpy.max(psun), numpy.average(psun))
print "Hr Precip [in ] MIN: %6.2f MAX: %6.2f YRAVG: %6.2f" % (numpy.min(p01i), 
  numpy.max(p01i), numpy.average(p01i)*365.0*24.0)
print "Rel Humid [ %% ] MIN: %6.2f MAX: %6.2f   AVG: %6.2f" % (numpy.min(relh),
  numpy.max(relh), numpy.average(relh))

"""
for yr in range(40):
  print yr+1970,
  for mo in range(12):
    i0 = int((mx.DateTime.DateTime(1970+yr,mo+1,1) - mx.DateTime.DateTime(1970,1,1)).days) * 24
    i1 = int((mx.DateTime.DateTime(1970+yr,mo+1,1) + mx.DateTime.RelativeDateTime(months=1) - mx.DateTime.DateTime(1970,1,1)).days) * 24
    print "%4.1f" % ( numpy.sum( p01i[i0:i1] ), ),
  print
"""
