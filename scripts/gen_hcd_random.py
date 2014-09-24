"""
 Generate .hcd file for Pavement Design Guide

YYYYMMDDHH,Temperature (F),Wind speed (mph),% Sun shine, Precipitation, 
    Relative humidity.
"""
import math
import numpy as np
import netCDF4
import mx.DateTime
import subprocess

periods = [
  [mx.DateTime.DateTime(1970,1,1), mx.DateTime.DateTime(1974,12,31)],
  [mx.DateTime.DateTime(1975,1,1), mx.DateTime.DateTime(1978,12,31)],
  [mx.DateTime.DateTime(1979,1,1), mx.DateTime.DateTime(1982,12,31)],
  [mx.DateTime.DateTime(1983,1,1), mx.DateTime.DateTime(1987,12,31)],
  [mx.DateTime.DateTime(1988,1,1), mx.DateTime.DateTime(1991,12,31)],
  [mx.DateTime.DateTime(1992,1,1), mx.DateTime.DateTime(1996,12,31)],
  [mx.DateTime.DateTime(1997,1,1), mx.DateTime.DateTime(2002,12,31)],
  [mx.DateTime.DateTime(2003,1,1), mx.DateTime.DateTime(2009,12,31)],
]
# Result of my random sorting [914, 12, 335, 654, 162, 644, 294, 903]
periods2 = [
  [mx.DateTime.DateTime(1975,1,1), mx.DateTime.DateTime(1978,12,31,23)],
  [mx.DateTime.DateTime(1988,1,1), mx.DateTime.DateTime(1991,12,31,23)],
  [mx.DateTime.DateTime(1997,1,1), mx.DateTime.DateTime(2002,12,31,23)],
  [mx.DateTime.DateTime(1979,1,1), mx.DateTime.DateTime(1982,12,31,23)],
  [mx.DateTime.DateTime(1992,1,1), mx.DateTime.DateTime(1996,12,31,23)],
  [mx.DateTime.DateTime(1983,1,1), mx.DateTime.DateTime(1987,12,31,23)],
  [mx.DateTime.DateTime(2003,1,1), mx.DateTime.DateTime(2009,12,31,23)],
  [mx.DateTime.DateTime(1970,1,1), mx.DateTime.DateTime(1974,12,31,23)],
]

onehour = mx.DateTime.RelativeDateTime(hours=1)

anc = netCDF4.Dataset("../data/asosgrid.nc", 'r')
atmpk = anc.variables['tmpk']
asmps = anc.variables['smps']
askyc = anc.variables['skyc']
arelh = anc.variables['relh']
ap01m = anc.variables['p01m']

cnc = netCDF4.Dataset("../data/coopgrid.nc", 'r')
ahigh = cnc.variables['high']
alow = cnc.variables['low']
ap01d = cnc.variables['p01d']

def hourly_fitter_temp(asos, base, trange):
    """
    Use the hourly fit of asos data to do something with the COOP data
    """
    weights = (asos - np.min(asos)) / (np.max(asos) - np.min(asos))
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
        return [0.]*len(asos)
    if np.sum(asos) == 0 and coop > 0:
        asos = [0]*len(asos)
        asos[15:19] = [1.,2.,1.,1.] # Fake Storm signature
    weights = asos / np.sum( asos )
    return coop * weights


def boundschk(val, pval, lower, upper):
    v = np.where( val >= lower, val, pval[:len(val)])
    v = np.where( v <= upper, v, pval[:len(val)])
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

def runner():
    for line in open('../data/station.dat'):
        tokens = line.split(",")
        stid = tokens[0]
        #if stid != '00070':
        #  continue
        lat = float(tokens[3])
        lon = float(tokens[4])
    
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
        fmt = ",%.1f,%.1f,%.1f,%.2f,%.1f\n"
        sts = mx.DateTime.DateTime(1970,1,1,7)
        ets = mx.DateTime.DateTime(2010,1,1,7)
        BASE = mx.DateTime.DateTime(1970,1,1,0)
        END = mx.DateTime.DateTime(2010,1,1,0)
        MAXSZ = int((END - BASE).hours )
        MAXDSZ = int((END - BASE).days )
        interval = mx.DateTime.RelativeDateTime(days=1)
    
        now = sts
        p_tmpf = [0]*24
        p_mph = [0]*24
        p_psun = [0]*24
        p_phour = [0]*24
        p_relh = [0]*24
    
        ds = {}
    
        # We need to bootstrap the first 7 hours
        tmpf = k2f( s_atmpk[:7] ) # Stored in K
        mph = s_asmps[:7] * 2.0
        psun = 100. - s_askyc[:7]
        phour = s_ap01m[:7] / 25.4 # Convert to inches
        relh = s_arelh[:7]
        for i in range(len(tmpf)):
            ts = now - mx.DateTime.RelativeDateTime(hours=(7-i))
            ds[ts] = fmt % ( tmpf[i], mph[i], 
                           psun[i], phour[i], relh[i]) 

 
        while now < ets:
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
            high = k2f( s_high[coffset] )
            low = k2f( s_low[coffset] )
            p01d = s_p01d[coffset] / 25.4 # Convert to inches
            
            # we smear the temperature data
            tmpf = hourly_fitter_temp(tmpf, low, high - low) 
            tmpf, p_tmpf = boundschk(tmpf, p_tmpf, -20., 120.)
            
            # we smear the precipitation data
            phour = hourly_fitter_precip(phour, p01d)
            phour, p_phour = boundschk(phour, p_phour, 0.0, 10.) # 0,10 inches
            #if p01d > 4:
            #  print phour, p01d
            
            # can't touch these
            mph, p_mph = boundschk(mph, p_mph, 0.0, 100.)
            psun, p_psun = boundschk(psun, p_psun, 0.0, 100.)
            relh, p_relh = boundschk(relh, p_relh, 0.0, 100.) # 0,100 %
            
            for i in range(len(tmpf)):
                ts = now + mx.DateTime.RelativeDateTime(hours=i)
                ds[ts] = fmt % (tmpf[i], mph[i], 
                               psun[i], phour[i], relh[i]) 

            now += interval

        # Okay, time for magic shifter...
        realts = mx.DateTime.DateTime(1970,1,1,0)
        for (sts,ets) in periods:
            now = sts
            while now <= ets:
                out.write( realts.strftime("%Y%m%d%H") + ds[now] )
                now += onehour
                realts += onehour

        out.close()
        subprocess.call("python qc_hcd.py %s.hcd" % (stid,), shell=True)

if __name__ == '__main__':
    runner()
