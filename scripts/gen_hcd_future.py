import subprocess
import pdb
import math
import numpy as np
import numpy.ma
import netCDF4
import psycopg2
import datetime
import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

DEBUG = (len(sys.argv) == 2)

POSTGIS = psycopg2.connect(database='postgis', host='iemdb', user='nobody')
pcursor = POSTGIS.cursor()

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

# Temperature periods
A = [datetime.datetime(1970, 1, 1), datetime.datetime(1974, 12, 31, 23)]
B = [datetime.datetime(1975, 1, 1), datetime.datetime(1978, 12, 31, 23)]
C = [datetime.datetime(1979, 1, 1), datetime.datetime(1982, 12, 31, 23)]
D = [datetime.datetime(1983, 1, 1), datetime.datetime(1987, 12, 31, 23)]
E = [datetime.datetime(1988, 1, 1), datetime.datetime(1991, 12, 31, 23)]
F = [datetime.datetime(1992, 1, 1), datetime.datetime(1996, 12, 31, 23)]
G = [datetime.datetime(1997, 1, 1), datetime.datetime(2002, 12, 31, 23)]
H = [datetime.datetime(2003, 1, 1), datetime.datetime(2009, 12, 31, 23)]

# Temperature Periods, Random Sorted
PERIODS = [B, E, G, C, F, D, H, A]
# Temperature Periods, Worst Case
# PERIODS = [G, H, C, D, F, B, E, A]
# Temperature Periods, Best Case (hot last)
# PERIODS.reverse()
"""
# Precipitation Periods
A = [datetime.datetime(1970, 1, 1), datetime.datetime(1975, 12, 31, 23)]
B = [datetime.datetime(1976, 1, 1), datetime.datetime(1979, 12, 31, 23)]
C = [datetime.datetime(1980, 1, 1), datetime.datetime(1985, 12, 31, 23)]
D = [datetime.datetime(1986, 1, 1), datetime.datetime(1989, 12, 31, 23)]
E = [datetime.datetime(1990, 1, 1), datetime.datetime(1992, 12, 31, 23)]
F = [datetime.datetime(1993, 1, 1), datetime.datetime(1999, 12, 31, 23)]
G = [datetime.datetime(2000, 1, 1), datetime.datetime(2004, 12, 31, 23)]
H = [datetime.datetime(2005, 1, 1), datetime.datetime(2009, 12, 31, 23)]

# There is no random sorted percipitation period!
# Precipitation Periods, Worst Case
PERIODS = [E, G, A, C, D, B, F, H]
# Precipitation Periods, Best Case
PERIODS.reverse()
"""
# OFFSET these into the future period 2010->2050
for i in range(len(PERIODS)):
    for j in range(2):
        newyear = PERIODS[i][j].year + 40
        PERIODS[i][j] = PERIODS[i][j].replace(year=newyear)


def hourly_fitter_temp(asos, base, trange):
    """
    Use the hourly fit of asos data to do something with the COOP data
    """
    weights = (asos - numpy.min(asos)) / (numpy.max(asos) - numpy.min(asos))
    # if (base + trange) > 100:
    #   print
    #   print trange, base
    #   print weights
    #   print base + ( trange * weights )

    return base + (trange * weights)


def hourly_fitter_precip(asos, coop):
    """
    Use the hourly fit of asos data to do something with the COOP data
    """
    if coop == 0:
        return numpy.ma.array([0.]*len(asos))
    if numpy.sum(asos) == 0 and coop > 0:
        asos = numpy.ma.zeros((len(asos)))
        asos[1:5] = [1, 2, 1, 1]  # Fake storm
    weights = asos / numpy.sum(asos)
    return coop * weights


def boundschk(vname, val, pval, lower, upper):
    # print vname, val
    v = numpy.ma.where(val >= lower, val, pval[:len(val)])
    v = numpy.ma.where(v <= upper, v, pval[:len(val)])
    return v, v


def k2f(thisk):
    return (9.00/5.00 * (thisk - 273.15)) + 32.00


def computeIJ(lon, lat):
    lats = anc.variables['lat'][:]
    lons = anc.variables['lon'][:]
    mindist = 100000
    for j in range(len(lats)):
        for i in range(len(lons)):
            dist = math.sqrt((lons[i] - lon)**2 + (lats[j] - lat)**2)
            if dist < mindist:
                mindist = dist
                mini = i
                minj = j

    return mini, minj


def find_delta_yx(lon, lat):
    nc = netCDF4.Dataset('../data/monthly_deltas.nc')
    lats = nc.variables['lat'][:]
    lons = nc.variables['lon'][:]
    nc.close()

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


def get_deltas(lon, lat):
    """
    Compute what our deltas are for this location, moved to a climate
    sector first!
    """
    pcursor.execute("""
        SELECT ST_x(ST_Centroid(the_geom)), ST_y(ST_Centroid(the_geom)),
        ST_distance(the_geom, ST_SetSRID(
        ST_GeomFromText('SRID=4326;POINT(%s %s)'), 4326))
        from climate_div WHERE st = 'LA' ORDER by ST_DISTANCE ASC
        """ % (lon, lat))
    row = pcursor.fetchone()
    if row is None:
        return None, None
    lon = row[0]
    lat = row[1]
    y, x = find_delta_yx(lon, lat)
    nc = netCDF4.Dataset('../data/monthly_deltas.nc')

    T = nc.variables['deltat'][:, y, x]
    P = nc.variables['deltap'][:, y, x] / 25.4
    nc.close()
    return T, P


def runner():
    for line in open('/mesonet/share/pickup/ladot/virtual_station.dat'):
        tokens = line.split(",")
        stid = tokens[0]
        # if stid not in ['00035', '00037']:
        #     continue
        lat = float(tokens[3])
        lon = float(tokens[4])

        DELTAT, DELTAP = get_deltas(lon, lat)
        if DELTAT is None:
            print("get_deltas failed for stid: %s" % (stid,))
            continue
        gridx, gridy = computeIJ(lon, lat)

        print("%s %s X:%s Y:%s" % (stid, tokens[1], gridx, gridy,))

        s_atmpk = atmpk[:, gridy, gridx]
        s_asmps = asmps[:, gridy, gridx]
        s_askyc = askyc[:, gridy, gridx]
        s_ap01m = ap01m[:, gridy, gridx]
        s_arelh = arelh[:, gridy, gridx]

        s_high = ahigh[:, gridy, gridx]
        s_low = alow[:, gridy, gridx]
        s_p01d = ap01d[:, gridy, gridx]

        fmt = ",%.1f,%.1f,%.1f,%.2f,%.1f\n"
        # So the timestamps stored in the ASOS grid netCDF file are in local
        # time (sigh), so the MEPDG file has timestamps in UTC, so we
        # bootstrap the first seven hours
        sts = datetime.datetime(2010, 1, 1, 7)
        ets = datetime.datetime(2050, 1, 1, 7)
        BASE = datetime.datetime(2010, 1, 1, 0)
        END = datetime.datetime(2050, 1, 1, 0)
        MAXSZ = int((END - BASE).total_seconds() / 3600.)
        MAXDSZ = int((END - BASE).total_seconds() / 86400.)
        interval = datetime.timedelta(days=1)

        p_tmpf = [0]*24
        p_mph = [0]*24
        p_psun = [0]*24
        p_phour = [0]*24
        p_relh = [0]*24

        ds = {}

        # We need to bootstrap the first 7 hours
        tmpf = k2f(s_atmpk[:7])  # Stored in K
        mph = s_asmps[:7] * 2.0
        psun = 100. - s_askyc[:7]
        phour = s_ap01m[:7] / 25.4  # Convert to inches
        relh = s_arelh[:7]
        for i in range(len(tmpf)):
            ts = sts - datetime.timedelta(hours=(7-i))
            ds[ts] = [tmpf[i], mph[i], psun[i], phour[i], relh[i]]

        now = sts
        while now < ets:
            high_off = DELTAT[now.month - 1]
            low_off = DELTAT[now.month - 1]

            aoffset1 = int((now - BASE).total_seconds() / 3600.) - 6
            aoffset2 = int(((now + datetime.timedelta(days=1)) -
                            BASE).total_seconds() / 3600.) - 6
            if aoffset2 >= MAXSZ:
                aoffset2 = MAXSZ
            coffset = int((now - BASE).days) + 1
            if coffset >= MAXDSZ:
                coffset = MAXDSZ - 1

            tmpf = k2f(s_atmpk[aoffset1:aoffset2])  # Stored in K
            mph = s_asmps[aoffset1:aoffset2] * 2.0
            psun = 100. - s_askyc[aoffset1:aoffset2]
            phour = s_ap01m[aoffset1:aoffset2] / 25.4  # Convert to inches
            relh = s_arelh[aoffset1:aoffset2]
            high = k2f(s_high[coffset]) + high_off
            low = k2f(s_low[coffset]) + low_off
            p01d = s_p01d[coffset] / 25.4  # Convert to inches
            # Apply the delta
            if p01d > 0.01:
                p01d += DELTAP[now.month - 1] * p01d

            # we smear the temperature data
            tmpf = hourly_fitter_temp(tmpf, low, high - low)
            tmpf, p_tmpf = boundschk('t', tmpf, p_tmpf, -20., 120.)

            # we smear the precipitation data
            phour = hourly_fitter_precip(phour, p01d)
            # 0,10 inches
            phour, p_phour = boundschk('p', phour, p_phour, 0.0, 10.)
            # if p01d > 4:
            #  print phour, p01d

            # can't touch these
            mph, p_mph = boundschk('m', mph, p_mph, 0.0, 100.)
            psun, p_psun = boundschk('s', psun, p_psun, 0.0, 100.)
            relh, p_relh = boundschk('r', relh, p_relh, 0.0, 100.)  # 0,100 %

            for i in range(len(tmpf)):
                ts = now + datetime.timedelta(hours=i)
                ds[ts] = [tmpf[i], mph[i],
                          psun[i], phour[i], relh[i]]
            now += interval
        # Now we need to do some smoothing between the periods to prevent
        # wild jumps in the data, due to oddities with leap day, we need to
        # merge both 1 January and 31 December, tricky!
        for hroffset in [0, -24]:
            for year in range(2011, 2050):
                sts = datetime.datetime(year, 1, 1, 7)
                sts += datetime.timedelta(hours=hroffset)
                t1 = sts + datetime.timedelta(hours=-1)
                deltat = (ds[sts][0] - ds[t1][0]) / 2.
                midpt = ds[sts][0] - deltat
                rdelta = midpt - ds[sts][0]
                ldelta = midpt - ds[t1][0]
                oldt = []
                newt = []
                deltas = []
                for hr in range(-36, 36):
                    ts = sts + datetime.timedelta(hours=hr)
                    delta = ldelta if hr < 0 else rdelta
                    newv = ((ds[ts][0] + delta * (1. / (abs(hr)+1)))
                            if hr != 0 else midpt)
                    if DEBUG:
                        oldt.append(ds[ts][0])
                        newt.append(newv)
                        deltas.append(newt[-1] - ds[ts][0])
                        print(("%s %5.2f [%5.2f] %5.2f %5.2f %5.2f %5.2f"
                               ) % (ts.strftime("%Y%m%d%H"), oldt[-1],
                                    newt[-1],
                                    ds[ts][1],
                                    ds[ts][2], ds[ts][3], ds[ts][4]))
                    # Assign the new value
                    ds[ts][0] = newv
                    # print("%s -> %s %s" % (ts, newv, ds[ts][0]))
                if DEBUG:
                    (fig, ax) = plt.subplots(2, 1, sharex=True)
                    ax[0].set_title("Example Temp Smoother to Merge Periods")
                    ax[0].set_ylabel("Air Temperature $^\circ$F")
                    ax[0].bar(np.arange(-36, 36)-0.4, oldt, width=0.8,
                              facecolor='tan',
                              label='Old', zorder=1, edgecolor='tan')
                    ax[0].scatter(range(-36, 36), newt, marker='o', color='b',
                                  s=40, label='New', zorder=2)
                    ax[0].set_xlim(-36, 36)
                    ax[0].set_xticks(np.arange(-36, 36, 6))
                    ax[0].grid(True)
                    ax[0].set_ylim(min(oldt)-3, max(oldt)+3)
                    ax[0].legend(loc='best')
                    bars = ax[1].bar(np.arange(-36, 36)-0.4, deltas, width=0.8,
                                     facecolor='b', edgecolor='b')
                    for i, bar in enumerate(bars):
                        if deltas[i] > 0:
                            bar.set_facecolor('r')
                            bar.set_edgecolor('r')
                    ax[1].grid(True)
                    ax[1].set_ylabel("Delta Applied $^\circ$F")
                    ax[1].set_xlabel("Hours around 1 January 12 AM")
                    fig.savefig('test-%s.png' % (sts.strftime("%Y%m%d"),))
                    plt.close()
                    print("--------------")
        # Now the magic to rewite these data
        out = open("%s.hcd" % (stid,), 'w')
        realts = datetime.datetime(2010, 1, 1)
        for (sts, ets) in PERIODS:
            now = sts
            print 'now is %s realts is %s' % (now, realts)
            while now <= ets:
                s = fmt % tuple(ds[now])
                if s.find('nan') > 0:
                    print('NaN Found in output: %s, debugging' % (s,))
                    pdb.set_trace()
                out.write(realts.strftime("%Y%m%d%H") + s)
                now += datetime.timedelta(hours=1)
                realts += datetime.timedelta(hours=1)
        out.close()
        subprocess.call("python qc_hcd.py %s.hcd" % (stid,), shell=True)
        # return

if __name__ == '__main__':
    runner()
