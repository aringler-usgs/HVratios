#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import math
import os
from obspy.core import read, UTCDateTime
from obspy.clients.fdsn import Client
from obspy.io.xseed import Parser
from time import gmtime, strftime
from obspy.geodetics.base import gps2dist_azimuth
from scipy.signal import hilbert
from obspy.signal.cross_correlation import xcorr_max, xcorr
import scipy.signal
from obspy.signal.detrend import polynomial
from obspy.signal.filter import envelope
from scipy.stats import pearsonr

#Font parameters from matplotlib
import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=22)


def filtergauss(st, f0):
    # f0 is the central frequency
    # alpha is our spread
    band = 2.5
    alpha = 3./band**2
    for tr in st:
        nfft = int(2**(math.ceil(math.log(tr.stats.npts,2))))
        data = np.fft.rfft(tr.data, n=nfft)
        w = 2.*np.pi*np.fft.rfftfreq(nfft, tr.stats.delta)
        w0=2.*np.pi*f0
        H=np.zeros(len(w)) 
        for idx in range(len(w)):
            if (w[idx] >= (1-band)*w0) and (w[idx] <= (1+band)*w0):
                H[idx] =np.exp(-alpha*(w[idx]-w0)**2/(w0**2))
        data *=H
        data[-1] = abs(data[-1]) + 0.0j
        tr.data = np.fft.irfft(data)[0:tr.stats.npts]
    return st

def getstalist(sp, etime, net):
    """ A function to get a station list. """
    stations = []
    for cursta in sp.stations:
# As we scan through blockettes we need to find blockettes 50 
        for blkt in cursta:
            if blkt.id == 50:
# Pull the station info for blockette 50
                stacall = blkt.station_call_letters.strip()
                if debug:
                    print "Here is a station in the dataless: " + stacall
                if type(blkt.end_effective_date) is str:
                    curdoy = strftime("%j", gmtime())
                    curyear = strftime("%Y", gmtime())
                    curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
                    
                    if blkt.start_effective_date <= etime:
                        stations.append(blkt.station_call_letters.strip())
                elif blkt.start_effective_date <= etime and blkt.end_effective_date >= etime:
                    stations.append(blkt.station_call_letters.strip())    
    return stations


def grabdata(astime, aetime, bazi, sp, sta, net, loc, debug = False):
    # Check if astime
    path = '/msd/' + net + '_' + sta + '/'
    path += str(astime.year) + '/' + str(astime.julday).zfill(3)
    path += '/' + loc + '_LH*'
    st = read(path)
    # In case we have a year or day boundary
    if (astime.year < aetime.year) or (astime.julday < aetime.julday):
        path = '/msd/' + net + '_' + sta + '/'
        path += str(astime.year) + '/' + str(astime.julday).zfill(3)
        path += '/' + loc + '_LH*'
        st += read(path)
    st.merge(fill_value=0.)
    st.trim(astime-20., aetime+20.)
    st.detrend('linear')
    if debug:
        print(st)
    for tr in st:
        paz = sp.get_paz(tr.id, astime)

        tr.simulate(paz_remove=paz, nfft_pow2=True, taper=True)
    
        # Deal with microsecond timing issues
        tcorr = round(float(tr.stats.starttime.microsecond)/10.**6,1)
        if debug:
            print(tcorr)
        tr.stats.starttime += tcorr
        tr.stats.starttime += -float(tr.stats.starttime.microsecond)/10.**6
        if debug:
            print(tr.stats.starttime)
    if debug:
        print(st)
    st = sp.rotate_to_zne(st)
    return st
    
def finalfilter(st,f0,bazi,astime,aetime, rot, debug=False):
    st2 = st.copy()
    if f0 > 0.:
        st2 = filtergauss(st2,f0)
    for tr in st2:
        polynomial(tr.data, 3)
    if debug:
        print(st2)
    if rot:
        st2.rotate('NE->RT', back_azimuth=bazi)
    st2.trim(astime,aetime)
    st2.taper(0.05)
    st2.sort()
    return st2


def HVscheme2(st, f0, bazi, astime, aetime, disdeg, eve):
    
    st2 = finalfilter(st,f0,bazi,astime,aetime,True)
    corrs =[]
    win = int(round(1./f0,0))
    for window in st2.slide(window_length=win, step=int(round(win/16.,0))):
        HilbertV = np.imag(hilbert(window.select(component="Z")[0].data))
        lag, corr = xcorr(HilbertV,window.select(component="R")[0].data, 5, full_xcorr=False)
        #corr = pearsonr(HilbertV, window.select(component="R")[0].data)
        corrs.append(corr)
    corr = corrs
    HilbertV = np.imag(hilbert(st2.select(component="Z")[0].data))
    oldx = np.asarray(range(len(corr)))/(float(len(corr))/float(len(HilbertV)))
    corr = np.interp(range(len(HilbertV)),oldx, corr)
    env = envelope(st2.select(component="R")[0].data)*envelope(HilbertV)
    HV = envelope(st2.select(component="R")[0].data)/envelope(HilbertV)
    env *= 1./np.max(np.abs(env))
    #corr *= env
    
    
    t = np.asarray(range(len(HilbertV)))
    
    lim = t[(corr >= .90)]
    if len(lim) == 0:
        return 0., 0., 0.
    
    
    
    HV2 = HV[(corr >= .90)]
    lim = lim[(HV2 <= np.mean(HV2) + 3.*np.std(HV2)) & (HV2 >= np.mean(HV2) - 3.*np.std(HV2)) ]
    HV2 = HV2[(HV2 <= np.mean(HV2) + 3.*np.std(HV2)) & (HV2 >= np.mean(HV2) - 3.*np.std(HV2)) ]
    
    mHV = np.mean(HV2)
    stdHV = np.std(HV2)
    
    
    
    
    return mHV, stdHV,  corr



stime = UTCDateTime('2017-001T00:00:00.0')
etime = UTCDateTime('2018-001T00:00:00.0')

net = 'IC'
debug = True
window= 60.*8.
f0s = [1./150., 1./100., 1./75., 1./50., 1./25]
#f0s = [1./25.]
plots = False
stalist = False


client = Client("IRIS")


# Grab the list of station
sp = Parser('/APPS/metadata/SEED/' + net + '.dataless')

def makecolocplot(eve, sta, debug=False):
    
    
    
    return
    
    




def proceve(eve, sta, debug=False):
    
    try:
        coords = sp.get_coordinates(net + '.' + sta + '.00.LHZ', eve.origins[0].time)
    except:
        return
    (dis, azi, bazi) = gps2dist_azimuth(coords['latitude'], coords['longitude'], 
                                            eve.origins[0].latitude,eve.origins[0].longitude)
        
    # Now in km
    dis *= 1./1000.
    disdeg = dis*0.0089932
    # Check for events way outside of our interested window
    if disdeg <=50. or disdeg >=120.:
        return
    fstring = 'NEWRESULTS/' + net + '_' + sta + '_' + str(eve.origins[0].time.year) + '_' + str(eve.origins[0].time.julday).zfill(3) + \
            '_' + str(eve.origins[0].time.hour).zfill(2) + '_' + str(eve.origins[0].time.minute).zfill(2) + \
            '_Results.csv'
    feve = open(fstring,'w')
    feve.write('sta, loc, year, day, f0, distance, azimuth, corr, mag, Log(mHV), Log(stdHV), pRV, pNV, pEV \n')
    if debug:
        print('Distance: ' + str(dis))
        print('Azimuth: ' + str(azi))
    # compute arrival start and end times 620 s window
    astime = eve.origins[0].time + int(dis/4.)-30.
    aetime = astime + window-30.

    # Grab the data now have trimmed RT data in velocity with filter
    locs = glob.glob('/msd/' + net + '_' + sta + '/' + str(astime.year) + '/'
                    + '/' + str(astime.julday).zfill(3) + '/*LHZ*')
    locs = [(loc.split('/')[-1]).split('_')[0] for loc in locs]
    for loc in locs:
        try:
        #if True:
            if debug:
                print('Grabbing the event data')
            st = grabdata(astime, aetime, bazi, sp, sta, net, loc)
            
        except:
            print('No data for: ' + sta)
            continue
        if debug:
            print(Noisest)
            print(st)
        
        for f0 in f0s:

            st2 = finalfilter(st,0.,bazi,astime,aetime,True)
            # Here is our data in the ZNE directions so no rotation
            st2ZNE = finalfilter(st,0.,0.,astime,aetime, False)
            st2 += st2ZNE
            st2.merge()
            
            if st2[0].stats.npts < .9*window:
                continue
            # We now have a good event with high SNR
            
            
            st2 = finalfilter(st,f0,bazi,astime,aetime,True)













            corrs =[]
            win = int(round(1./f0,0))
            for window2 in st2.slide(window_length=win, step=int(round(win/16.,0))):
                HilbertV = np.imag(hilbert(window2.select(component="Z")[0].data))
                lag, corr = xcorr(HilbertV,window2.select(component="R")[0].data, 5, full_xcorr=False)
                #corr = pearsonr(HilbertV, window.select(component="R")[0].data)
                corrs.append(corr)
            corr = corrs
            HilbertV = np.imag(hilbert(st2.select(component="Z")[0].data))
            oldx = np.asarray(range(len(corr)))/(float(len(corr))/float(len(HilbertV)))
            corr = np.interp(range(len(HilbertV)),oldx, corr)

            HV = envelope(st2.select(component="R")[0].data)/envelope(HilbertV)
            TV = envelope(st2.select(component="T")[0].data)/envelope(HilbertV)
            
            t = np.asarray(range(len(HilbertV)))
            #def angfit()
            lim = t[(corr >= .90)]
            HV2 = HV[(corr >= .90)]
            TV2 = TV[(corr >= .90)]
            

            lim = lim[(HV2 <= np.mean(HV2) + 3.*np.std(HV2)) & (HV2 >= np.mean(HV2) - 3.*np.std(HV2)) ]
            TV2 = TV2[(HV2 <= np.mean(HV2) + 3.*np.std(HV2)) & (HV2 >= np.mean(HV2) - 3.*np.std(HV2)) ]
            HV2 = HV2[(HV2 <= np.mean(HV2) + 3.*np.std(HV2)) & (HV2 >= np.mean(HV2) - 3.*np.std(HV2)) ]
            
            #nm = max([max(HV2), max(TV2)])
            




            mHV = np.mean(HV2)
            stdHV = np.std(HV2)
            stdHVL = np.std(np.log10(HV2))
            try:
                fig = plt.figure(1, figsize=(12,12))
                plt.subplots_adjust(hspace=0.001)
                plt.subplot(211)
                plt.title(st[0].stats.network + ' ' + st[0].stats.station + ' ' + ' Period: ' + str(int(round(1./f0,0))) + ' s Distance: ' + str(round(disdeg,0)) + ' degrees')
                plt.plot(t, HilbertV*10**9,linewidth=2.5, label=' Shifted Vertical: ' + loc)
                plt.xlim((min(t),max(t)))
                plt.plot(t,st2.select(component="R")[0].data*10**9, linewidth=2.5, label='Radial: ' + loc)
                plt.ylabel('Velocity (nm/s)')
                plt.xticks([])
                plt.legend(loc=1)
                plt.subplot(212)
                plt.plot(t, HV, label=loc + ' Log(HV=' + str(round(np.log10(mHV),2)) + '$\pm$' + str(round(stdHVL,2)), linewidth=2.5)
                #plt.plot(t, corr, label=loc + ' Characteristic Function')
                plt.ylim((0., 2.))
                plt.yticks([0.,  1., 2.])
                plt.ylabel('HV Ratio')
                plt.axvspan(min(lim), max(lim), 0.,2.,alpha=.3, color='.5')
                plt.xlim((min(t),max(t)))
                plt.xlabel('Time (s)')
                plt.legend(loc=1)
            except:
                print('Problem continue')
            #plt.show()
            #plt.clf()
    
            
            
            
            
            
            #pV = np.fft.rfft(st2.select(component="Z")[0].data)
            #freqmin = f0/np.sqrt(2.)
            #freqmax = f0*np.sqrt(2.)
            #freqs = np.fft.rfftfreq(st2[0].stats.npts)
            #pV = pV[(freqs <= freqmax) & (freqs >= freqmin)]
            #feve.write(sta + ', ' + loc + ', ' + str(eve.origins[0].time.year) +', ' + 
                            #str(eve.origins[0].time.julday).zfill(3) + ', ' + str(f0) + ', ' + 
                            #str(disdeg) + ', ' + str(azi) + ', ' + str(eve.magnitudes[0].mag)
                            #+ ', ' + str(mHV) + ', ' + str(stdHV))
            #for comp in ['R', 'N', 'E']:
                #pC = np.fft.rfft(st2.select(component=comp)[0].data)
                #pC = pC[(freqs <= freqmax) & (freqs >= freqmin)]
                #pCV = np.mean(np.abs(pC)/np.abs(pV))
                #feve.write( ', ' + str(pCV) )
            #feve.write(' \n')
    plt.savefig('BOTHPLT_' + st[0].stats.network + '_' + st[0].stats.station + '_' + st[0].stats.location + '_' + str(eve.origins[0].time.year) + '_' + str(eve.origins[0].time.julday).zfill(3) + \
            '_' + str(eve.origins[0].time.hour).zfill(2) + '_' + str(eve.origins[0].time.minute).zfill(2) + '_' + str(int(round(1./f0,0))) + '.pdf', format='PDF', dpi=400)
    plt.clf()                        
    feve.close()        
    num_lines = sum(1 for line in open(fstring))
    if num_lines <= 1:
        os.remove(fstring)
    return


def multifun(double):
    proceve(double[0],double[1])
    return


cat = client.get_events(starttime=stime, endtime=etime, 
                                minmagnitude=6.5, maxdepth=50.)


stations = getstalist(sp, stime, net)
stations = ['HIA']

from multiprocessing import Pool

pool = Pool(20)

#cat = cat[:2]
for idx, eve in enumerate(cat):
    doubles = []
    print('One event: ' + str(idx) + ' of ' + str(len(cat)))
    for sta in stations:
        doubles.append([eve,sta])
    for double in doubles:
        multifun(double)
    #pool.map(multifun,doubles)
    





