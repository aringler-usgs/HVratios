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

def filtergauss(st, f0):
    # f0 is the central frequency
    # alpha is our spread
    band = 2.5
    alpha = 3./band**2
    for tr in st:
        nfft = 2**(math.ceil(math.log(tr.stats.npts,2)))
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

def calcSNR(st, Noisest):
    snrR = np.sum((st.select(component="R")[0].data)**2)/np.sum((Noisest.select(component="R")[0].data)**2)
    snrV = np.sum((st.select(component="Z")[0].data)**2)/np.sum((Noisest.select(component="Z")[0].data)**2)
    SNR = min([snrR, snrV])
    return 10.*np.log10(SNR)

stime = UTCDateTime('2001-001T00:00:00.0')
etime = UTCDateTime('2018-001T00:00:00.0')

net = 'IC'
debug = True
window= 60.*8.
#f0s = [1./150., 1./100., 1./75., 1./50., 1./25]
f0s = [1./150.]
plots = False
stalist = False


client = Client("IRIS")


# Grab the list of station
sp = Parser('/APPS/metadata/SEED/' + net + '.dataless')



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
    fstring = 'Results/' + net + '_' + sta + '_' + str(eve.origins[0].time.year) + '_' + str(eve.origins[0].time.julday).zfill(3) + \
            '_' + str(eve.origins[0].time.hour).zfill(2) + '_' + str(eve.origins[0].time.minute).zfill(2) + \
            '_Results.csv'
    feve = open(fstring,'w')
    feve.write('sta, loc, year, day, f0, distance, azimuth, SNR, lag, corr, mag, HV, HVnolag, pNV, pNVstd, pEV, pEVstd, MeanZetaNR, STDZetaNR, MeanZetaER, STDZetaER, DeltaPhase, NV, EV, pRV, pTV \n')
    if debug:
        print('Distance: ' + str(dis))
        print('Azimuth: ' + str(azi))
    # compute arrival start and end times 620 s window
    astime = eve.origins[0].time + int(dis/4.)-30.
    aetime = astime + window-30.
    if debug:
        print('Here is astime:' + str(astime))
        print('Here is aetime:' + str(aetime))
        print('Here is our time window:' + str(window))
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
            if debug:
                print('Grabbing the noise data')
            # Here is data from an hour before
            Noisest = grabdata(eve.origins[0].time-60.*60., eve.origins[0].time-60.*60.+window, bazi, sp, sta, net, loc)
        except:
            print('No data for: ' + sta)
            continue
        if debug:
            print(Noisest)
            print(st)
        
        for f0 in f0s:

            st2 = finalfilter(st,f0,bazi,astime,aetime,True)
            # Here is our data in the ZNE directions so no rotation
            st2ZNE = finalfilter(st,f0,0.,astime,aetime, False)
            
            Noisest2 = finalfilter(Noisest,f0,bazi,eve.origins[0].time-60.*60.,eve.origins[0].time-60.*60.+window, True)
            if debug:
                print(st2)
            SNR = calcSNR(st2, Noisest2)
            del Noisest2
            if debug:
                print('Here is the SNR:' + str(SNR))
            if st2[0].stats.npts < .9*window:
                continue
            # We now have a good event with high SNR
            # Here is the phase shifted Vertical
            HilbertV = np.imag(hilbert(st2.select(component="Z")[0].data))
            lag, maxc, corr = xcorr(HilbertV,st2.select(component="R")[0].data, 100, full_xcorr=True)
            oldx = np.asarray(range(len(corr)))/(float(len(corr))/float(len(HilbertV)))

            corr = np.interp(range(len(HilbertV)),oldx, corr)

            deltaphi = 360.*lag*f0

            env = envelope(st2.select(component="R")[0].data)*envelope(HilbertV)
            env *= 1./np.max(np.abs(env))
            fig = plt.figure(1)
            #plt.subplot(211)
            #plt.plot(HilbertV, label='Vertical')
            #plt.plot(st2.select(component="R")[0].data, label='Radial')
            #plt.plot(st2.select(component="T")[0].data, label='Transverse')
            #plt.subplot(212)
            #plt.plot(env, label='Envelope')
            #plt.plot(corr, label='Coorelation')
            #plt.plot(env*corr, label='Function')
            #plt.legend()
            #plt.show()
            #plt.clf()
            
            
            
            
            # Now we want to compute the frequency domain pieces
            pV = np.fft.rfft(st2ZNE.select(component="Z")[0].data)
            pN = np.fft.rfft(st2ZNE.select(component="N")[0].data)
            pE = np.fft.rfft(st2ZNE.select(component="E")[0].data)
            pR = np.fft.rfft(st2.select(component="R")[0].data)
            pT = np.fft.rfft(st2.select(component="T")[0].data)
            
            
            freqs = np.fft.rfftfreq(st2[0].stats.npts)
            # Swtich to acceleration 
            pV *= (2.*np.pi*freqs*1j)**2
            pN *= (2.*np.pi*freqs*1j)**2
            pE *= (2.*np.pi*freqs*1j)**2
            pR *= (2.*np.pi*freqs*1j)**2
            pT *= (2.*np.pi*freqs*1j)**2
            # One-octave band freqmin
            freqmin = f0/np.sqrt(2.)
            freqmax = f0*np.sqrt(2.)
            pV = pV[(freqs <= freqmax) & (freqs >= freqmin)]
            pN = pN[(freqs <= freqmax) & (freqs >= freqmin)]
            pE = pE[(freqs <= freqmax) & (freqs>= freqmin)]
            pT = pT[(freqs <= freqmax) & (freqs >= freqmin)]
            pR = pR[(freqs <= freqmax) & (freqs>= freqmin)]
            freqs = freqs[(freqs <= freqmax) & (freqs >= freqmin)]
            ENV = env*corr
            corr = maxc
            #HV = np.sqrt(sum(st2.select(component="R")[0].data**2)/sum(HilbertV**2))
            H = st2.select(component="R")[0].data
            
            H = H[(ENV > .9)]**2
            V = HilbertV[(ENV >.9)]**2
            print(H)
            print(V)
            HV = np.sqrt(np.mean(H/V))
            
            
            
            
            
            
            
            
            trtemp = st2.select(component="R")[0].copy()
            trtemp.stats.starttime += -lag
            HVnolag = np.sqrt(sum(trtemp.data**2)/sum(HilbertV**2))

            pNV = 20.*np.log10(np.abs(pN)/np.abs(pV))
            pNVstd = np.std(pNV)
            pNV = np.mean(pNV)
            pEV = 20.*np.log10(np.abs(pE)/np.abs(pV))
            pEVstd = np.std(pEV)
            pEV = np.mean(pEV)
            pTV = 20.*np.log10(np.abs(pR)/np.abs(pV))
            pRV = 20.*np.log10(np.abs(pR)/np.abs(pV))
            pRV = np.mean(pRV)
            pTV = np.mean(pTV)
            
            
            
            N = ((st2ZNE.select(component='N')[0].data)[(ENV>.9)])**2
            V = ((st2ZNE.select(component='Z')[0].data)[(ENV>.9)])**2
            E = ((st2ZNE.select(component='E')[0].data)[(ENV>.9)])**2
            NV= np.sqrt(np.mean(N/V))
            EV = np.sqrt(np.mean(E/V))
            
            #NV = np.sqrt(sum(st2ZNE.select(component="N")[0].data**2)/sum(HilbertV**2))
            #EV = np.sqrt(sum(st2ZNE.select(component="E")[0].data**2)/sum(HilbertV**2))
            
            del trtemp
        
            # Now we compute zeta for NS and EW
            zetaN = pN*np.conj(pV)/(np.abs(pV)**2)
            zetaNR = -zetaN.real*(2.*np.pi*freqs)/9.81
            zetaE = pE*np.conj(pV)/(np.abs(pV)**2)
            zetaER = -zetaE.real*(2.*np.pi*freqs)/9.81
            
            feve.write(sta + ', ' + loc + ', ' + str(eve.origins[0].time.year) +', ' + 
                            str(eve.origins[0].time.julday).zfill(3) + ', ' + str(f0) + ', ' + 
                            str(disdeg) + ', ' + str(azi) + ', ' + str(SNR)
                            + ', ' + str(lag) + ', ' + str(corr) + ', ' + str(eve.magnitudes[0].mag)
                            + ', ' + str(HV) + ', ' + str(HVnolag) + ', '  
                            + str(pNV) + ', ' + str(pNVstd) + ', ' + str(pEV) + ', ' + str(pEVstd)
                            + ', ' + str(np.mean(zetaNR)) + ', ' + str(np.std(zetaNR))
                            + ', ' + str(np.mean(zetaER)) + ', ' + str(np.std(zetaER))
                            + ', ' + str(deltaphi) + ', ' + str(NV) + ', ' + str(EV) + ', ' + str(pRV) + ', ' + str(pTV) + ' \n')
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

#from multiprocessing import Pool

#pool = Pool(20)

#cat = cat[:2]
for idx, eve in enumerate(cat):
    doubles = []
    print('One event: ' + str(idx) + ' of ' + str(len(cat)))
    for sta in stations:
        doubles.append([eve,sta])
    for double in doubles:
        multifun(double)
    #pool.map(multifun,doubles)
    





        
