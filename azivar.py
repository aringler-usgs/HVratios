#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys

stas = glob.glob('Results/*.csv')
stas = [sta.split('_')[1] for sta in stas]
stas = list(set(stas))
stas.sort()

print(stas)

import matplotlib as mpl
#mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=12)


def computeA(freq, loc, station, debug = False):
    files = glob.glob('Results/*' + station + '*Results.csv')
    azis = []
    HV = []
    EV =[]
    NV =[]
    for curfile in files:
        nl = sum(1 for line in open(curfile))
        if True:
        #if nl >= 11:
            with open(curfile,'r') as f:
                for line in f:
           
                    try:
                        if (loc in line.split(',')[1])  and (freq == float(line.split(',')[4])):
                            data = line
                            data = data.split(',')
                            #if float(data[7]) > .8:
                            if True:
                                azis.append(float(data[6]))
                                NV.append(float(data[10]))
                                EV.append(float(data[11]))
                                HV.append(float(data[8]))
                    except:
                        pass
    
    return HV, azis



    


letts = ['(a)', '(b)', '(c)', '(d)', '(e)']
fig =plt.figure(1, figsize=(12,12))
plt.subplots_adjust(wspace=0.001)
for idx, freq in enumerate([0.00666666666667, 1./100., 0.0133333333333, 1./50., 1./25]):

    HV00ang =[]
    HV10ang =[]
    sta00idx =[]
    sta10idx =[]
    angs00=[]
    angs10=[]
    for gidx, sta in enumerate(stas):
        HV00, azis00 = computeA(freq, '00', sta)
        HV10, azis10 = computeA(freq, '10', sta)
        HV00 = np.asarray(HV00)
        HV10 = np.asarray(HV10)
        azis00 = np.asarray(azis00)
        azis10 = np.asarray(azis10)
        azis00 = azis00[(HV00>0.) & (HV00 <= 2.)]
        azis10= azis10[(HV10>0.) & (HV10 <= 2.)]
        HV00 = HV00[(HV00 > 0.) & (HV00 <=2.)]
        HV10 = np.asarray(HV10)
        HV10 = HV10[(HV10 > 0.) & (HV10 <= 2.)]
        #print(HV10)
        #print(azis10)
        for ang in range(0,360, 30):

            HV00t = HV00[(azis00 > ang) & ( azis00 <= (ang+ 30.))]  
            HV10t = HV10[(azis10 > ang) & ( azis10 <= (ang+ 30.))] 

            HV00ang.append(np.mean(HV00t))
            HV10ang.append(np.mean(HV10t))
            sta00idx.append(gidx)
            sta10idx.append(gidx)
            angs00.append(ang)
            angs10.append(ang)
    sta00idx = np.asarray(sta00idx)
    sta10idx = np.asarray(sta10idx)
    plt.subplot(151+idx)
    #sc =plt.scatter(angs00, sta00idx,c=HV00ang)
    #
    sc =plt.scatter(angs10, sta10idx,c=HV10ang)
    plt.text(0.25, len(stas)+.5, letts[idx] + ' ' + str(int(round(1/freq,0))) + ' s')
    plt.ylim((-0.5, len(stas)+2.5))
    #plt.xlim((0.,2.))
    plt.xticks([0., 120., 240.])
    if idx == 0:
        plt.yticks(range(len(stas)), stas)
    else:
        plt.yticks([])
    if idx == 2:
        plt.xlabel('Back-Azimuth (degrees)')
    if idx == 4:
        cb=plt.colorbar(sc, ticks=[0.5, 1., 1.5])
        cb.set_label('H/V Ratio')

plt.subplot(153)
#plt.legend(ncol=2, loc='center', bbox_to_anchor=(0.5, -.1))
#plt.show()

plt.savefig('AllHVratiosByAngle10.jpg',format='JPEG',dpi=400)
