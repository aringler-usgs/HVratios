#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import glob

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
                            if float(data[7]) > .8:
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
    HVS00 = []
    HVS10 =[]
    HVSstd00=[]
    HVSstd10=[]
    for sta in stas:
        HV00 = computeA(freq, '00', sta)
        HV10 = computeA(freq, '10', sta)
        HV00 = np.asarray(HV00)
        HV00 = HV00[(HV00 > 0.) & (HV00 <=3.)]
        HV10 = np.asarray(HV10)
        HV10 = HV10[(HV10 > 0.) & (HV10 <= 3.)]
        HVS00.append(np.mean(HV00))
        HVS10.append(np.mean(HV10))
        HVSstd00.append(np.std(HV00))
        HVSstd10.append(np.std(HV10))
    plt.subplot(151+idx)
    plt.errorbar(HVS00,range(len(stas)),fmt='o',xerr=HVSstd00,label='00')
    plt.errorbar(HVS10, range(len(stas)),fmt='o',xerr=HVSstd10, label='10')
    plt.text(0.25, len(stas)+.5, letts[idx] + ' ' + str(int(round(1/freq,0))) + ' s')
    plt.ylim((-0.5, len(stas)+2.5))
    plt.xlim((0.,2.))
    plt.xticks([0.5, 1., 1.5])
    if idx == 0:
        plt.yticks(range(len(stas)), stas)
    else:
        plt.yticks([])
    if idx == 2:
        plt.xlabel('H/V Ratio')


plt.subplot(153)
plt.legend(ncol=2, loc='center', bbox_to_anchor=(0.5, -.1))


plt.savefig('AllHVratios.jpg',format='JPEG',dpi=400)
