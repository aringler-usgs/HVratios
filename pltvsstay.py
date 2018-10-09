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
mpl.rc('font',size=20)


def computeA(freq, loc, station, debug = False):
    files = glob.glob('NEWRESULTS/*' + station + '*Results.csv')
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
                            if ((float(data[8]) != 0.) & (float(data[9]) != 0.)):
                                azis.append(float(data[6]))
                                NV.append(float(data[10]))
                                EV.append(float(data[11]))
                                HV.append(float(data[8]))
                    except:
                        pass
    
    return HV



    


letts = ['(a)', '(b)', '(c)', '(d)', '(e)']
fig =plt.figure(1, figsize=(12,18))
plt.subplots_adjust(wspace=0.001)
for idx, freq in enumerate([0.00666666666667, 1./100., 0.0133333333333, 1./50., 1./25]):
    diff=[]
    HVS00 = []
    HVS10 =[]
    HVSstd00=[]
    HVSstd10=[]
    for sta in stas:
        HV00 = computeA(freq, '00', sta)
        HV10 = computeA(freq, '10', sta)
        HV00 = np.asarray(HV00)
        HV00 = HV00[(HV00 >= -0.6) & (HV00 <=0.6)]
        HV10 = np.asarray(HV10)
        HV10 = HV10[(HV10 >= -0.6) & (HV10 <= 0.6)]
        HVS00.append(np.mean(HV00))
        HVS10.append(np.mean(HV10))
        HVSstd00.append(np.std(HV00))
        HVSstd10.append(np.std(HV10))
        diff.append(np.mean(HV00)- np.mean(HV10))
    diff = np.asarray(diff)
    print('Here is the mean ' + str(np.nanmean(diff)) + ' here is the std ' + str(np.nanstd(diff)))
    print(freq)    
    plt.subplot(151+idx)
    plt.errorbar(HVS00,np.asarray(range(len(stas)))+0.1,fmt='o',xerr=HVSstd00,label='00', alpha=0.7)
    plt.errorbar(HVS10, np.asarray(range(len(stas)))-0.1,fmt='o',xerr=HVSstd10, label='10',alpha=.7)
    plt.text(-0.45, len(stas)+.5, letts[idx] + ' ' + str(int(round(1/freq,0))) + ' s')
    plt.ylim((-0.5, len(stas)+2.5))
    plt.xticks([-0.3, 0.0 , 0.3])
    plt.xlim((-0.5,0.5))
    #plt.xlim((0.,2.))
    #plt.xticks([0.5, 1., 1.5])
    if idx == 0:
        plt.yticks(range(len(stas)), stas, fontsize=12)
    else:
        plt.yticks([])
    if idx == 2:
        plt.xlabel('H/V Ratio')


plt.subplot(153)
plt.legend(ncol=2, loc='center', bbox_to_anchor=(0.5, -.1))


plt.savefig('AllHVratios.jpg',format='JPEG',dpi=400)
#plt.show()
