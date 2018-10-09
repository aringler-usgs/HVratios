#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys

stas = glob.glob('Results/*.csv')
stas = [sta.split('_')[1] for sta in stas]
stas = list(set(stas))
stas.sort()
#stas =['ANMO','HIA']
#print(stas)

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
    
    return HV, azis



    


letts = ['(a)', '(b)', '(c)', '(d)', '(e)']
fig =plt.figure(1, figsize=(12,18))
plt.subplots_adjust(wspace=0.001)
for idx, freq in enumerate([0.00666666666667, 1./100., 0.0133333333333, 1./50., 1./25]):

    HV00max =[]
    HV00min =[]
    HV10max =[]
    HV10min =[]
    stadiffs =[]
    for gidx, sta in enumerate(stas):
        HV00, azis00 = computeA(freq, '00', sta)
        HV10, azis10 = computeA(freq, '10', sta)
        HV00 = np.asarray(HV00)
        HV10 = np.asarray(HV10)
        azis00 = np.asarray(azis00)
        azis10 = np.asarray(azis10)
        azis00 = azis00[(HV00>=-0.6) & (HV00 <= .6)]
        azis10= azis10[(HV10>=-0.6) & (HV10 <= .6)]
        HV00 = HV00[(HV00 >= -0.6) & (HV00 <=.6)]
        HV10 = np.asarray(HV10)
        HV10 = HV10[(HV10 >= -0.6) & (HV10 <= .6)]
        #print(HV10)
        #print(azis10)
        HV00ang=[]
        HV10ang=[]

        for ang in range(0,360, 30):

            HV00t = HV00[(azis00 > ang) & ( azis00 <= (ang+ 30.))]  
            HV10t = HV10[(azis10 > ang) & ( azis10 <= (ang+ 30.))] 
            
            HV00ang.append(np.mean(HV00t))
            HV10ang.append(np.mean(HV10t))
            stadiffs.append(np.mean(HV00t)-np.mean(HV10t))    
        HV00ang = np.asarray(HV00ang)
        HV10ang = np.asarray(HV10ang)
        HV00max.append(np.nanmax(HV00ang))
        HV00min.append(np.nanmin(HV00ang))
        HV10max.append(np.nanmax(HV10ang))
        HV10min.append(np.nanmin(HV10ang))
    stadiffs=np.asarray(stadiffs)
    
    print('Here is stadiffs:' + str(np.nanmean(stadiffs)) + ' ' + str(np.nanstd(stadiffs)) + ' ' + str(np.nanmax(np.abs(stadiffs))) + ' ' + str(np.nanmin(stadiffs)))
    print(freq)
    plt.subplot(151+idx)
    #sc =plt.scatter(angs00, sta00idx,c=HV00ang)
    #
    for hidx, ele in enumerate(stas):
        if hidx == 0:
            plt.plot([HV00min[hidx],HV00max[hidx]] , [hidx+0.1, hidx+0.1], label='00',alpha = 0.7, color='C0', marker='.',linewidth=3)
            plt.plot([HV10min[hidx],HV10max[hidx]] , [hidx-0.1, hidx-0.1], label='10', alpha=.7, color='C1', marker='.', linewidth=3)
        else:
            plt.plot([HV00min[hidx],HV00max[hidx]] , [hidx+0.1, hidx+0.1],  color='C0', alpha =0.7, marker='.', linewidth=3)
            plt.plot([HV10min[hidx],HV10max[hidx]] , [hidx-0.1, hidx-0.1],  alpha=.7, color='C1', marker='.', linewidth=3)
        #plt.plot(HV00min,range(len(stas)),'.',label='00')
        #plt.plot(HV10max,range(len(stas)),'.',label='10')
        #plt.plot(HV10min,range(len(stas)),'.',label='10')
    plt.text(-0.45, len(stas)+.5, letts[idx] + ' ' + str(int(round(1/freq,0))) + ' s')
    plt.ylim((-0.5, len(stas)+2.5))
    #plt.xlim((0.,2.))
    plt.xticks([-0.3, 0.0 , 0.3])
    plt.xlim((-0.5,0.5))
    if idx == 0:
        plt.yticks(range(len(stas)), stas, fontsize=12)
    else:
        plt.yticks([])
    if idx == 2:
        plt.xlabel('H/V Ratio')
    #if idx == 4:
        #cb=plt.colorbar(sc, ticks=[0.5, 1., 1.5])
        #cb.set_label('H/V Ratio')

plt.subplot(153)
plt.legend(ncol=2, loc='center', bbox_to_anchor=(0.5, -.1))
#plt.show()

plt.savefig('NEWminmaxazi.jpg',format='JPEG',dpi=400)
