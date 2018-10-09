#!/usr/bin/env python
import matplotlib.pyplot as plt
import glob
import numpy as np

import matplotlib as mpl
#Set font parameters using matplotlib
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

import glob


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
                            if ((float(data[8]) == 0.) & (float(data[9]) == 0.)):
                                continue
                            if (float(data[7]) > .8):
                                azis.append(float(data[6]))
                                NV.append(float(data[10]))
                                EV.append(float(data[11]))
                                HV.append(float(data[8]))
                    except:
                        pass
    
    return HV, azis
    
stas = glob.glob('NEWRESULTS/*')
stas = [i.split('_')[1] for i in stas]
stas = list(set(stas))

m00HV =[]
std00HV=[]
m10HV = []
std10HV = []
stas = ['POHA']
f0s = [0.00666666666667, 1./100., 0.0133333333333, 1./50., 1./25]
alpha = ['a','b','c','d','e']
#f0s = [1./150.]
print('Here is the length of stas:' + str(len(stas)))

for sta in stas:
    binval = np.linspace(-0.6,0.6,40)
    fig = plt.figure(1,figsize=(12,12))
    plt.subplots_adjust(hspace=0.001)
    for idx, f in enumerate(f0s):
        plt.subplot(len(f0s),1,idx +1)
        if idx == 0:
            plt.title(sta + ' H/V Ratio Histograms')
        HV00,azis = computeA(f, '00', sta)
        HV00 = np.asarray(HV00)
        sigma = np.std(HV00)
        #HV00 = HV00[(HV00 >= (-sigma + np.mean(HV00))) & (HV00 <= (sigma + np.mean(HV00)))]
        
        HV10,azis = computeA(f, '10', sta)
        
        HV10 = np.asarray(HV10)
        sigma = np.std(HV10)
        #HV10 = HV10[(HV10 >= (-sigma + np.mean(HV10))) & (HV10 <= (sigma + np.mean(HV10)))]
        plt.hist(HV00, bins=binval, label='00 Mean=' + str(round(np.mean(HV00),2)) + ' Median=' + str(round(np.median(HV00),2)), normed=True, edgecolor='black', color='C0')
        plt.hist(HV10, bins=binval, label='10 Mean=' + str(round(np.mean(HV10),2)) + ' Median=' + str(round(np.median(HV10),2)), alpha=.5, normed=True, edgecolor='black', color='C1')
        plt.text(-0.57, 8., '(' + alpha[idx] + ') ' + str(int(round(1/f,0))) + ' s')
        plt.xlim((-0.6,0.6))
        plt.ylim((0,10.))
        plt.yticks([2.0, 4.0, 6.0, 8.0])
        #plt.ylabel('Counts')
        if idx < len(f0s)-1:
            plt.xticks([])
        else:
            plt.xlabel('HV')
        if idx == 2:
            plt.ylabel('Distribution (Normalized)')
        plt.legend(loc=1)
    plt.savefig('Histogram_' + sta + '.jpg', format='JPEG', dpi=400)
    plt.show()
    plt.clf()
    #plt.show()
    
    
    
    #if ((not np.isnan((np.mean(HV00)))) and ( not np.isnan((np.mean(HV10))))):
        #m00HV.append(np.mean(HV00))
        #std00HV.append(np.std(HV00))
        #m10HV.append(np.mean(HV10))
        #std10HV.append(np.std(HV10))

    ##fig = plt.figure(1,figsize=(12,12))
    ##binval = np.linspace(0.,1.,40)
    ##plt.hist(HV00, bins=binval, label='00 Mean=' + str(round(np.mean(HV00),2)))

    
    
    ##binval = np.linspace(0.,1.,40)
    ##plt.hist(HV10, bins=binval, label='10 Mean=' + str(round(np.mean(HV10),2)), alpha=.5)
    ##plt.xlabel('HV')
    ##plt.ylabel('Counts')
    ##plt.title(sta)
    ##plt.legend()
    ##plt.show()

#m00HV = np.asarray(m00HV)
#m10HV = np.asarray(m10HV)

#m00HV[(m00HV >= 1.5)] = 1.5
#m10HV[(m10HV >= 1.5)] = 1.5

    
#fig = plt.figure(2)
#plt.plot(np.asarray(m00HV), np.asarray(m10HV),'.')
#plt.show()
#print(m00HV)
#print(m10HV)
#print(np.mean(np.asarray(m00HV)-np.asarray(m10HV)))
