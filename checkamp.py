#!/usr/bin/env python

import glob
import matplotlib.pyplot as plt
import numpy as np

#Font parameters from matplotlib
import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

def procfile(curfile):
    res =[]
    with open(curfile,'r') as f:
        next(f)
        
        for line in f:
            temp = {}
            line = line.split(',')
            temp['loc'] = line[1]
            temp['f0'] = float(line[4])
            temp['theta'] = float(line[8])
            temp['AR'] = float(line[10])
            temp['AT'] = float(line[12])
            temp['AZ'] = float(line[11])
            res.append(temp)
    return res
    

sta = 'WVT'

files = glob.glob('AZIS/*' + sta +'*')
res = []
for curfile in files:
   res += procfile(curfile)
   
print(res)

fig = plt.figure(1,figsize=(12,12))
plt.subplots_adjust(wspace=0.001)
letters = ['a','b','c','d','e']

for idx, f in enumerate([0.00666666666667, 0.01,  0.0133333333333, 0.02, 0.04]):
    ARs =[]
    AZs =[]
    ATs =[]
    thetas =[]
    for ele in res:
        if ele['f0'] == f:
            ARs.append(np.log10(ele['AR']))
            ATs.append(np.log10(ele['AT']))
            AZs.append(np.log10(ele['AZ']))
            thetas.append(ele['theta'])
    ARs = np.asarray(ARs)
    AZs = np.asarray(AZs)
    ATs = np.asarray(ATs)
    ARs[ARs >= .6] = .6
    ARs[ARs <= -0.6] = -.6
    AZs[AZs >= .6] = .6
    AZs[AZs <= -.3] = -.3
    ATs[ATs >= .6] = .6
    ATs[ATs <= -.6] = -.6
    print(thetas)
    plt.subplot(1,5,idx+1)
    
    
    plt.ylim((-.3,.6))
    plt.xlim((-.6,.6))
    plt.xticks([-.4, 0, .4])
    if idx > 0:
        plt.tick_params(labelleft=False)   
    if idx == 2:
        plt.xlabel('Horizontal 00 to 10 Ratio')
        plt.plot(ARs, AZs,'.', label='Radial',markersize= 10, alpha=0.3)
        plt.plot(ATs, AZs,'.', label='Transverse', markersize=10, alpha=0.3)
        
    else:
        plt.plot(ARs, AZs,'.', markersize=10, alpha=0.3)
        plt.plot(ATs, AZs,'.', markersize=10, alpha=0.3)
    if idx == 0:
        plt.ylabel('Vertical 00 to 10 Ratio')
    plt.text(-0.55, 0.57,'(' + letters[idx] + ') ' + str(int(round(1./f,0))) + ' s')
ax=plt.subplot(1,5,3)        
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -.14),
        fancybox=False, shadow=False, ncol=2,fontsize=18)
plt.title('Relative Amplitude Ratios ' + sta)
plt.savefig(sta + 'REL.jpg',format='JPEG',dpi=400)
#plt.show()    
