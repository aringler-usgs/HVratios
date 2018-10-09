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
    

sta = 'ANMO'

files = glob.glob('AZIS2/*' + sta +'*')
res = []
for curfile in files:
   res += procfile(curfile)
   
print(res)

fig = plt.figure(1,figsize=(12,12))
plt.subplots_adjust(hspace=0.001)
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
            thetas.append(ele['theta']-0.1)
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
    plt.subplot(5,1,idx+1)
    
    
    plt.ylim((-.3,.6))
    plt.xlim((-5,5))
    #plt.xticks([-.4, 0, .4])
    if idx <= 4:
        plt.tick_params(labelbottom=False)   
    if idx == 4:
        plt.xlabel('Horizontal 00 to 10 Ratio')
        plt.plot(thetas, ARs,'.', label='Radial',markersize= 10, alpha=0.3)
        plt.plot(thetas, ATs,'.', label='Transverse',markersize= 10, alpha=0.3)
        plt.plot(thetas, AZs,'.', label='Vertical', markersize=10, alpha=0.3)
        
        
        
    else:
        plt.plot(thetas, ARs,'.', markersize=10, alpha=0.3)
        plt.plot(thetas, AZs,'.', markersize=10, alpha=0.3)
        plt.plot(thetas, ATs,'.', markersize= 10, alpha=0.3)
    if idx == 2:
        plt.ylabel('00 to 10 Ratio')
    plt.text(-4.8, 0.48,'(' + letters[idx] + ') ' + str(int(round(1./f,0))) + ' s')
ax=plt.subplot(5,1,5)        
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -.14),
        fancybox=False, shadow=False, ncol=2,fontsize=18)
plt.title('Relative Amplitude Ratios ' + sta)
#plt.savefig(sta + 'REL.jpg',format='JPEG',dpi=400)
plt.show()    
