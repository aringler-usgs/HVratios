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
                            if float(data[7]) > .8:
                                azis.append(float(data[6]))
                                NV.append(float(data[10]))
                                EV.append(float(data[11]))
                                HV.append(float(data[8]))
                    except:
                        pass
    
    return HV, azis
    
stas = glob.glob('Results/*')
stas = [i.split('_')[1] for i in stas]
stas = list(set(stas))


#stas = ['HIA']
f0s = [0.00666666666667, 1./100., 0.0133333333333, 1./50., 1./25]
alpha = ['a','b','c','d','e']
#f0s = [1./150.]
print('Here is the length of stas:' + str(len(stas)))



from scipy import odr

def  fg(B, x):
    return B[0]*x + B[1]

linear = odr.Model(fg)

fig = plt.figure(1, figsize=(16,10))
for idx, f in enumerate(f0s):
    m00HV =[]
    std00HV=[]
    m10HV = []
    std10HV = []
    
    if idx == 0:    
        ax1 = plt.subplot2grid(shape =(2,6), loc =(0,0), colspan=2, aspect='equal')
    if idx ==1:
        ax1 = plt.subplot2grid(shape =(2,6), loc =(0,2), colspan=2, aspect='equal')
    if idx ==2:
        ax1 = plt.subplot2grid(shape =(2,6), loc =(0,4), colspan=2, aspect='equal')
    if idx ==3:
        ax1 = plt.subplot2grid(shape =(2,6), loc =(1,1), colspan=2, aspect='equal')
    if idx ==4:
        ax1 = plt.subplot2grid(shape =(2,6), loc =(1,3), colspan=2, aspect='equal')
    for sta in stas:
        HV00,azis = computeA(f, '00', sta)
        HV00 = np.asarray(HV00)
        #HV00 = HV00[(HV00 > 0.) & (HV00 <=3.)]
        
        #HV00 = HV00[(HV00>= 2.*np.mean(HV00))]
        HV10,azis = computeA(f, '10', sta)
        HV10 = np.asarray(HV10)
        #HV10 = HV10[(HV10 > 0.) & (HV10 <= 3.)]
        #HV10 = HV10[(HV10>= 2.*np.mean(HV10))]
        
        if ((not np.isnan((np.mean(HV00)))) and ( not np.isnan((np.mean(HV10))))):
            m00HV.append(np.mean(HV00))
            std00HV.append(np.std(HV00))
            m10HV.append(np.mean(HV10))
            std10HV.append(np.std(HV10))
    std10HV = np.asarray(std10HV)
    std10HV *= 0.5
    std00HV = np.asarray(std00HV)
    std00HV *= 0.5
    plt.errorbar(m00HV,m10HV, std00HV, std10HV,fmt='o', color='k', label='00 vs 10 Ratio')
    mydata = odr.Data(m00HV, m10HV)
    myodr = odr.ODR(mydata, linear, beta0=[1., 0.])
    myoutput = myodr.run()
    pfit = fg(myoutput.beta, np.asarray([-2.,2.]))
    ppara = myoutput.beta
    print('Here is the frequency band:' + str(f))
    print(pfit)
    print('Mean:' + str(np.mean(np.asarray(m00HV)-np.asarray(m10HV))))
    print('std:' + str(np.std(np.asarray(m00HV)-np.asarray(m10HV))))
    myoutput.pprint()
    print('New table 2 info')

    #ppara = np.polyfit(m00HV,m10HV,1)
    #pfit = np.poly1d(np.polyfit(m00HV,m10HV,1))
    plt.plot([-2.,2.], pfit, color ='0.7', label='Slope: ' + str(round(ppara[0],2)) + ' Intercept: ' + str(round(ppara[1],2)), linewidth=2.5)
    plt.title(str(int(round(1./f,0))) + ' s Period')
    #plt.legend()
    #plt.xlabel('Primary Sensor')
    #plt.ylabel('Secondary Sensor')
    #plt.title('H/V Ratio Uncertainty ' + str(int(round(1./f,0))) + ' s Period')
    plt.ylim((-0.6,.6))
    plt.xlim((-0.6,.6))
    plt.yticks((-0.3,0.0, 0.3))
    plt.xticks((-0.3,0.0,0.3))

    #plt.show()
    #plt.clf()
    del m00HV, std00HV, m10HV, std10HV

    #plt.savefig('Difference_' + str(1./f) + '.jpg', format='JPEG', dpi=400)
    #plt.clf()

#plt.show()

fig.text(0.08, 0.5, 'Sensor 10 H/V Ratio', ha = 'center', va = 'center', rotation = 'vertical', fontsize = 20)
fig.text(0.5, 0.05, 'Sensor 00 H/V Ratio', ha = 'center', va = 'center', rotation = 'horizontal', fontsize = 20)

#plt.tight_layout()

xs =[.12, 0.385, 0.655, 0.255, 0.52]
ys=[.89 ,.89, .89, 0.47, 0.47 ]
letters =['a','b','c','d','e']

for triple in zip(xs, ys, letters):
    plt.text(triple[0], triple[1], '(' + triple[2] + ')', fontsize=24, transform=plt.gcf().transFigure)





plt.savefig('Difference_PLOT.jpg', dpi= 400, format='JPEG')    
    
    
    #

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
print(np.mean(np.asarray(m00HV)-np.asarray(m10HV)))
