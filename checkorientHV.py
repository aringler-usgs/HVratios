#!/usr/bin/env python

import glob
import matplotlib.pyplot as plt
import numpy as np
import sys

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
            with open(curfile.replace('AZIS2/','NEWRESULTS/'),'r') as f2:
                next(f2)
                for line in f2:
                    if float(line.split(',')[4]) == temp['f0']:
                        temp['HV' + (line.split(',')[1]).replace(' ','')] = float(line.split(',')[8])
            
            res.append(temp)
    
    
    return res
    



#def computeA(freq, loc, station, debug = False):
    #files = glob.glob('Results/*' + station + '*Results.csv')
    #azis = []
    #HV = []
    #EV =[]
    #NV =[]
    #for curfile in files:
        #nl = sum(1 for line in open(curfile))
        #if True:
        ##if nl >= 11:
            #with open(curfile,'r') as f:
                #for line in f:
           
                    #try:
                        #if (loc in line.split(',')[1])  and (freq == float(line.split(',')[4])):
                            #data = line
                            #data = data.split(',')
                            ##if float(data[7]) > .8:
                            #if True:
                                #azis.append(float(data[6]))
                                #NV.append(float(data[10]))
                                #EV.append(float(data[11]))
                                #HV.append(float(data[8]))
                    #except:
                        #pass
    
    #return HV, azis













stas = glob.glob('AZIS2/*')
stas = list(set([sta.split('_')[1] for sta in stas]))

print(stas)

#sta = 'ANMO'

for sta in stas:
    if (sta == 'WAKE') or (sta =='RSSD'):
        continue
    print(sta)
    files = glob.glob('AZIS2/*' + sta +'*')
    res = []
    for curfile in files:
       res += procfile(curfile)
       
    #print(res)
    
    
    
    
    
    
    
    fig = plt.figure(1,figsize=(16,10))
    #plt.subplots_adjust(hspace=0.001)
    letters = ['a','b','c','d','e']
    
    for idx, f in enumerate([0.00666666666667, 0.01,  0.0133333333333, 0.02, 0.04]):
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
        
        
        HV00s = []
        HV10s = []
        thetas = []
        for ele in res:
            if ele['f0'] == f:
                HV00s.append(ele['HV00'])
                HV10s.append(ele['HV10'])
                #AZs.append(np.log10(ele['AZ']))
                thetas.append(ele['theta']-0.1)
        HV00s = np.asarray(HV00s)
        HV10s = np.asarray(HV10s)
        HV00s[HV00s >= 0.6] = 0.6
        HV00s[HV00s <= -0.6] = -.6
        HV10s[HV10s >= 0.6] = 0.6
        HV10s[HV10s <= -.6] = -.6
        thetas = np.asarray(thetas)
        plt.ylim((-.6,.6))
        plt.xlim((-0.6,0.6))
        plt.xticks([-.4, 0, .4])
        plt.yticks([-.4,0., .4])

        s = plt.scatter(HV00s, HV10s,c=thetas)

        plt.text(-0.55, 0.5,'(' + letters[idx] + ') ' + str(int(round(1./f,0))) + ' s')
    
    fig.text(0.08, 0.5, 'Sensor 10 H/V Ratio', ha = 'center', va = 'center', rotation = 'vertical', fontsize = 20)
    fig.text(0.5, 0.05, 'Sensor 00 H/V Ratio', ha = 'center', va = 'center', rotation = 'horizontal', fontsize = 20)    

    cbar_ax=fig.add_axes([.90, 0.35, 0.04, 0.3])    
    cbar = fig.colorbar(s, ticks=[-4,  -2,  0,  2,  4], cax=cbar_ax, orientation='vertical', label='Relative Orientation (degrees)')
    fig.suptitle('H/V Ratio and Relative Orientation ' + sta)
    plt.savefig(sta + '_Relative.pdf',format='PDF', dpi=400)
    plt.clf()
    plt.close()
  
    
