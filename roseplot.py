#!/usr/bin/env python
import matplotlib.pyplot as plt
import glob
import matplotlib as mpl
import numpy as np
import sys
from scipy.optimize import fmin
#from scipy.spatial import Delaunay
from matplotlib.collections import PolyCollection
##Import font parameters from matplotlib
#mpl.rc('font',family='serif')
#mpl.rc('font',serif='Times') 
#mpl.rc('text', usetex=True)
#mpl.rc('font',size=14)


import matplotlib as mpl
#mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=22)

letts = ['(a)', '(b)', '(c)', '(d)', '(e)']
fRes = open('Final_Results.csv','w')
stas = glob.glob('NEWRESULTS/*.csv')
stas = [sta.split('_')[1] for sta in stas]
stas = list(set(stas))
stas=['MAJO']
debug = True

plotheader = False
for gidx, sta in enumerate(stas):
    print('On sta ' + sta + ' ' + str(gidx) + ' of ' + str(len(stas)))
    #try:
    if True:


        def checkstdAzi(azis,HVs,loc,p0):
            # Get variability as a gunction of angle
            HVsT = HVs[(locs == loc) & (p0s == p0)]
            azisT = np.rad2deg(azis[(locs == loc) & (p0s == p0)])
            var=[]
            for ang in range(0, 360, 10):
                temp = HVsT[(azisT >= float(ang)) & (azisT < float(ang) + 10.)]
                if len(temp) > 0:
                    var.append(np.std(temp))
            return np.mean(var)
                
        def complocs(HVs,p0, times):
            t10=times[(locs=='10') & (p0s == p0)]
            HV00s= HVs[(locs=='00') & (p0s == p0)]
            t00=times[(locs=='00') & (p0s == p0)]
            HV10s = HVs[(locs=='10') & (p0s == p0)]
            tGs = [time for time in t10 if time in t00]
            HV00G=[]
            HV10G=[]
            Ts =[]
            for tG in tGs:
                Ts.append(tG)
                HV00G.append(HV00s[list(t00).index(tG)])
                HV10G.append(HV10s[list(t10).index(tG)])
            return np.asarray(HV00G), np.asarray(HV10G), np.asarray(Ts)
             

        strprint = ''
        strele = 'sta, '


        strprint += sta + ', '

        files = glob.glob('NEWRESULTS/*_' + sta + '*.csv')

        vals =[]
        numevens = 0
        for curfile in files:
            net= curfile.split('_')[0]
            net = net.split('/')[1]
            f = open(curfile,'r')
            if sum(1 for line in f) < 3:
                f.close()
                continue
            else:
                numevens += 1
                f.close()
                f = open(curfile,'r')
                for line in f:
                    if sta in line:
                        line="".join(line.split())
                        line = line.split(",")
                        if (float(line[8]) >= -2.) & (float(line[8]) <= 2.):
                            val = {'loc': line[1], 'year': int(line[2]), 'doy': int(line[3]), 'hour': int(curfile.split('_')[4]),  
                                    'minute': int(curfile.split('_')[5]), 'p0': int(round(1./float(line[4]))),
                                    'dis': float(line[5]), 'azi': float(line[6]), 
                                    'mag': float(line[7]), 'HV': float(line[8]), 'stdHV': float(line[9]), 'phase': float(line[1])}
                            if ((val['stdHV'] == 0.) & (val['HV'] == 0.)):
                                continue
                            else:
                                vals.append(val)
                f.close()
        if debug:
            print('Got all the data into a dictionary')
        getValues = lambda key,inputData: np.asarray([subVal[key] for subVal in inputData if key in subVal])

        strele += 'Num eve, '
        strprint += str(numevens) + ', '

        locs = getValues('loc',vals)
        
        p0s = getValues('p0', vals)
        
    

        strele += ' Good eve 00, Good eve 10, '
        


        locs = getValues('loc',vals)
        times = getValues('year', vals) + np.asarray([float(day)/365.2 for day in getValues('doy',vals)])      
        times += np.asarray([float(hour)/(24.*365.2) for hour in getValues('hour',vals)]) 
        times += np.asarray([float(minute)/(60.*24.*365.2) for minute in getValues('minute',vals)]) 


        HVs = getValues('HV',vals)
        azis = np.deg2rad(getValues('azi', vals))
        p0s = getValues('p0', vals)

        ps = [150, 100, 75, 50, 25]
        for p in ps: 
            v00=checkstdAzi(azis, HVs,'00', p)
            v10 =checkstdAzi(azis, HVs,'10', p)
            strele += '00 Var ' + str(p) + ', 10 Var ' + str(p) + ', '
            strprint += str(v00) + ', ' + str(v10) + ', '


        
            
            
        fig=plt.figure(1,figsize=(16,16))
        #plt.suptitle(net +' ' + sta)
        #plt.subplots_adjust(hspace=0.001)
        #for idx, per in enumerate([150]):
        for idx, per in enumerate([150, 100, 75, 50, 25]):

            ax1 = plt.subplot(111, projection="polar")
            ax1.set_theta_zero_location("N")
            ax1.set_theta_direction(-1)
            #plt.plot(azis[(locs == '00') & (p0s == per) ], HVs[(locs == '00') & (p0s == per)],'.', markersize=12, alpha=.1, color='blue', label = '00 Measurement')


            plt.ylim((-0.6,0.6))
            plt.yticks([-0.3,0.0, 0.3])
            angs=[]
            m00=[]
            m10=[]
            std00=[]
            std10=[]
            for ang in range(0,360, 30):
                ang *= np.pi/180.
                cu00 = HVs[(locs == '00') & (p0s == per)]
                azis00 = azis[(locs == '00') & (p0s == per) ]
                cu00 = cu00[(azis00 > ang) & ( azis00 <= (ang+ 30.*np.pi/180.))]
                #print('Here is azis00')
                #print(azis00)
                cu10 = HVs[(locs == '10') & (p0s == per)]
                azis10 = azis[(locs == '10') & (p0s == per) ]
                cu10 = cu10[(azis10 > ang) & ( azis10 <= (ang+ 30.*np.pi/180.))]
                #print('HEre is cu00')
                #print(cu00)
                if((~np.isnan(np.mean(cu00))) & (~np.isnan(np.mean(cu10)))):
                    m00.append(np.mean(cu00))
                    m10.append(np.mean(cu10))
                    std00.append(np.std(cu00))
                    std10.append(np.std(cu10))
                    angs.append(float(ang))
            
            angs.append(angs[0])
            m00.append(m00[0])
            m10.append(m10[0])
            #print(angs)
            #print(m00)
            angs = np.asarray(angs)
            std00 = np.asarray(std00)
            std10=np.asarray(std10)
            m00 = np.asarray(m00)
            m10 = np.asarray(m10)
            std00 = np.nan_to_num(std00)
            std10 = np.nan_to_num(std10)
            print(m00)
            print(std00)

            


            points = np.concatenate([angs,m00]).reshape((2,len(m00))).T
            #hull = Delaunay(points)
            #plt.subplot(211)
            plt.polar(angs, m00 , markersize=18, color='C' + str(idx), linewidth=3., label='00: ' + str(per) + ' s', alpha=.5)

            #plt.subplot(212)
            plt.polar(angs, m10,  markersize=18, color='C' + str(idx), label='10: ' + str(per) + ' s', linewidth=3., alpha=.5, linestyle=":")
            plt.xticks([0., 90.*np.pi/180., 180.*np.pi/180., 270*np.pi/180.], ['N', 'E', 'S', 'W'])
            plt.ylim((-0.6,0.6))
            plt.yticks([-0.3,0.0, 0.3])
            #plt.text(10., 0.3, letts[idx] + ' ' + str(per) + ' s')
            #if idx < 4:
            #    plt.xticks([])
            #if idx == 0:
                
            #    plt.title(sta + ' H/V Ratio and Azimuth')
        #plt.subplot(515)
        #plt.subplot(212)
        plt.legend(ncol=5, loc=8, bbox_to_anchor=(0.5, -0.15))
        #plt.xlabel('Back-Azimuth (degrees)')
        #plt.tight_layout()
        plt.title('H/V Ratios by Azimuth ' + sta)
        plt.savefig( sta + '_ROSE.jpg',format='JPEG', dpi=400)
        #plt.show()
        plt.clf()
        plt.close()
