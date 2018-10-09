#!/usr/bin/env python

import glob
import numpy as np
import math
from scipy.optimize import nnls, lsq_linear

def computeA(sta, freq, loc, quad, debug = True):
    files = glob.glob('Results/*' + sta + '*Results.csv')
    azis = []
    amps = []
    hvs = []
    EV =[]
    NV =[]
    for curfile in files:

        with open(curfile,'r') as f:
            for line in f:
                if (loc in line.split(',')[1]) and (freq == float(line.split(',')[4])):
                    data = line
                    data = data.split(',')
                    if math.isnan(float(data[-2])):
                       continue
                    if math.isnan(float(data[-1])):
                        continue
                    if math.isnan(float(data[11])):
                        continue
                        
                    

                    azis.append(float(data[6]))
                    NV.append(float(data[-2]))
                    EV.append(float(data[-1]))
                    hvs.append(float(data[-3]))
    #print(hvs)
    HV = np.mean(hvs) 
    
    
    
    
    
    
    
    
    
    
    #print('Here is HV:' + str(HV))               
    for idx, pair in enumerate(zip(azis, NV, EV, hvs)):
        theta1 = np.deg2rad(pair[0])
        
        theta2 = np.deg2rad(pair[0]+ 90.)
        if pair[0] >= quad[0] and pair[0] <= quad[1]:
            #ampN = pair[1]/pair[3]
            ampN = pair[3]/(pair[1]*np.cos(theta1))
            ampE = pair[3]/(pair[2]*np.cos(theta2))
            #ampE = pair[2]/pair[3]
            #angs1 = np.asarray([[np.cos(theta1)**2, np.sin(theta1)**2, np.sin(2*theta1)]])
            #angs2 = np.asarray([[np.cos(theta2)**2, np.sin(theta2)**2, np.sin(2*theta2)]])
            angs1 = np.asarray([[np.cos(theta1), np.sin(theta1)]])
            angs2 = np.asarray([[np.cos(theta2), np.sin(theta2)]])
            if 'A1' not in vars():
                
                A1=angs1
                A2 = angs2
                x1 = np.asarray([[ampN]])
                x2 = np.asarray([[ampE]])
                x1G = np.asarray([ampN])
                x2G = np.asarray([ampE])
            else:
                A1=np.append(A1,angs1, axis=0)
                A2=np.append(A2,angs2, axis=0)
                x1=np.append(x1,[[ampN]],axis=0)
                x2 = np.append(x2,[[ampE]],axis=0)
                x1G=np.append(x1G,[ampN],axis=0)
                x2G = np.append(x2G,[ampE],axis=0)
    
    sol1 = lsq_linear(A1,x1G, bounds=(0.,1.))
    sol1 = sol1.x
    print(sol1)
    sol2 =lsq_linear(A2,x2G, bounds = (0,1.))
    sol2 = sol2.x
    print('Here is sol:' + str(sol1[0]))
    pinv1 = np.linalg.pinv(A1)
    corrs1 = np.matmul(pinv1,x1)
    resi1 = np.matmul(A1,corrs1) - x1
    neve = len(resi1)
    resi1 = np.std(resi1)
    resi1 = (np.var(np.matmul(A1,corrs1)-x1)/np.var(x1))
    resi1G = 100.*(1.-(np.var(np.matmul(A1,sol1)-x1G)/np.var(x1G)))
    print('Here is resi1G' + str(resi1G))
    pinv2 = np.linalg.pinv(A2)
    corrs2 = np.matmul(pinv2,x2)
    resi2 = ((np.var(np.matmul(A2,corrs2)-x2)/np.var(x2)))
    
    print('Here is vari1: ' + str(np.var(x1G)))
    print('Here is reduced: '+ str((np.var(np.matmul(A2,corrs2)-x2))))
    
    
    resi2G = 100.*(1.-(np.var(np.matmul(A2,sol2)-x2G)/np.var(x2G)))
    del A1
    del A2
    return sol1, sol2, neve, resi1G, resi2G
    
stas = glob.glob('Results/*.csv')

goodstas = []
for sta in stas:
    goodstas.append(sta.split('_')[1])
goodstas = list(set(goodstas))

print(goodstas)

#goodstas = ['HIA']

f2 = open('Goodstuff_NEW','w')
for sta in goodstas:
    print(sta)
    for loc in ['00','10']:
        for freq in [0.00666666666667, 0.01, 0.0133333333333, 0.02, 0.04]:
            try:
            #if True:
                A1, A2, ln, resi1, resi2 = computeA(sta, freq, loc, [0., 360.])
                f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad 1, ' + str(ln) + ', A1, ' + str(A1[0]) +  ', ' + str(A1[1]) + ', ' + str(resi1) + ' \n')
                f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad 1, ' + str(ln) + ', A2, ' + str(A2[0]) +  ', ' + str(A2[1]) + ', ' + str(resi2) + ' \n')
                #f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad 1, ' + str(ln) + ', A1, ' + str(A1[0][0]) +  ', ' + str(A1[1][0]) + ', ' + str(resi1) + ' \n')
                #f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad 1, ' + str(ln) + ', A2, ' + str(A2[0][0]) +  ', ' + str(A2[1][0]) + ', ' + str(resi2) + ' \n')
                #A1, A2, ln = computeA(sta, freq, loc, [90., 180.])
                #f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad 2, ' + str(ln) + ', A1, ' + str(A1[0][0]) +  ', ' + str(A1[1][0]) + ', ' + str(A1[2][0]) + ' \n')
                #f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad 2, ' + str(ln) + ', A2, ' + str(A2[0][0]) +  ', ' + str(A2[1][0]) + ', ' + str(A2[2][0]) + ' \n')
                #A1, A2, ln = computeA(sta, freq, loc, [180., 270.])
                #f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad 3, ' + str(ln) + ', A1, ' + str(A1[0][0]) +  ', ' + str(A1[1][0]) + ', ' + str(A1[2][0]) + ' \n')
                #f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad 3, ' + str(ln) + ', A2, ' + str(A2[0][0]) +  ', ' + str(A2[1][0]) + ', ' + str(A2[2][0]) + ' \n')
                #A1, A2, ln = computeA(sta, freq, loc, [270., 360.])
                #f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad 4, ' + str(ln) + ', A1, ' + str(A1[0][0]) +  ', ' + str(A1[1][0]) + ', ' + str(A1[2][0]) + ' \n')
                #f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad 4, ' + str(ln) + ', A2, ' + str(A2[0][0]) +  ', ' + str(A2[1][0]) + ', ' + str(A2[2][0]) + ' \n')
            except:
                print('Can not compute: ' + sta + ' ' + loc)
f2.close()
