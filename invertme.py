#!/usr/bin/env python

import glob
import numpy as np
import math
from scipy.optimize import lsq_linear
import matplotlib.pyplot as plt

def computeA(sta, freq, loc, quad, debug = True):
    files = glob.glob('Results/*' + sta + '*Results.csv')
    azis = []
    HV = []
    EV =[]
    NV =[]
    for curfile in files:
        with open(curfile,'r') as f:
            for line in f:
                if (loc in line.split(',')[1]) and (freq == float(line.split(',')[4])):
                    data = line
                    data = data.split(',')
                    if float(data[-6]) > 0.:
                        azis.append(float(data[6]))
                        NV.append(float(data[-2]))
                        EV.append(float(data[-1]))
                        HV.append(float(data[-6]))
    
    
    
    
    for comp in ['N', 'E']:
        if comp == 'N':
            CV = NV
        else:
            CV = EV
        for idx, pair in enumerate(zip(azis, CV, HV)):


            if comp == 'N':
                theta = np.deg2rad(pair[0])
            else:
                theta = np.deg2rad(pair[0]+ 90.)
            
            if (pair[0] >= quad[0]) and (pair[0] <= quad[1]):
                #print('Here is azi: ' + str(pair[0]))
                if pair[2] == 0.:
                    continue
                amp = pair[1]/pair[2]
                angs = np.asarray([[np.cos(theta)**2, np.sin(theta)**2, np.sin(2*theta)]])
                if 'A' not in vars():   
        del A
    return corrs1, corrs2, ln, reduction, reduction2 
    
stas = glob.glob('Results/*.csv')

goodstas = []
for sta in stas:
    goodstas.append(sta.split('_')[1])
goodstas = list(set(goodstas))

print(goodstas)

#goodstas = ['HIA']

f2 = open('Goodstuff','w')
for sta in goodstas:
    print(sta)
    for loc in ['00','10']:
        for freq in [0.00666666666667, 0.01, 0.0133333333333, 0.02, 0.04]:

                for idx, quad in enumerate([[0.,360.]]):
                    #try:
                    if True:
                        A1, A2, ln, reduction, reduction2 = computeA(sta, freq, loc, quad)
                     
                        if ln > 5:
                            f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad ' + str(idx + 1) + ', ' + str(ln) + ', A1, ' + str(round(A1[0],2)) +  ', ' + str(round(A1[1],2)) + ', ' + str(round(A1[2],2)) + ', ' + str(round(reduction,2)) + ' \n')
                            f2.write(sta + ', ' + loc + ', ' + str(freq) + ', Quad ' + str(idx + 1) + ', ' + str(ln) + ', A2, ' + str(round(A2[0],2)) +  ', ' + str(round(A2[1],2)) + ', ' + str(round(A2[2],2)) + ', ' + str(round(reduction2,2)) + ' \n')
                    #except:
                    #    print('Have a problem on: ' + str(idx))
                    #    continue

f2.close()
