#!/usr/bin/env python

import glob
import numpy as np
import matplotlib.pyplot as plt


def getHV(sta, loc, freq):
    files = glob.glob('Results/*' + sta + '*Results.csv')
    azis = []
    HV = []
    for curfile in files:
        with open(curfile,'r') as f:
            for line in f:
                if (loc in line.split(',')[1]) and (freq == float(line.split(',')[4])):
                    data = line
                    data = data.split(',')
                    if float(data[7]) > .8:
                        azis.append(float(data[6]))
                        HV.append(float(data[9]))
    HV = np.asarray(HV)
    azis = np.asarray(azis)
    return HV, azis




stas = glob.glob('Results/*.csv')

goodstas = []
for sta in stas:
    goodstas.append(sta.split('_')[1])
goodstas = list(set(goodstas))



for sta in goodstas:
    fig = plt.figure(1)
    for loc in ['00','10']:
        for freq in [0.00666666666667, 0.01, 0.0133333333333, 0.02, 0.04]:
            HV, azis = getHV(sta, loc, freq)
            azis = azis[(HV <= 2.*np.std(HV))]
            HV = HV[(HV <= 2.*np.std(HV))]
            print('Here is HV: ' + str(np.mean(HV)) + ' +/- ' + str(np.std(HV)))
            plt.plot(azis, HV, '.', label = sta + ' ' + loc + ' ' + str(int(round(1./freq,0))))
    plt.legend()
    plt.xlabel('Azimuth (deg)')
    plt.ylabel('HV Ratio')
    plt.show()
