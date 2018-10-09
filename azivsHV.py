#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np
debug = True



stas = ['TUC']

def readHVfile(filename):
    results =[]
    with open(filename,'r') as f:
        f.readline()
        for line in f:
            dic = {}
            dic['f0'] = float(line.split(',')[4])
            dic['loc'] = (line.split(',')[1]).replace(' ','')
            dic['HV'] = float(line.split(',')[8])
            results.append(dic)
    with open(filename.replace('NEWRESULTS','AZIS'),'r') as f:    
        f0s =[]
        locs = []
        azisR = []
        azi =[]
        amps = []
        f.readline()
        for line in f:
            f0s.append(float(line.split(',')[4]))
            locs.append((line.split(',')[1]).replace(' ', ''))
            azisR.append(float(line.split(',')[8]))
            amps.append(float(line.split(',')[9]))
            azi.append(float(line.split(',')[6]))
    f0s = np.asarray(f0s)
    locs = np.asarray(locs)
    amps = np.asarray(amps)
    azisR = np.asarray(azisR)
    azi = np.asarray(azi)
    # for each dictionary in the results find the azimuth and add it
    for dic in results:
        
        if dic['loc'] == '00':
            dic['aziR'] = azisR[('10' == locs) & (dic['f0'] == f0s)][0]
            dic['ampR'] = amps[('10' == locs) & (dic['f0'] == f0s)][0]
            dic['azi'] = azi[('10' == locs) & (dic['f0'] == f0s)][0]
        else:
            print(dic['loc'])
            print((dic['loc'] == locs) & (dic['f0'] == f0s))
            dic['aziR'] = azisR[(dic['loc'] == locs) & (dic['f0'] == f0s)][0]
            dic['ampR'] = amps[('10' == locs) & (dic['f0'] == f0s)][0]
            dic['azi'] = azi[('10' == locs) & (dic['f0'] == f0s)][0]
    return results






pair = ['azi','HV','ampR']



for sta in stas:
    HVSfiles = glob.glob('NEWRESULTS/*' + sta + '*')
    estimate =[]
    for HVfile in HVSfiles:
        try:
        #if True:
            estimate += readHVfile(HVfile)
        except:
            continue
    #print(estimate)
    
    fig = plt.figure(1)
    locs = list(set([ dic['loc'] for dic in estimate]))
    print(locs)
    cols = {}
    for idx, loc in enumerate(locs):
        cols['loc'] = 'C' + str(idx)
        first =[]
        second =[]
        third = []
        for dic in estimate:
            if (dic['f0'] == 0.01) & (dic['loc'] == loc):
                first.append(dic[pair[0]])
                second.append(dic[pair[1]])
                third.append(dic[pair[2]])
        plt.plot(first, second, '.', c=cols['loc'])
        plt.plot(first, np.log10(np.asarray(third)), 'o', c=cols['loc'])
    plt.show()

