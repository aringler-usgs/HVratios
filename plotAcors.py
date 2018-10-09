#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import glob
import matplotlib as mpl
#mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=14)





stas = glob.glob('Results/*.csv')

goodstas = []
for sta in stas:
    goodstas.append(sta.split('_')[1])
goodstas = list(set(goodstas))

#goodstas =['PAYG']

for sta in goodstas:

    rans = np.asarray([25., 50., 75., 100., 150.])
    f = open('Goodstuff_NEW','r')
    fig = plt.figure(1,figsize=(12,12))
    plt.subplots_adjust(hspace=0.001)
    for line in f:
        if sta not in line:
            continue
        line = line.split(',')
        idx = np.argmin(np.abs(1./float(line[2])- rans))
        if line[1].replace(' ','') == '00':
            c='C0'
        else:
            c='C1'
        if line[5].replace(' ','') == 'A1':
            symbol = 'o'
        else:
            symbol ='v'
        plt.subplot(311)
        plt.title('Local Corrections ' + sta)
        if idx == 0:
            if line[5].replace(' ','') == 'A1':
                lbstr = line[1].replace(' ','') + ' N/S'
            else:
                lbstr = line[1].replace(' ','') + ' E/W'
                
            plt.plot(idx,float(line[6]),symbol, color=c, markersize=17, label=lbstr, alpha=0.5)
        else:
            plt.plot(idx,float(line[6]),symbol, color=c, markersize=17, alpha=0.5)
        plt.ylim((-0.2, 0.05))
        plt.yticks([-0.15, -0.1 ,-0.05, 0.0])
        plt.xlim([-0.25, 4.5])
        plt.text(0., 0.025, '(a) $A_{i11}$')
        plt.xticks([])
        plt.legend(ncol=4, loc=8)
        plt.subplot(312)
        plt.plot(idx,float(line[7]),symbol,color=c, markersize=17, alpha=0.5) 
        #plt.yticks([ 1., 2.])
        plt.text(-0.25, 0.025, '(a) $A_{i22}$')

        plt.ylim((-0.2, 0.05))
        plt.yticks([-0.15, -0.1 ,-0.05, 0.0])
        plt.ylabel('Coefficient Value') 
        plt.xlim([-0.5, 4.5])
        plt.xticks([])

        plt.subplot(313)
        
        plt.plot(idx,float(line[8]),symbol,color=c, markersize=17, alpha=0.5)
        plt.text(-.25, 0.025, '(a) $A_{i12}$')
        #plt.yticks((-2., -1., 0., 1., 2.))
        plt.xlim([-0.5, 4.5])
        plt.ylim((-0.2, 0.05))
        plt.yticks([-0.15, -0.1 ,-0.05, 0.0])

    plt.xlabel('Period (s)')
    plt.xticks([0., 1., 2., 3., 4.], ['25', '50', '75', '100', '150']) 
    #plt.xticks(rans)
    #plt.show()
    
    plt.savefig(sta + '_coeff.pdf', format='PDF',dpi =400)
    plt.clf()
    f.close()


