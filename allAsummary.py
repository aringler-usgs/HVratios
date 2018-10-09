#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys

stas = glob.glob('Results/*.csv')
stas = [sta.split('_')[1] for sta in stas]
stas = list(set(stas))
stas.sort()
#stas =['ANMO','HIA']
#print(stas)

import matplotlib as mpl
#mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=12)

stas = glob.glob('Results/*.csv')
stas = [sta.split('_')[1] for sta in stas]
stas = list(set(stas))
stas.sort()
    


letts = ['(a)', '(b)', '(c)', '(d)', '(e)']
fig =plt.figure(1, figsize=(12,18))
plt.subplots_adjust(wspace=0.001)
for idx, freq in enumerate([0.00666666666667, 1./100., 0.0133333333333, 1./50., 1./25]):
    As00=[]
    As10=[]
    stas00=[]
    stas10=[]
    As00E=[]
    As10E=[]
    stas00E=[]
    stas10E=[]
    f = open('Goodstuff_NEW','r')
    for line in f:
        
        line = line.replace(' ','')
        line = line.split(',')
        if float(line[2]) == freq:
            if line[1] == '00':
                if line[5] == 'A1':
                    stas00.append(line[0])
                    val = np.sqrt(float(line[6])**2 + float(line[7])**2 + float(line[8])**2)
                    if val >= .12:
                        val = .12
                    As00.append(val)
                else:
                    stas00E.append(line[0])
                    val = np.sqrt(float(line[6])**2 + float(line[7])**2 + float(line[8])**2)
                    if val >= 0.12:
                        val = .12
                    As00E.append(val)
            if line[1] == '10':
                if line[5] == 'A2':
                    stas10.append(line[0])
                    val = np.sqrt(float(line[6])**2 + float(line[7])**2 + float(line[8])**2)
                    if val >= 0.12:
                        val = 0.12
                    As10.append(val)
                else:
                
                    stas10E.append(line[0])
                    val = np.sqrt(float(line[6])**2 + float(line[7])**2 + float(line[8])**2)
                    if val >= 0.12:
                        val = 0.12
                    As10E.append(val)
    f.close()
                

    plt.subplot(151+idx)
    
    for hidx, ele in enumerate(stas):
        try:
            aidx00 = stas00.index(ele)
            aidx00E = stas00E.index(ele)
            if hidx == 0:
                plt.plot([As00[aidx00]], [hidx], label='00 North-South', color='C0', marker='.', linestyle='None')
                plt.plot([As00E[aidx00E]], [hidx], label='00 East-West', color='C0', marker='v', linestyle='None')
            else:
                plt.plot([As00[aidx00]], [hidx], color='C0', marker='.', linestyle='None')
                plt.plot([As00E[aidx00E]], [hidx],  color='C0', marker='v', linestyle='None')
        except:
            pass
        try:
            aidx10 = stas10.index(ele)
            aidx10E = stas10E.index(ele)
            if hidx == 0:
                plt.plot([As10[aidx10] ] , [hidx], label='10 North-South', alpha=.5, color='C1', marker='.', linestyle='None')
                plt.plot([As10E[aidx10E] ] , [hidx], label='10 East-West', alpha=.5, color='C1', marker='v', linestyle='None')
            else:
                plt.plot([As10[aidx10] ] , [hidx],  alpha=.5, color='C1', marker='.', linestyle='None')
                plt.plot([As10E[aidx10E] ] , [hidx], alpha=.5, color='C1', marker='v', linestyle='None')
        except:
            pass
        #plt.plot(HV00min,range(len(stas)),'.',label='00')
        #plt.plot(HV10max,range(len(stas)),'.',label='10')
        #plt.plot(HV10min,range(len(stas)),'.',label='10')
    plt.text(0.01, len(stas)+.5, letts[idx] + ' ' + str(int(round(1/freq,0))) + ' s')
    plt.ylim((-0.5, len(stas)+2.5))
    plt.xlim((0.,0.12))
    plt.xticks([0.03, 0.06, 0.09])
    if idx == 0:
        plt.yticks(range(len(stas)), stas)
    else:
        plt.yticks([])
    if idx == 2:
        plt.xlabel('$\sqrt{A_{i11} ^2 + A_{i22} ^2 + A_{i12} ^2}$')
    #if idx == 4:
        #cb=plt.colorbar(sc, ticks=[0.5, 1., 1.5])
        #cb.set_label('H/V Ratio')

plt.subplot(153)
plt.legend(ncol=2, loc='center', bbox_to_anchor=(0.5, -.1))
#plt.show()

plt.savefig('allatermsNEW.pdf',format='PDF',dpi=400)
