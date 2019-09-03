#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 08:58:27 2019

@author: Dartoon

Calculate the host flux ratio
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../../../py_tools')
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'

ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597', 'CID1281','CID255']
from load_result import load_host_ratio
WFC3_ratio = np.array(load_host_ratio(ID, folder='../../../'))

siml_mags_nicola = np.loadtxt('lagn_lhost_Vband.dat')
host_f = 10.**(-0.4*(siml_mags_nicola[:,0]))
agn_f = 10.**(-0.4*(siml_mags_nicola[:,1]))
Nicola_ratio = host_f/(agn_f+host_f)

Aklant_ratio = np.loadtxt('../../Aklant/host_ratio/host_total_ratio.txt')


plt.figure(figsize=(10,8))
high1, x1, _ = plt.hist(WFC3_ratio[:,0],normed=True, histtype=u'step',
         label=('Observed sample'), linewidth = 2, color='green', zorder = 100)
high2, x2, _ = plt.hist(Aklant_ratio*100.,normed=True, histtype=u'step',
         label=('MBII simulation'), linewidth = 2, color='steelblue')
high0, x0, _ = plt.hist(Nicola_ratio*100.,normed=True, histtype=u'step',
         label=('SAM model'), linewidth = 2, color='orange')

x0_m = np.median(Nicola_ratio*100.)
high_m0 = high0[np.where(abs(x0_m-x0) == abs(x0_m-x0).min())[0][0]-1]
x1_m = np.median(WFC3_ratio[:,0])
high_m1 = high1[np.where(abs(x1_m-x1) == abs(x1_m-x1).min())[0][0]-1]
x2_m = np.median(Aklant_ratio*100.)
high_m2 = high2[np.where(abs(x2_m-x2) == abs(x2_m-x2).min())[0][0]]


plt.plot(np.linspace(0,high_m0)*0+np.median(x0_m) , np.linspace(0,high_m0),'--' , linewidth = 2,color='orange')
plt.plot(np.linspace(0,high_m1)*0+np.median(x1_m) , np.linspace(0,high_m1),'--' , linewidth = 2, color='green')
plt.plot(np.linspace(0,high_m2)*0+np.median(x2_m) , np.linspace(0,high_m2),'--' , linewidth = 2, color='steelblue')

plt.text(np.median(x0_m)-0.2, high_m0*0.25, '{0}%'.format(round(np.median(x0_m),1)), color='orange',fontsize=25)
plt.text(np.median(x1_m)-0.2, high_m1*1.03, '{0}%'.format(round(np.median(x1_m),1)), color='green',fontsize=25)
plt.text(np.median(x2_m)-12, high_m2*0.35, '{0}%'.format(round(np.median(x2_m),1)), color='steelblue',fontsize=25)

plt.xlabel("Host flux ratio (%)",fontsize=27)
plt.ylabel("Number Density",fontsize=27)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})

plt.ylim(0,0.033)
plt.yticks([])
#plt.savefig('comp_host_ratio.pdf')
plt.show()

from scipy import stats
pvalue_0 = stats.ks_2samp(WFC3_ratio[:,0], Aklant_ratio*100).pvalue
pvalue_1 = stats.ks_2samp(WFC3_ratio[:,0], Nicola_ratio*100).pvalue
pvalue_2 = stats.ks_2samp(Aklant_ratio, Nicola_ratio).pvalue
print "p-value:", round(pvalue_0,3), round(pvalue_1,3), round(pvalue_2,3)
