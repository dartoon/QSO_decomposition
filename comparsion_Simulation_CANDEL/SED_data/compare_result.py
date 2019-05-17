#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:00:30 2019

@author: Dartoon

Comparing the fitting between Xuheng and Federica
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../../py_tools')
from load_result import load_host_p, load_err

ID = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242',
      'CID3570', 'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607',
      'CID70', 'LID1273', 'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396',
      'CDFS-1', 'CDFS-229', 'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136',
      'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']

Mstar = load_host_p(ID=ID, folder='../../')[1]
Mstar_err = load_err(prop = 'Mstar', ID=ID)
LR = load_host_p(ID=ID, folder='../../', dm = 0)[0] #!!! This dm is important 
LR_err = load_err(prop = 'LR', ID=ID)
Fede = np.loadtxt('Summary.txt') #0ID 1M*_SED 2M*_IMAGEDEC 3LR_SED 4LR_IMAGEDEC 5agreement

bool = [Fede[:,3]!=-99]  #exclude CID255 at this moment
#%%
plt.figure(figsize=(10, 10))
x = np.linspace(8., 12, 20)
y = x
plt.plot(x,y, 'gray', alpha=0.5)
#plt.plot(LR[bool], Fede[:,3][bool], 'bo', label='SED only')
plt.errorbar(LR[bool], Fede[:,3][bool], xerr=[np.abs(LR_err)[:,0][bool], np.abs(LR_err)[:,1][bool]],yerr=0.2 + np.zeros(len(Mstar[bool])),fmt='.',color='blue',markersize=15, label='SED only')

#plt.plot(LR[bool], Fede[:,4][bool], 'r^', label='fix HST result')
plt.xlim([8.8,11.8])
plt.ylim([8.8,11.8])
plt.title("Comparsion of LR",fontsize=35)
plt.xlabel("Xuheng log$(L_R/L_{\odot})$",fontsize=35)
plt.ylabel("Federica log$(L_R/L_{\odot})$",fontsize=35)
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)
#plt.legend(prop={'size':20})
plt.show()

#%%
plt.figure(figsize=(10, 10))
x = np.linspace(8.5, 12.5, 20)
y = x
plt.plot(x,y, 'gray',  alpha=0.5)
plt.errorbar(Mstar[bool], Fede[:,1][bool], xerr=[np.abs(Mstar_err)[:,0][bool], np.abs(Mstar_err)[:,1][bool]],yerr=0.3 + np.zeros(len(Mstar[bool])),fmt='.',color='blue',markersize=15, label='SED only')
#plt.plot(Mstar[bool], Fede[:,2][bool], 'r^', label='fix HST result')
plt.xlim([8.5,12.5])
plt.ylim([8.5,12.5])
plt.title("Comparsion of M*",fontsize=35)
plt.xlabel("Xuheng log$(M_*/M_{\odot})$",fontsize=35)
plt.ylabel("Federica log$(M_*/M_{\odot})$", fontsize=35)
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)
#plt.legend(prop={'size':20})
plt.show()