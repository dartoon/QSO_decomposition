#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 18:40:15 2019

@author: Dartoon

Measure the FWHM of the QSO with host fitting Reff<0.2
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'../py_tools')
from load_result import load_n, load_re
from flux_profile import SB_profile
from adjustText import adjust_text

pix_s = 0.0642
IDs_all = np.array(['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255'])

Re_results = np.asarray(load_re(IDs_all))
IDs = IDs_all[Re_results[:,0]<0.2]

Re_s = (Re_results[Re_results[:,0]<0.2])

FWHM_list = []
for i in range(len(IDs)):
    data= pyfits.getdata('{0}/analysis/{0}_cutout.fits'.format(IDs[i]))
    center = [len(data)/2,len(data)/2]
    r_SB, r_grids = SB_profile(data, center,ifplot=False, radius=35, start_p=0.5, grids=50, gridspace='log')
    idx = np.sum([r_SB > r_SB[0]/2])-1
    FWHM = r_SB[idx]* pix_s
    FWHM_list.append(FWHM)
    
plt.figure(figsize=(10, 10))
plt.errorbar(Re_s[:,0], FWHM_list,xerr= Re_s[:,1], color='red', fmt='o',ecolor='gray' )
texts = []
for i in range(len(IDs)):
    texts.append(plt.text(Re_s[:,0][i], FWHM_list[i], IDs[i], fontsize=17))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))
x=np.linspace(0,2,15)
y = x
plt.plot(x,y,'b')
plt.title('Compare host_Reff to QSO FWHM', fontsize=30)
plt.xlabel('Reff (arcsec)',fontsize=25)
plt.ylabel('FWHM (arcsec)', fontsize=25)
plt.xlim(0,1.3)
plt.ylim(0,1.3)
plt.tick_params(labelsize=25)     
#plt.savefig('/Users/Dartoon/Desktop/Comp_gtemp_{0}.pdf'.format(galay_temp))
plt.show()

#%%
data= pyfits.getdata('{0}/analysis/{0}_cutout.fits'.format('CID50'))
center = [len(data)/2,len(data)/2]
r_SB, r_grids = SB_profile(data, center,ifplot=False, radius=35, 
                           start_p=0.5, grids=50, gridspace='log', 
                           fits_plot=True)
idx = np.sum([r_SB > r_SB[0]/2])-1
print "FWHM", r_SB[idx]* pix_s
from matplotlib.ticker import AutoMinorLocator
minorLocator = AutoMinorLocator()
fig, ax = plt.subplots()
plt.plot(r_grids, r_SB, 'x-')
plt.scatter(r_grids[idx], r_SB[idx], c='k', marker="s", edgecolors='k',zorder=100)
ax.xaxis.set_minor_locator(minorLocator)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4, color='r')
plt.grid()
ax.set_ylabel("Surface Brightness")
ax.set_xlabel("Pixels")
ax.set_xscale('log')
plt.xlim(0.5*0.7,) 
plt.grid(which="minor")
plt.show()

