#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 14:02:54 2019

@author: Dartoon

Plot the two figures that Tommaso asked.
1.Color Vs z
2.Color Vs. M/L
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
import matplotlib as mpl
mpl.rc('image', cmap='jet')

import glob
import sys
sys.path.insert(0,'../../py_tools')
from dmag import pass_dmag
from filter_info import filt_info
from load_result import load_zs, load_mag, load_re, load_n, load_flux

ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']

filt = [filt_info[ID[i]] for i in range(len(ID))]

zs = np.asarray(load_zs(ID))
from load_result import load_host_p

plt.figure(figsize=(15, 11))
host_mag_WFC3 = np.array(load_mag(ID, folder = '../../', flt = 'WFC3'))[0]
host_mag_ACS = []
for i in range(len(ID)):
    ifexit = glob.glob('../../analysis_ACS/{0}'.format(ID[i]))
    if ifexit!= []:
        host_mag_ACS.append(load_mag([ID[i]], folder = '../../', flt = 'ACS')[0][0])
    else:
        host_mag_ACS.append(-99)
host_mag_ACS = np.asarray(host_mag_ACS)
two_bands = [host_mag_ACS>0]

R_lumi, Mstar, _ = load_host_p(ID, dm = 0, folder = '../../') # This dm is important
zs_c, host_mag_ACS_c, host_mag_WFC3_c, R_lumi_c, Mstar_c  = \
zs[two_bands],  host_mag_ACS[two_bands], host_mag_WFC3[two_bands], R_lumi[two_bands], Mstar[two_bands]

z_break = 1.44
plt.scatter(zs_c[zs_c>1.44], (host_mag_ACS_c - host_mag_WFC3_c)[zs_c>1.44],s=680, c ='green',
            marker="*",zorder=100, vmin=0, vmax=7, edgecolors='white', label='WFC/F140W')
plt.scatter(zs_c[zs_c<1.44], (host_mag_ACS_c - host_mag_WFC3_c)[zs_c<1.44],s=280, c ='green',
            marker="o",zorder=100, vmin=0, vmax=7, edgecolors='white', label='WFC/F125W')
plt.plot(np.arange(1,3.4,0.2) * 0 + z_break, np.arange(1,3.4,0.2))
plt.ylabel("mag_$_{ACS}$ - mag_$_{WFC3}$",fontsize=35)
plt.xlabel("redshift",fontsize=35)
plt.tick_params(labelsize=25)
plt.ylim([1.2, 2.8])
plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':24},ncol=1)
plt.show()

#Mstar = load_host_p(ID, folder = '../../')[1]

plt.figure(figsize=(15, 11))
plt.scatter((10**Mstar_c/10**R_lumi_c)[zs_c>1.44], (host_mag_ACS_c - host_mag_WFC3_c)[zs_c>1.44],s=680, c ='green',
            marker="*",zorder=100, vmin=0, vmax=7, edgecolors='white', label='WFC/F140W')
plt.scatter((10**Mstar_c/10**R_lumi_c)[zs_c<1.44], (host_mag_ACS_c - host_mag_WFC3_c)[[zs_c<1.44]],s=280, c ='green',
            marker="o",zorder=100, vmin=0, vmax=7, edgecolors='white', label='WFC/F125W')
plt.ylabel("mag_$_{ACS}$ - mag_$_{WFC3}$",fontsize=35)
plt.xlabel("(M_*/M_sun)/(L_R/L_Rsun)",fontsize=35)
plt.tick_params(labelsize=25)
plt.ylim([1.2, 2.8])
plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':24},ncol=1)
plt.show()