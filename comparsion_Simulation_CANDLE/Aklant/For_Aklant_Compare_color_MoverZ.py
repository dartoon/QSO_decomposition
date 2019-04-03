#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 14:02:54 2019

@author: Xuheng

Plot the two figures that Tommaso asked.
1.Color Vs z
2.Color Vs. M/L
"""
import numpy as np
import matplotlib.pyplot as plt

zs = np.array([1.301, 1.667, 1.552, 1.567, 1.618, 1.532, 1.244, 1.407, 1.478,
       1.239, 1.294, 1.617, 1.527, 1.579, 1.551, 1.516, 1.6  , 1.483,
       1.272, 1.445])  # The redshift of the targerts, z > 1.44 observed by F140w, z < 1.44 by F125w,

host_mag_WFC3 = np.array([21.994280325991404, 21.862018827793285, 21.48186850807134,
       21.510426015577707, 21.279055086951168, 21.160935881370094,
       21.164294483113114, 21.17569181223492, 21.195529570494177,
       20.93361752085037, 21.18894227543559, 20.938781271397527,
       21.25086447413113, 21.454856873262585, 21.86934793003787,
       21.160337845783772, 21.404246870835454, 21.821946001190145,
       21.868092263533335, 22.87919564399801]) # The host galaxy magnitudes by WFC3/F140w or F125w, depends on the redshift.

host_mag_ACS = np.array([23.76505005, 24.62530304, 23.213549  , 23.44615606, 23.7181948 ,
       23.5953285 , 22.97053958, 22.73207916, 23.35202394, 22.5043342 ,
       23.56621055, 23.28922372, 23.08836815, 23.25392999, 23.89780312,
       22.59159894, 23.36364747, 23.67104111, 23.56335223, 24.82799557]) # The host galaxy magnitudes by WFC3/814W. The -99 are the 

plt.figure(figsize=(15, 11))

z_break = 1.44
plt.scatter(zs[zs>1.44], (host_mag_ACS - host_mag_WFC3)[zs>1.44],s=680, c ='green',
            marker="*",zorder=100, vmin=0, vmax=7, edgecolors='white', label='WFC/F140W')
plt.scatter(zs[zs<1.44], (host_mag_ACS - host_mag_WFC3)[zs<1.44],s=280, c ='green',
            marker="o",zorder=100, vmin=0, vmax=7, edgecolors='white', label='WFC/F125W')
plt.plot(np.arange(1,3.4,0.2) * 0 + z_break, np.arange(1,3.4,0.2))
plt.ylabel("mag_$_{ACS}$ - mag_$_{WFC3}$",fontsize=35)
plt.xlabel("redshift",fontsize=35)
plt.tick_params(labelsize=25)
plt.ylim([1.2, 2.8])
plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':24},ncol=1)
plt.show()

R_lumi = np.array([10.70304515, 10.98416944, 11.05007582, 11.04991576, 11.18053661,
       11.16344494, 10.98195418, 11.1263415 , 11.10772781, 11.06953281,
       11.01874675, 11.31589169, 11.12369835, 11.0810814 , 10.89432957,
       11.15156228, 11.11681513, 10.86116271, 10.72672663, 10.40771812]) # The R band Luminosity in AB system, in log10(L_R/L_Rsun) unit.

Mstar = np.array([10.43450801, 10.7156323 , 10.78153868, 10.78137863, 10.91199948,
       10.8949078 , 10.71341704, 10.85780437, 10.83919067, 10.80099567,
       10.75020961, 11.04735455, 10.85516122, 10.81254426, 10.62579244,
       10.88302514, 10.84827799, 10.59262557, 10.4581895 , 10.13918098]) # The M_stellar, in log10(M_*/M_sun) unit.

plt.figure(figsize=(15, 11))
plt.scatter((10**Mstar/10**R_lumi)[zs>1.44], (host_mag_ACS - host_mag_WFC3)[zs>1.44],s=680, c ='green',
            marker="*",zorder=100, vmin=0, vmax=7, edgecolors='white', label='WFC/F140W')
plt.scatter((10**Mstar/10**R_lumi)[zs<1.44], (host_mag_ACS - host_mag_WFC3)[[zs<1.44]],s=280, c ='green',
            marker="o",zorder=100, vmin=0, vmax=7, edgecolors='white', label='WFC/F125W')
plt.ylabel("mag_$_{ACS}$ - mag_$_{WFC3}$",fontsize=35)
plt.xlabel("(M_*/M_sun)/(L_R/L_Rsun)",fontsize=35)
plt.tick_params(labelsize=25)
plt.ylim([1.2, 2.8])
plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':24},ncol=1)
plt.show()