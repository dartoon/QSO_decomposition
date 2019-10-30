#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 16:01:48 2019

@author: Dartoon

Test if I can find the color information.p
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
cmap = matplotlib.cm.get_cmap('viridis')
import matplotlib as mpl
mpl.rc('image', cmap='jet')
import glob


##The COSMOS files in http://www.mpia.de/homes/vdwel/candels.html by van der Wel et al. (2012).
##   NUMBER         RA        DEC          f        mag       dmag         re        dre          n         dn          q         dq         pa        dpa          sn
#file_galfit_COSMOS =  np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/cos_2epoch_wfc3_f160w_060mas_v1.0_galfit.cat')
#galfit_COSMOS_loc = file_galfit_COSMOS[:,[1,2]]  #RA, DEC
#galfit_COSMOS = file_galfit_COSMOS[:,[4,6,8,3]] # mag, re, n, flag
###The data from 3D HST: https://3dhst.research.yale.edu/Data.php  (PHOTOMETRY)
#file_stellar_COSMOS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/cosmos_3dhst.v4.1.cat')
#stellar_COSMOS_loc = file_stellar_COSMOS[:,[3,4]]  # RA, DEC
#stellar_COSMOS_flux_ap = file_stellar_COSMOS[:,[69,51,39]]  # flux, F140w, F125w, F814w.
#stellar_COSMOS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/cosmos_3dhst.v4.1.fout')[:,[1,6,8]]  # redshift, stellar mass, star formation rate

color_COSMOS_file = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/cosmos_3dhst.v4.1.cats/RF_colors/cosmos_3dhst.v4.1.master.RF') #U, V, J flux
color_COSMOS = color_COSMOS_file[:,[3, 7, 9]]  #Johnson_U, Johnson_V, 2MASS/J, SDSS/u, REST_FRAME/Bessel_V,  #For cosmos
## Check if one can recover the figure 1.
ColorUV_COSMOS = -(2.5* np.log10(color_COSMOS[:,0])-2.5* np.log10(color_COSMOS[:,1]))
ColorVJ_COSMOS = -(2.5* np.log10(color_COSMOS[:,1])-2.5* np.log10(color_COSMOS[:,2]))


#%%
#Define the region for blue and red galaxy:
#The three point: 
if_blue = []
for i in range(len(ColorVJ_COSMOS)):
    if ColorVJ_COSMOS[i] < 0.8425:
        if_blue.append((ColorUV_COSMOS[i] < 1.286))
    else:
        k = 1.17
        b = 0.3
        line_p = k*ColorVJ_COSMOS[i]+b
        if_blue.append((ColorUV_COSMOS[i] -line_p < 0))    
if_blue = np.asarray(if_blue)
fig, ax = plt.subplots(figsize=(7,7))
plt.plot(ColorVJ_COSMOS[if_blue], ColorUV_COSMOS[if_blue],'b.')
plt.plot(ColorVJ_COSMOS[np.logical_not(if_blue)], ColorUV_COSMOS[np.logical_not(if_blue)],'r.')
plt.xlim([-0.3, 2.3])
plt.ylim([-0.2, 2.695])
plt.tick_params(labelsize=25)
plt.xlabel("rest-frame V-J",fontsize=35)
plt.ylabel("rest-frame U-V",fontsize=35)
plt.show()
