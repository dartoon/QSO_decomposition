#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 18:25:20 2018

@author: Dartoon

Compaer the PSFs Profile
"""
import sys
sys.path.insert(0,'../../py_tools/')
from flux_profile import QSO_psfs_compare, profiles_compare
import astropy.io.fits as pyfits
import glob
import numpy as np
PSF_name_list = glob.glob('*_PSF.fits')
psf_list = []
for i in range(len(PSF_name_list)):
#    if i ==5 or i==8:
       psf_get = pyfits.getdata('{0}'.format(PSF_name_list[i]))
       psf_list.append(psf_get)

fig_pro_compare = profiles_compare(psf_list, np.ones(len(psf_list)), prf_name_list=PSF_name_list,
                                   gridspace = None,if_annuli=True,norm_pix=3.0)#, y_log=True)
#plt.show()
