#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 15:56:58 2018

@author: Dartoon

PSF analysis
"""

import sys
sys.path.insert(0,'../../../py_tools')
from cut_image import cut_center_bright, save_loc_png
import astropy.io.fits as pyfits
import numpy as np
from flux_profile import QSO_psfs_compare
import glob
ID = 'CID1174'

psf_NO=7 # The number of the psf.
for i in range(psf_NO):
    fitsFile = pyfits.open('PSF{0}.fits'.format(i))
    PSF = fitsFile[0].data
    if i == 0 :
        psf_list = np.empty([psf_NO, PSF.shape[0], PSF.shape[1]])
        psf_list[0] = PSF
    else:
        psf_list[i] = PSF
#    PSFs= PSFs.append(PSF)

fitsFile = pyfits.open('{0}_cutout.fits'.format(ID))
QSO = fitsFile[0].data
              
mask_list = glob.glob("PSF*.reg")   # Read *.reg files in a list.
              
#PSF_SB_compare(psf_list, radius=20, grids=20)
QSO_psfs_compare(QSO, psf_list, 
           mask_list=mask_list,
           plt_which_PSF=(0,1,2,3,4,5,6), include_QSO=True, radius=20, grids=20)
#QSO_psfs_compare(QSO, psf_list, include_QSO=True, radius=20, grids=20)

