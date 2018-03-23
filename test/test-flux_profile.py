#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 10:44:08 2018

@author: Dartoon

test flux_profile.py
"""
import sys
sys.path.insert(0,'../py_tools')
import astropy.io.fits as pyfits
from flux_profile import pix_region, flux_in_region, cr_mask, SB_profile, QSO_psfs_compare
import matplotlib.pylab as plt
import numpy as np

# =============================================================================
# Example:
# =============================================================================
fitsFile = pyfits.open('psf.fits')
img = fitsFile[0].data 
region = pix_region(center=([49,49]), radius=5)
flux_in_region(img, region)
SB_profile(image=img, center=([49,49]),ifplot=False, fits_plot=False)

mask = cr_mask(image=img, filename='test_circle.reg')
plt.imshow((mask),origin='lower')
plt.show()

fitsFile = pyfits.open('psf.fits')
img = fitsFile[0].data 
region = pix_region(center=([49,49]), radius=5)
flux_in_region(img, region)
SB_profile(image=img, center=([49,49]),ifplot=True,
           fits_plot=True)

SB_profile(image=img, center=([49,49]),ifplot=True,
           fits_plot=True, mask_list=['test_circle.reg', 'test_box.reg'])


SB_profile(image=img, center=([40,40]),ifplot=True,
           fits_plot=True)

from flux_profile import text_in_string_list
print text_in_string_list("QSq2O1", ['QSO0asdw', '2adQSO2asdw', 'QaskQSO1]', '2adQSO2asdw', 'QaskQSO1]'])

#==============================================================================
#test .total_compare
#==============================================================================
from flux_profile import total_compare
data = pyfits.getdata('files/CID1174_cutout.fits')
QSO = pyfits.getdata('files/model_QSO.fits') 
host = pyfits.getdata('files/model_host.fits')
flux_list = [data, QSO, host]
label = ['data', 'QSO', 'host', 'model', 'residual']
import glob
mask_list = glob.glob("files/QSO*.reg")   # Read *.reg files in a list.
total_compare(label_list = label, flux_list = flux_list, target_ID = 'CID1174', data_mask_list = mask_list)
