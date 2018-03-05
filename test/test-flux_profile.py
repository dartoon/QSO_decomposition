#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 10:44:08 2018

@author: Dartoon

test flux_profile.py
"""
import astropy.io.fits as pyfits
from flux_profile import SB_profile
    
# =============================================================================
# Example:
# =============================================================================
#fitsFile = pyfits.open('psf.fits')
#img = fitsFile[0].data 
#region = pix_region(center=([49,49]), radius=5)
#flux_in_region(img, region)
#SB_profile(image=img, center=([49,49]),ifplot=False, fits_plot=False)

#mask = cr_mask(image=img, filename='test_circle.reg')
#plt.imshow((mask),origin='lower')
#plt.show()


fitsFile = pyfits.open('psf.fits')
img = fitsFile[0].data 
region = pix_region(center=([49,49]), radius=5)
flux_in_region(img, region)
#SB_profile(image=img, center=([49,49]),ifplot=True,
#           fits_plot=True)

SB_profile(image=img, center=([49,49]),ifplot=True,
           fits_plot=True, if_mask = True, mask_NO=2, mask_reg=['test_circle.reg', 'test_box.reg'])

#from flux_profile import text_in_string_list
#print text_in_string_list("QSq2O1", ['QSO0asdw', '2adQSO2asdw', 'QaskQSO1]', '2adQSO2asdw', 'QaskQSO1]'])

