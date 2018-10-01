#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 12:27:43 2018

@author: dartoon

Calculate the noise map level from the undrizzled image.
"""
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pylab as plt
#import glob
#from matplotlib.colors import LogNorm

ID = 'xxx'

wht = pyfits.getdata('wht_map.fits')
exp = 2395.399
mean_wht = exp * (0.0642/0.135)**2
exp_map = exp * wht/mean_wht

QSO_outer = pyfits.getdata('{0}_cutout_outer.fits'.format(ID))
from photutils import make_source_mask
mask = make_source_mask(QSO_outer, snr=2, npixels=5, dilate_size=11)
plt.imshow(QSO_outer* (1-mask*1), origin='low')
plt.show()
stdd = np.std(QSO_outer* (1-mask*1))
print "stdd:",stdd
#stdd = 0.0076

filename = '{0}_cutout.fits'.format(ID)
image = pyfits.open(filename)
flux = image[0].data.copy()
rms =(flux/exp_map+stdd**2)**0.5
pyfits.PrimaryHDU(exp_map).writeto('exp_map.fits',overwrite=True) 
pyfits.PrimaryHDU(rms).writeto('wht_err.fits',overwrite=True) 