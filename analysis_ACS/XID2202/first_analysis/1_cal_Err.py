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

import os
path = os.getcwd()
ID = path.split('/')[-1]

exp = 2028

QSO_outer = pyfits.getdata('{0}_cutout_outer.fits'.format(ID))
from photutils import make_source_mask
mask = make_source_mask(QSO_outer, snr=2, npixels=5, dilate_size=11)
plt.imshow(QSO_outer* (1-mask*1), origin='low')
plt.show()
stdd = np.std(QSO_outer* (1-mask*1))
print "stdd:",stdd

filename = '{0}_cutout.fits'.format(ID)
image = pyfits.open(filename)
flux = image[0].data.copy()
rms =(flux/exp+stdd**2)**0.5
pyfits.PrimaryHDU(rms).writeto('wht_err.fits',overwrite=True) 