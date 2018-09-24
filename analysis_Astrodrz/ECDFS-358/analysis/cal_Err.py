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

ID = 'ECDFS-358'

wht = pyfits.getdata('wht_map.fits')
exp = 2400.
mean_wht = 2395.399 * (0.0642/0.135)**2
exp_map = exp * wht/mean_wht


stdd = 0.0050
filename = '{0}_cutout.fits'.format(ID)
image = pyfits.open(filename)
flux = image[0].data.copy()
rms =(flux/exp_map+stdd**2)**0.5
pyfits.PrimaryHDU(exp_map).writeto('exp_map.fits',overwrite=True) 
pyfits.PrimaryHDU(rms).writeto('wht_err.fits',overwrite=True) 