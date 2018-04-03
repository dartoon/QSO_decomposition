#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 12:35:45 2018

@author: Dartoon

test cut_image
"""

import numpy as np
import sys
sys.path.insert(0,'../py_tools')
import matplotlib.pylab as plt
#
#import astropy.io.fits as pyfits
#fitsFile = pyfits.open('psf.fits')
#img = fitsFile[0].data 
#
#from cut_image import cut_image, cut_center_bright
#box = cut_image(image=img, center=(49, 55), radius=10)
#plt.imshow(box, origin='low')
#plt.show()
#
#box = cut_center_bright(image=img, center=(49, 55), radius=10)
#plt.imshow(box, origin='low')
#plt.show()
#
#from cut_image import grab_pos
#filename= 'stars_and_QSO.reg'
#test = grab_pos(filename)