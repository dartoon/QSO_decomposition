#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 16:45:53 2018

@author: Dartoon
"""

import numpy as np
#import os
import astropy.io.fits as pyfits
#import matplotlib.pylab as plt
import copy

noise_map = pyfits.open('noise_level.fits')[0].data.copy()

noise_bost = copy.deepcopy(noise_map)
noise_bost[noise_bost>0.019] = 10**5
plt.imshow(noise_bost)          
pyfits.PrimaryHDU(noise_bost).writeto('noise_boost.fits',overwrite=True)
