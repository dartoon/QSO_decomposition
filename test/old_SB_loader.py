#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:47:59 2018

@author: dartoon

Read the surface brightness
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

minorLocator = AutoMinorLocator()

def radial_profile(data, center):
    x, y = np.indices((data.shape))  # derive (build) each x and y for the each pixel. 
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)  # derive the distance of each pixel to the center.
    r = r.astype(np.int) #remove the decimal part.

    tbin = np.bincount(r.ravel(), weights=data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 

    
fitsFile = fits.open('psf.fits')
img = fitsFile[0].data
#img[np.isnan(img)] = 0

#center = np.unravel_index(img.argmax(), img.shape)
center = (49, 49)
rad_profile = radial_profile(img, center)

fig, ax = plt.subplots()
plt.plot(rad_profile[0:35], 'x-')

ax.xaxis.set_minor_locator(minorLocator)

plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4, color='r')
plt.grid()
ax.set_ylabel(fitsFile[0].header['Label'] + " (" + fitsFile[0].header['BUNIT'] + ")")
ax.set_xlabel("Pixels")
plt.grid(which="minor")
plt.show()