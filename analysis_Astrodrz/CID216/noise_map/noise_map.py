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
import glob
#from matplotlib.colors import LogNorm

explong=400.
stddlong=2.57

files = glob.glob('../astrodrz/OrIg_files/*')
n = 6 # In principle, six for dirzzling.

for i in range(n):
    image = pyfits.open(files[i])
    flux = image[1].data.copy()
#    plt.imshow(np.log10(flux), origin='lower',vmin=-1,vmax=2)
#    plt.colorbar()
#    plt.show()
    rms_sq =(flux/explong+stddlong**2)**0.5
#    hdu = pyfits.PrimaryHDU(rms_sq[1].data,header=input_header)
    new_hdul = pyfits.HDUList()
    new_hdul.append(pyfits.PrimaryHDU(header=image[0].header))
    new_hdul.append(pyfits.ImageHDU(rms_sq, header=image[1].header, name='SCI'))
    new_hdul.append(pyfits.ImageHDU(image[2].data.copy(), header=image[2].header, name='ERR'))
    new_hdul.append(pyfits.ImageHDU(image[3].data.copy(), header=image[3].header, name='DQ'))
    new_hdul.append(pyfits.ImageHDU(image[4].data.copy(), header=image[4].header, name='SAMP'))
    new_hdul.append(pyfits.ImageHDU(image[5].data.copy(), header=image[5].header, name='TIME'))
#    new_hdul.append(pyfits.ImageHDU(image[6].data.copy(), header=image[6].header, name='WCSCORR'))
    new_hdul.writeto('rms_sq_{0}_flt.fits'.format(i),overwrite=True)
#    image.close()