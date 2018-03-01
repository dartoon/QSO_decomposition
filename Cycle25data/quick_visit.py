#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 09:46:58 2018

@author: Dartoon
"""

import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import numpy as np
filename='CID452/idnl06jlq_ima.fits'
image = pyfits.open(filename)#[1].data.copy()


input_header= image[1].header
flux= image[1].data.copy()

#plt.imshow(np.log10(flux),origin='lower',vmax=2)
#plt.colorbar()
#plt.show()

hdu = pyfits.PrimaryHDU(flux,header=input_header)
hdu.writeto('/Users/Dartoon/Desktop/dither/idnl06jlq_ima.fits',overwrite=True)
