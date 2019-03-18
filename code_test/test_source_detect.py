#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 15:59:14 2018

@author: dxh
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from photutils.datasets import make_100gaussians_image
data = make_100gaussians_image()

from photutils import detect_threshold
threshold = detect_threshold(data, snr=2.)

#threshold = bkg + (2.0 * bkg_rms)    

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources
sigma = 3.0 * gaussian_fwhm_to_sigma    # FWHM = 3.
kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
kernel.normalize()
segm = detect_sources(data, threshold, npixels=5, filter_kernel=kernel)

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
norm = ImageNormalize(stretch=SqrtStretch())
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
ax1.imshow(data, origin='lower', cmap='Greys_r', norm=norm)
ax1.set_title('Data')
ax2.imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345))
ax2.set_title('Segmentation Image')
plt.show()

from photutils import deblend_sources
segm_deblend = deblend_sources(data, segm, npixels=5,
                                filter_kernel=kernel, nlevels=32,
                                contrast=0.001)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
ax1.imshow(data, origin='lower', cmap='Greys_r', norm=norm)
ax1.set_title('Data')
ax2.imshow(segm_deblend, origin='lower', cmap=segm.cmap(random_state=12345))
ax2.set_title('Segmentation Image')
plt.show()
