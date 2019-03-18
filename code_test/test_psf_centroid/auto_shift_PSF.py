#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 23:18:20 2018

@author: Dartoon

PSF centorid
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
#import photutils

PSF = pyfits.open('PSF1.fits')[0].data.copy()
scl = len(PSF)
cut = 50
data = PSF[cut:-cut,cut:-cut]
plt.imshow(data,origin='low')
plt.show()


from photutils import centroid_com, centroid_1dg, centroid_2dg

x1, y1 = centroid_com(data)
x2, y2 = centroid_1dg(data)
x3, y3 = centroid_2dg(data)

fig, ax = plt.subplots(1, 1)
ax.imshow(data, origin='lower', interpolation='nearest', cmap='viridis')
marker = '+'
ms, mew = 30, 2.
plt.plot(x1, y1, color='y', marker=marker, ms=ms, mew=mew)
plt.plot(x2, y2, color='r', marker=marker, ms=ms, mew=mew)
plt.plot(x3, y3, color='b', marker=marker, ms=ms, mew=mew)
plt.show()

shift_x, shift_y = len(data)/2-x3, len(data)/2-y3
from lenstronomy.Util.kernel_util import de_shift_kernel
shift_PSF_it20 = de_shift_kernel(data, shift_x, shift_y, iterations=20)
#shift_PSF_it2 = de_shift_kernel(data, shift_x, 0, iterations=2)
plt.imshow(shift_PSF_it20,origin='low')
plt.show()