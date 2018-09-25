#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 09:48:34 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

filt = ['F125w', 'F140w', 'F140w', 'F140w', 'F140w', 'F125w', 'F140w', 'F125w',
        'F125w', 'F140w', 'F140w', 'F120w', 'F140w', 'F140w', 'F140w', 'F140w', 'F140w', 'F140w']

mag= np.array([23.574, 22.772, 23.099, 21.506, 21.258, 20.992, 21.547, 21.427,
               21.366, 21.669, 21.111, 21.198, 21.447, 20.609, 21.06, 21.21, 20.609, 20.129])

host_total_ratio = np.array([1.50,8.87,11.10,95.86,30.57,90.61,26.15,52.28,37.42,34.86,48.50,75.00,18.48,70.89,51.99,71.60,57.51,78.48])
host_total_ratio /= 100.

host_flux = np.zeros_like(mag)
for i in range(len(filt)):
    if filt[i] == 'F140w':
        zp = 26.4524
    elif filt[i] == 'F125w':
        zp = 26.2303
    host_flux[i] = 10.**(-0.4*(mag[i]-zp))

AGN_flux = host_flux/host_total_ratio - host_flux