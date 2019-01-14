#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:37:04 2018

@author: Dartoon

For doing the k-correction
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from scipy import interpolate

def pass_dmag(z): 
    '''
    Calculate the changes of magntidues for the passive mag evolution.
    '''
    return 3.70*np.log10(1+z)

def k_corr_R(z, filt = 'F140w', galaxy_age = '5Gyrs'):
    filename = 'K_{0}.dat'.format(galaxy_age)
    with open('../material/{0}'.format(filename)) as f:
        content = f.readlines()
    lines = [x.strip() for x in content] 
    arr_lines = []
    arr_lines = [line.split( ) for line in lines if '#' not in line]
    K_array = np.array(arr_lines, dtype=float)
    if filt == 'F140w':
        k_grid = np.stack([K_array[:,0],K_array[:,1]])
    elif filt == 'F125w':
        k_grid = np.stack([K_array[:,0],K_array[:,2]])
    elif filt == 'F814w':
        k_grid = np.stack([K_array[:,0],K_array[:,3]])
    x = k_grid[0]
    y = k_grid[1]
    f = interpolate.interp1d(x, y)
    dm = f(z)
    return dm

##Plot TT's figure:
#z_d = np.linspace(1,2,20)
#k_c_140_5gy = k_corr_R(z_d, filt = 'F140w', galaxy_age = '5Gyrs')
#k_c_814_5gy = k_corr_R(z_d, filt = 'F814w', galaxy_age = '5Gyrs')
#k_c_140_1gy = k_corr_R(z_d, filt = 'F140w', galaxy_age = '1Gyrs')
#k_c_814_1gy = k_corr_R(z_d, filt = 'F814w', galaxy_age = '1Gyrs')
#plt.plot(z_d,k_c_140_5gy-k_c_814_5gy, label = '5gy', c='r')
#plt.plot(z_d,k_c_140_1gy-k_c_814_1gy, label = '1gy', c='b')
#plt.legend()
#plt.show()
