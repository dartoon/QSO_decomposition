#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 13:55:40 2018

@author: Dartoon

test the PSF average
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

psf_name_list = glob.glob("PSF*.fits")   # Read *.reg files in a list.
psf_list = []
#if not count PSF?, just name the file to not_count_PSF?.fits and +1 in the following line.
for i in range(len(psf_name_list)+7):
    if 'PSF{0}.fits'.format(i) in psf_name_list:
        psf_get = pyfits.getdata('PSF{0}.fits'.format(i))
        psf_list.append(psf_get)
    else:
        psf_list.append(None)
frame_size = 61
frame = '{0}'.format(frame_size)
ct = (121-frame_size)/2     # If want to cut to 61, i.e. (121-61)/2=30
mask_list = glob.glob("PSF*.reg")   # Read *.reg files in a list.

import sys
sys.path.insert(0,'../../py_tools')
from psfs_average import psf_ave
psf_ave, psf_std=psf_ave(psf_list,mode = 'CI', not_count=(0,1,2,3,4,6,7),
                  mask_list=mask_list)

box_c = len(psf_list[5])/2
box_r = len(psf_list[5])/6
weights = np.zeros(len(psf_list))
weights = [np.sqrt(np.sum(psf_list[i][box_c-box_r:box_c+box_r,box_c-box_r:box_c+box_r])) for i in range(10) if psf_list[i][61,61]!=0]
weights_blank = [np.sqrt(np.sum(psf_list[i][box_c-box_r:box_c+box_r,box_c-box_r:box_c+box_r])) for i in range(10)]

PSFs = [psf_list[i]/weights_blank[i]**2 for i in range(10) if weights_blank[i]!=0]
weights = np.asarray(weights)
PSFs = np.asarray(PSFs)
print [PSFs[i][61,61] for i in range(3)]
print psf_ave[61,61]/psf_ave[60,60]
print np.sum(PSFs[:,61,61]*weights/(weights.sum()))/np.sum(PSFs[:,60,60]*weights/(weights.sum())) # the PSF is not normalized yet.

print psf_std[61,61]/psf_std[60,60] 

PSF_weight = np.sum([PSFs[i]*weights[i] for i in range(3)]/(weights.sum()),axis=0)
#print PSF_weight.shape
var = (PSFs - PSF_weight)**2
var_weighted = np.sum([var[i]*weights[i] for i in range(3)]/(weights.sum()),axis=0)
std = np.sqrt(var_weighted)
print std[61,61]/std[60,60]