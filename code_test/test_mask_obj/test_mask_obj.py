#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 15:07:38 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../../py_tools')
from mask_objects import mask_obj
QSO_img = pyfits.getdata('../../analysis_Astrodrz/CID216/analysis/CID216_cutout_outer.fits')
#QSO_img = pyfits.getdata('../../analysis_Astrodrz/CID1174/analysis/PSF1.fits')

target_mask, obj_mask = mask_obj(img=QSO_img, snr=2.5, exp_sz= 1.2)
