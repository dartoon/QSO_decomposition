#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 13:40:56 2018

@author: Dartoon

test the flux profile for QSO
"""

import numpy as np
import sys
sys.path.insert(0,'../../../py_tools')
import astropy.io.fits as pyfits
from flux_profile import flux_profile, SB_profile
import glob

ID='xxx'

test_n = 10
img = pyfits.getdata('{0}_cutout_outer.fits'.format(ID)
center = np.asarray(img.shape) /2
mask_list = glob.glob("QSO_outer_msk*.reg")   # Read *.reg files in a list.
r_flux, r_grids, regions=flux_profile(img, center,radius=center.min(), grids=50, ifplot=True, fits_plot= True, mask_list=mask_list)
print r_flux[-test_n:],'\n',r_flux[-test_n:][1:]-r_flux[-test_n:][:-1], '\n', (r_flux[-test_n:][1:]-r_flux[-test_n:][:-1]).sum()


#img = pyfits.getdata('PSF2.fits')
#center = np.asarray(img.shape) /2
#mask_list = glob.glob("PSF2*.reg")   # Read *.reg files in a list.
#r_flux, r_grids, regions=flux_profile(img, center,radius=center.min(), grids=30, ifplot=True, fits_plot= True, mask_list=mask_list)
#print r_flux[-test_n:],'\n',r_flux[-test_n:][1:]-r_flux[-test_n:][:-1], '\n', (r_flux[-test_n:][1:]-r_flux[-test_n:][:-1]).sum()

#img = pyfits.getdata('PSF3.fits')
#center = np.asarray(img.shape) /2
#mask_list = glob.glob("PSF3*.reg")   # Read *.reg files in a list.
#r_flux, r_grids, regions=flux_profile(img+0.00012, center,radius=center.min(), grids=30, ifplot=True, fits_plot= True, mask_list=mask_list)
#print r_flux[-test_n:],'\n',r_flux[-test_n:][1:]-r_flux[-test_n:][:-1], '\n', (r_flux[-test_n:][1:]-r_flux[-test_n:][:-1]).sum()
#
#img = pyfits.getdata('PSF4.fits')
#center = np.asarray(img.shape) /2
#mask_list = glob.glob("PSF4*.reg")   # Read *.reg files in a list.
#r_flux, r_grids, regions=flux_profile(img+0.00025, center,radius=center.min(), grids=30, ifplot=True, fits_plot= True, mask_list=mask_list)
#print r_flux[-test_n:],'\n',r_flux[-test_n:][1:]-r_flux[-test_n:][:-1], '\n', (r_flux[-test_n:][1:]-r_flux[-test_n:][:-1]).sum()
#

#img = psf_ave_pb
#center = np.asarray(img.shape) /2
#mask_list = glob.glob('')   # Read *.reg files in a list.
#r_flux, r_grids, regions=flux_profile(img, center,radius=center.min(), grids=30, ifplot=True, fits_plot= True, mask_list=mask_list)
#print r_flux[-test_n:],'\n',r_flux[-test_n:][1:]-r_flux[-test_n:][:-1], '\n', (r_flux[-test_n:][1:]-r_flux[-test_n:][:-1]).sum()
