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
from flux_profile import flux_profile, min_sub, fit_bkg_as_gaussian
import astropy.io.fits as pyfits
import glob
import matplotlib.pylab as plt

ID='CID454'

#==============================================================================
# Measure QSO
#==============================================================================
#img_outer = pyfits.getdata('{0}_cutout_outer.fits'.format(ID))
#plt.imshow(np.log10(img_outer), origin='low') 
#plt.show()   
#QSO_bkg_value = fit_bkg_as_gaussian(img_outer)
#
#img = pyfits.getdata('{0}_cutout.fits'.format(ID))
#center = np.asarray(img.shape) /2
#mask_list = glob.glob("QSO_msk*.reg")   # Read *.reg files in a list.
#r_flux, r_grids, regions=flux_profile(img-QSO_bkg_value, center,radius=center.min(), grids=50, ifplot=True, fits_plot= True, mask_list=mask_list)
     

##==============================================================================
## Measure PSF
##==============================================================================
PSF_list = glob.glob("PSF_outer_*.fits")
n = len(PSF_list)
#PSF_bkg_values = np.zeros(n)
for i in range(5,7):
    print "FIT PSF{0}".format(i)
    img_outer = pyfits.getdata('PSF_outer_{0}.fits'.format(i))
    plt.imshow(np.log10(img_outer), origin='low') 
    plt.show()
#    PSF_bkg_values[i] = fit_bkg_as_gaussian(img_outer)
    name_PSF = glob.glob('*PSF{0}.fits'.format(i))
    img = pyfits.getdata(name_PSF[0])
    center = np.asarray(img.shape) /2
    mask_list = glob.glob("PSF{0}*.reg".format(i))   # Read *.reg files in a list.
    r_flux, r_grids, regions=flux_profile(img- PSF_bkg_values[i], center,
                                          radius=center.min(), grids=50,
                                          ifplot=True, fits_plot= False,
                                          mask_list=mask_list)
print QSO_bkg_value, PSF_bkg_values
#PSF_bkg_values[2] = np.average([PSF_bkg_values[i] for i in range(len(PSF_bkg_values)) if i!=2])
#PSF_bkg_values[5] /= 2
#PSF_bkg_values[6] /= 2