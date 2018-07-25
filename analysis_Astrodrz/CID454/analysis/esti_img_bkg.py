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
from flux_profile import flux_profile, SB_profile, min_sub
import glob

ID='CID454'

#test_n=20
#img = pyfits.getdata('{0}_cutout.fits'.format(ID))
#center = np.asarray(img.shape) /2
#mask_list = glob.glob("QSO_msk*.reg")   # Read *.reg files in a list.
#sub, value = min_sub(0.001, img, mask_list=mask_list, test_n=test_n)
#r_flux, r_grids, regions=flux_profile(img-sub, center,radius=center.min(), grids=50, ifplot=True, fits_plot= True, mask_list=mask_list)
##print r_flux[-test_n:],'\n',r_flux[-test_n:][1:]-r_flux[-test_n:][:-1], '\n', 
#print sub
#print (r_flux[-test_n:][1:]-r_flux[-test_n:][:-1]).sum(), value

sub_list = np.zeros(10)
value = np.zeros(10)
for i in (1,3,4,5,6,7,8,9):
    print 'PSF{0}'.format(i)
    img = pyfits.getdata('PSF{0}.fits'.format(i))
    center = np.asarray(img.shape) /2
    mask_list = glob.glob("PSF{0}*.reg".format(i))   # Read *.reg files in a list.
    sub_list[i], value[i] = min_sub(0.001, img, mask_list=mask_list)
    print sub_list[i]
    r_flux, r_grids, regions=flux_profile(img- sub_list[i], center,
                                          radius=center.min(), grids=50,
                                          ifplot=True, fits_plot= True,
                                          mask_list=mask_list)
#    print r_flux[-test_n:],'\n',r_flux[-test_n:][1:]-r_flux[-test_n:][:-1], '\n'
#    print (r_flux[-test_n:][1:]-r_flux[-test_n:][:-1]).sum()
print sub_list, value