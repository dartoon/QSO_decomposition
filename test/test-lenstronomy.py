#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 10:19:08 2018

@author: Dartoon

test lenstronomy center and sharper PSF
"""

import astropy.io.fits as pyfits
import numpy as np
import glob
import copy
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

ID = 'CID1174'

psf_NO=7 # The number of the psf.
for i in range(psf_NO):
    fitsFile = pyfits.open('files/PSF{0}.fits'.format(i))
    PSF = fitsFile[0].data
    if i == 0 :
        psf_list = np.empty([psf_NO, PSF.shape[0], PSF.shape[1]])
        psf_list[0] = PSF
    else:
        psf_list[i] = PSF
#    PSFs= PSFs.append(PSF)

#mask_list = glob.glob("files/PSF*.reg")   # Read *.reg files in a list.
#psf_final=psf_ave(psf_list,mode = 'CI', not_count=(0,1,4,5,6),
#                  mask_list=mask_list)
#

#==============================================================================
# Test shift 
#==============================================================================
test_PSF = copy.deepcopy(psf_list[6])
from lenstronomy.Util.kernel_util import de_shift_kernel
test_shift_PSF = de_shift_kernel(test_PSF, 0.1 ,0)
plt.imshow(test_shift_PSF, origin = 'low', norm=LogNorm())
plt.show()

#==============================================================================
# Test subgrid
#==============================================================================
from lenstronomy.Util.kernel_util import subgrid_kernel
import sys
sys.path.insert(0,'../py_tools')
from flux_profile import profiles_compare
subgrid_PSF = subgrid_kernel(test_PSF, subgrid_res=11)
plt.imshow(subgrid_PSF, origin = 'low', norm=LogNorm())
plt.show()

prf_list = [test_PSF,subgrid_PSF]
scal_list = [1,11]
profiles_compare(prf_list, scal_list)

print np.where(subgrid_PSF == subgrid_PSF.max()), len(subgrid_PSF)/2