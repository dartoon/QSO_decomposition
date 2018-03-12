#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 11:34:57 2018

@author: Dartoon
"""
import sys
sys.path.insert(0,'../py_tools')
import astropy.io.fits as pyfits
from psfs_average import psf_ave, im_2_high_res
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from flux_profile import SB_profile, SB_compare


psf_NO=6 # The number of the psf.
for i in range(psf_NO):
    folder = '../fits_image/CID70/psf_cut_analysis/'
    fitsFile = pyfits.open(folder+'PSF{0}.fits'.format(i))
    PSF = fitsFile[0].data
    if i == 0 :
        psf_list = np.empty([psf_NO, PSF.shape[0], PSF.shape[1]])
        psf_list[0] = PSF
    else:
        psf_list[i] = PSF
#    PSFs= PSFs.append(PSF)

psf_final1=psf_ave(psf_list,mode = 'direct', not_count=(1,),
                  mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'])
psf_final2=psf_ave(psf_list,mode = 'direct', not_count=(1,4),
                  mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'])


plt.matshow(psf_final2, origin= 'low', norm=LogNorm())
plt.colorbar()
plt.show()

#SB_profile(image=img, center=np.asarray(psf_final1.shape)/2,ifplot=True)
#           fits_plot=True, if_mask = True, mask_NO=2, mask_reg=['test_circle.reg', 'test_box.reg'])

loadpsf = np.empty([2, psf_final1.shape[0], psf_final1.shape[1]])
loadpsf[0] = psf_final1
loadpsf[1] = psf_final2
SB_compare(QSO=psf_final2,psfs=psf_list, include_QSO=True, radius=30, grids=20)


#plt.matshow(psf_list[2], origin= 'low', norm=LogNorm())
#plt.show()

#psf_high = im_2_high_res(psf_list[1], shift_center = False)
#plt.matshow(psf_list[2], origin= 'low', norm=LogNorm())
#plt.show()
#plt.matshow(psf_high, origin= 'low', norm=LogNorm())
#plt.show()
