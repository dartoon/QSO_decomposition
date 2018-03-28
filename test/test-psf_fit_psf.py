#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 11:22:50 2018

@author: Dartoon

fit psf with the initial final psf and see the shifting of the center.
"""

import sys
sys.path.insert(0,'../py_tools')
import astropy.io.fits as pyfits
from psfs_average import psf_ave
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from flux_profile import SB_profile, QSO_psfs_compare, profiles_compare
import time
from lenstronomy.Util.kernel_util import de_shift_kernel

psf_NO=6 # The number of the psf.
for i in range(psf_NO):
    folder = '../fits_image/CID70/analysis/'
    fitsFile = pyfits.open(folder+'PSF{0}.fits'.format(i))
    PSF = fitsFile[0].data
    if i == 0 :
        psf_list = np.empty([psf_NO, PSF.shape[0], PSF.shape[1]])
        psf_list[0] = PSF
    else:
        psf_list[i] = PSF

not_count = (1,4)
psf_init_ave, psf_std=psf_ave(psf_list,mode = 'CI', not_count=not_count,
                  mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'])

QSO_psfs_compare(QSO=psf_init_ave, psfs=psf_list,
                 plt_which_PSF=(4,),
                 mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'],
                 include_QSO=True, radius=20, grids=20)

shifted_psf_list = np.zeros_like(psf_list)
for i in range(psf_NO):
    fitted_PSF = psf_list[i]
    runfile('load_4_psf_fit_psf.py', wdir='../py_tools')
    start_time = time.time()
    lens_result, source_result, lens_light_result, ps_result, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(fitting_kwargs_list)
    end_time = time.time()
#    print(end_time - start_time, 'total time needed for computation')
#   print('============ CONGRATULATION, YOUR JOB WAS SUCCESSFUL ================ ')
    image_reconstructed, _, _, _ = imageModel.image_linear_solve(kwargs_source=source_result, kwargs_ps=ps_result)
    from lenstronomy.Util.kernel_util import de_shift_kernel
    shifted_psf_list[i] = de_shift_kernel(fitted_PSF, -ps_result[0]['ra_image'][0] , -ps_result[0]['dec_image'][0])
psf_final, psf_std=psf_ave(shifted_psf_list,mode = 'CI', not_count=not_count,
                  mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'])

plt.imshow(psf_final, origin = 'low', norm=LogNorm())
plt.show()
prf_list = [psf_init_ave,psf_final]
scal_list = [1,1]
from flux_profile import QSO_psfs_compare, profiles_compare

profiles_compare(prf_list, scal_list)

QSO_psfs_compare(QSO=psf_final, psfs=psf_list,
                 plt_which_PSF=(4,),
                 mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'],
                 include_QSO=True, radius=20, grids=20)