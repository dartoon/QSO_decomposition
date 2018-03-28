#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 11:34:57 2018

@author: Dartoon
"""
import sys
sys.path.insert(0,'../py_tools')
import astropy.io.fits as pyfits
from psfs_average import psf_ave
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from flux_profile import SB_profile, QSO_psfs_compare, profiles_compare

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
#    PSFs= PSFs.append(PSF)

psf_ave1, psf_std1 =psf_ave(psf_list,mode = 'direct', not_count=(1,),
                  mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'])

plt.imshow(psf_ave1, origin='low')
plt.colorbar()
plt.show()

plt.imshow(psf_std1, origin='low')
plt.colorbar()
plt.show()

psf_ave2, psf_std2=psf_ave(psf_list,mode = 'CI', not_count=(1,4),
                  mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'])

plt.imshow(psf_ave2, origin='low')
plt.colorbar()
plt.show()

plt.imshow(psf_std2, origin='low')
plt.colorbar()
plt.show()

loadpsf = np.empty([2, psf_ave1.shape[0], psf_ave1.shape[1]])
loadpsf[0] = psf_ave1
loadpsf[1] = psf_ave2
QSO_psfs_compare(QSO=psf_ave2, psfs=psf_list,
                 plt_which_PSF=(4,),
                 mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'],
                 include_QSO=True, radius=20, grids=20)

prf_list = [psf_ave1,psf_ave2]
scal_list = [1,1]
profiles_compare(prf_list, scal_list)

#==============================================================================
# test psf_average.psf_shift_ave
#==============================================================================
from fit_psf_pos import fit_psf_pos
from lenstronomy.Util.kernel_util import de_shift_kernel 
#print fit_psf_pos(psf_list[3],psf_ave2,psf_std =psf_std2)
#print fit_psf_pos(psf_list[3],psf_ave2) # Not include the PSF std
ra_image, dec_image = fit_psf_pos(psf_list[3],psf_ave2)
shifted_psf = de_shift_kernel(psf_list[3], -ra_image, -dec_image)
print fit_psf_pos(shifted_psf,psf_ave2) # The shifted PSF are at the center now


from psfs_average import psf_shift_ave
psf_final1, psf_final_std1 =  psf_shift_ave(psf_list,mode = 'direct', not_count=(1,4),
                    mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'])

psf_final2, psf_final_std2 =  psf_shift_ave(psf_list,mode = 'CI', not_count=(1,4),
                    mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'])

prf_list = [psf_ave1,psf_ave2,psf_final1, psf_final2]
scal_list = [1,1,1,1]
profiles_compare(prf_list, scal_list)

ra_image, dec_image = fit_psf_pos(psf_list[0],psf_final1)
shifted_psf = de_shift_kernel(psf_list[0], -ra_image, -dec_image)
print fit_psf_pos(shifted_psf,psf_final2)

psf_it_1, psf_std_it_1 =  psf_shift_ave(psf_list,mode = 'CI', not_count=(1,4),
                    mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg'])

psf_it_5, psf_std_it_5 =  psf_shift_ave(psf_list,mode = 'CI', not_count=(1,4),
                    mask_list=[folder+'PSF1_1.reg', folder+'PSF1_2.reg', folder+'PSF2_1.reg', folder+'PSF4_1.reg', folder+'PSF5_1.reg', folder+'PSF5_2.reg', folder+'PSF5_3.reg']
                    ,num_iter = 5)
prf_list = [psf_it_1,psf_it_5]
scal_list = [1,1]
profiles_compare(prf_list, scal_list)
plt.imshow(psf_it_5, norm = LogNorm(),origin='low')
plt.show() # See that the profile seem alike. But the psf_it_5 is a mass.

#from psfs_average import psf_shift_ave
#print psf_shift_ave(psf_list[0],psf_ave2, psf_std2)
#print psf_shift_ave(psf_list[0],psf_ave2, psf_std2)

#
#plt.matshow(psf_list[2], origin= 'low', norm=LogNorm())
#plt.show()
#
#psf_high = im_2_high_res(psf_list[1], shift_center = False)
#plt.matshow(psf_list[2], origin= 'low', norm=LogNorm())
#plt.show()
#plt.matshow(psf_high, origin= 'low', norm=LogNorm())
#plt.show()
