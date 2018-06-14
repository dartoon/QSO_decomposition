#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 20:17:41 2018

@author: Dartoon

On the analysis of xxx
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pylab as plt
import glob
import sys
sys.path.insert(0,'../../../py_tools')
from psfs_average import psf_ave
from flux_profile import QSO_psfs_compare, profiles_compare
from matplotlib.colors import LogNorm

ID = 'xxx'
filt = 'xxx'

# =============================================================================
# Read PSF and QSO image
# =============================================================================
psf_name_list = glob.glob("PSF*.fits")   # Read *.reg files in a list.
psf_list = []

#if not count PSF?, just name the file to not_count_PSF?.fits and +1 in the following line.
for i in range(len(psf_name_list)):
    if 'PSF{0}.fits'.format(i) in psf_name_list:
        psf_get = pyfits.getdata('PSF{0}.fits'.format(i))
        psf_list.append(psf_get)
    else:
        psf_list.append(None)

mask_list = glob.glob("PSF*.reg")   # Read *.reg files in a list.
QSO_im = pyfits.getdata('{0}_cutout.fits'.format(ID))

#==============================================================================
# Compare the profile and derive the Average image
#==============================================================================
cut = 40      #cut_range
if_QSO_l = [False, True]
gridsp_l = ['log', None]
if_annuli_l = [False, True] 
for i in range(2):
    for j in range(2):
        for k in range(2):
            plt_which_PSF = None
            plt_QSO = False
#            if i+k+j == 0:
#                plt_which_PSF = (0,1,2,3,4,5,6)
            if i==1 and j+k ==0:
                plt_QSO = True
            fig_psf_com = QSO_psfs_compare(QSO=QSO_im[cut:-cut,cut:-cut], psfs=psf_list,
                                               plt_which_PSF=plt_which_PSF,
                                               mask_list=mask_list, grids=40,
                                               include_QSO=if_QSO_l[i], 
                                               plt_QSO = plt_QSO, norm_pix = 3.0, astrodrz = True,
                                               gridspace= gridsp_l[j], if_annuli=if_annuli_l[k])
            fig_psf_com.savefig('PSFvsQSO{0}_{1}_{2}.pdf'.format(i,['xlog','xlin'][j],['circ','annu'][k]))
            if i==1 and k==1:
                plt.show()
            else:
                plt.close()

'''
psf_a, psf_a_std=psf_ave(psf_list,mode = 'CI', not_count=(4,5),
                  mask_list=mask_list)

psf_b, psf_b_std=psf_ave(psf_list,mode = 'CI', not_count=(4,5,6),
                  mask_list=mask_list)
#
#prf_list = [QSO_im,psf_a, psf_b]
#scal_list = [1,1,1]
#prf_name_list = ['QSO', 'Plan a', 'Plan b']
#fig_pro_compare = profiles_compare(prf_list, scal_list, prf_name_list=prf_name_list,norm_pix = 3.0,
#                                   gridspace = 'log',if_annuli=True,astrodrz=True)
#fig_pro_compare.savefig('PSFavd_vs_QSO_xlin_annu1.pdf')
#plt.show()

#pyfits.PrimaryHDU(psf_b).writeto('../../PSF_legacy/{0}_PSF.fits'.format(ID),overwrite=True)
#pyfits.PrimaryHDU(psf_b_std).writeto('../../PSF_legacy/{0}_PSF_std.fits'.format(ID),overwrite=True)

# =============================================================================
# Doing the fitting
# =============================================================================
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result
#from flux_profile import cr_mask_img
#mask_list = glob.glob("QSO_msk*.reg")   # Read *.reg files in a list.
#QSO_msk = cr_mask_img(QSO_im[cut:-cut,cut:-cut], mask_list, mask_reg_cut=cut)
QSO_msk =None
fit_result = open('fit_result.txt','w') 
background_rms = 0.01

##############################Fit
print "by psf_a"
fixcenter = True
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_a, psf_std = psf_a_std, QSO_msk=QSO_msk,
                                                       source_params=None, image_plot = True, corner_plot=True, flux_ratio_plot=True,
                                                       deep_seed = False, fixcenter=fixcenter,background_rms=background_rms, pix_sz = 'drz06')
result = transfer_to_result(data=QSO_im[cut:-cut,cut:-cut], pix_sz = 'drz06',
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=data_C_D,
                            cut=cut, filt=filt, fixcenter=fixcenter,ID=ID)
fit_result.write("#fit with PSF by Plan a: \n")
fit_result.write(repr(result) + "\n")

##############################Fit
print "by psf_a, relax center"
fixcenter = False
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_a, psf_std = psf_a_std, QSO_msk=QSO_msk,
                                                       source_params=None, image_plot = True, corner_plot=False, flux_ratio_plot=True,
                                                       deep_seed = False, fixcenter= fixcenter,background_rms=background_rms, pix_sz = 'drz06')
result = transfer_to_result(data=QSO_im[cut:-cut,cut:-cut],pix_sz = 'drz06',
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=data_C_D,
                            cut=cut, filt=filt, fixcenter=fixcenter,ID=ID, plot_compare= False)
fit_result.write("#fit with PSF by Plan a, relax center: \n")
fit_result.write(repr(result)+ "\n")
##############################Fit
print "by psf_b"
fixcenter = True
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_b, psf_std = psf_b_std, QSO_msk=QSO_msk,
                                                       source_params=None, image_plot = True, corner_plot=False, flux_ratio_plot=True,
                                                       deep_seed = False, fixcenter= fixcenter,background_rms=background_rms, pix_sz = 'drz06')
result = transfer_to_result(data=QSO_im[cut:-cut,cut:-cut],pix_sz = 'drz06',
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=data_C_D,
                            cut=cut, filt=filt, fixcenter=fixcenter,ID=ID, savepng=True, plot_compare= True)
fit_result.write("#fit with PSF by Plan b: \n")
fit_result.write(repr(result) + "\n")
##############################Fit
print "by psf_b, relax center"
fixcenter = False
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_b, psf_std = psf_b_std, QSO_msk=QSO_msk,
                                                       source_params=None, image_plot = True, corner_plot=False, flux_ratio_plot=True,
                                                       deep_seed = False, fixcenter= fixcenter,background_rms=background_rms, pix_sz = 'drz06')
result = transfer_to_result(data=QSO_im[cut:-cut,cut:-cut],pix_sz = 'drz06',
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=data_C_D,
                            cut=cut, filt=filt, fixcenter=fixcenter,ID=ID, plot_compare= False)
fit_result.write("#fit with PSF by Plan b, relax center: \n")
fit_result.write(repr(result)+ "\n")

fit_result.close()
'''
