#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 20:17:41 2018

@author: Dartoon

On the analysis of CID50
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

ID = 'CID50'
filt = 'F125w'

# =============================================================================
# Read PSF and QSO image
# =============================================================================
psf_name_list = glob.glob("PSF*.fits")   # Read *.reg files in a list.
psf_list = []
for i in range(len(psf_name_list)):
    psf_get = pyfits.getdata('PSF{0}.fits'.format(i))
    psf_list.append(psf_get)
mask_list = glob.glob("PSF*.reg")   # Read *.reg files in a list.
QSO_im = pyfits.getdata('{0}_cutout.fits'.format(ID))

#==============================================================================
# Compare the profile and derive the Average image
#==============================================================================
cut = 20      #cut_range
if_QSO_l = [False, True]
gridsp_l = ['log', None]
if_annuli_l = [False, True] 
for i in range(2):
    for j in range(2):
        for k in range(2):
            fig_psf_com = QSO_psfs_compare(QSO=QSO_im[cut:-cut,cut:-cut], psfs=psf_list,
#                                               plt_which_PSF=(0,1,2,3,4,5,6,7,8,9),
                                               mask_list=mask_list, grids=30,
                                               include_QSO=if_QSO_l[i], 
                                               gridspace= gridsp_l[j], if_annuli=if_annuli_l[k])
            fig_psf_com.savefig('PSFsvsQSO{0}_{1}_annu{2}.pdf'.format(i,['xlog','xlin'][j],k))

psf_ave_pa, psf_std_pa=psf_ave(psf_list,mode = 'CI', not_count=(9,6,7,3),
                  mask_list=mask_list)

psf_ave_pb, psf_std_pb=psf_ave(psf_list,mode = 'CI', not_count=(0,9,6,7,3,2,1,8),
                  mask_list=mask_list)

prf_list = [QSO_im,psf_ave_pa, psf_ave_pb]
scal_list = [1,1,1]
prf_name_list = ['QSO', 'Plan a', 'Plan b']
fig_pro_compare = profiles_compare(prf_list, scal_list, prf_name_list=prf_name_list, gridspace = None,if_annuli=True)
fig_pro_compare.savefig('PSFavd_vs_QSO_xlin_annu1.pdf')

pyfits.PrimaryHDU(psf_ave_pb).writeto('../../PSF_legacy/{0}_PSF.fits'.format(ID),overwrite=True)
pyfits.PrimaryHDU(psf_std_pb).writeto('../../PSF_legacy/{0}_PSF_std.fits'.format(ID),overwrite=True)
'''
# =============================================================================
# Doing the fitting
# =============================================================================
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result

fit_result = open('fit_result.txt','w') 
###############################Fit
#print "by plan a"
#fixcenter = True
#source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_ave_pa, psf_std = psf_std_pa,background_rms=0.037,
#                                                       source_params=None, image_plot = True, corner_plot=False, flux_ratio_plot=True,
#                                                       fixcenter=fixcenter)
#result = transfer_to_result(data=QSO_im[cut:-cut,cut:-cut],
#                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=data_C_D,
#                            cut=cut, filt=filt, fixcenter=fixcenter,ID=ID, savepng=True, plot_compare= True)  #If want to save this image
#fit_result.write("#fit with plan a, ave with 'CI': \n")
#fit_result.write(repr(result) + "\n")
##############################Fit
#print "by plan b"
#fixcenter = True
#source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_ave_pb, psf_std = psf_std_pb,background_rms=0.037,
#                                                       source_params=None, image_plot = True, corner_plot=True, flux_ratio_plot=True,
#                                                       deep_seed = False, fixcenter= fixcenter)
#result = transfer_to_result(data=QSO_im[cut:-cut,cut:-cut],
#                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=data_C_D,
#                            cut=cut, filt=filt, fixcenter=fixcenter,ID=ID)
#fit_result.write("#fit with plan b, ave with 'CI': \n")
#fit_result.write(repr(result) + "\n")

###############################Fit
print "by plan b, relax center, fix n == 1:"
fixcenter = False
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_ave_pb,background_rms=0.037,
                                                       source_params=None, image_plot = True, corner_plot=True, flux_ratio_plot=False,
                                                       deep_seed = False, fixcenter= fixcenter, fix_n = 1.)
result = transfer_to_result(data=QSO_im[cut:-cut,cut:-cut],
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=data_C_D,
                            cut=cut, filt=filt, fixcenter=fixcenter,ID=ID)
fit_result.write("#fit with plan b, ave with 'CI', fix n == 1: \n")
fit_result.write(repr(result) + "\n")

##############################Fit
print "by plan b, relax center"
fixcenter = False
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_ave_pb,background_rms=0.037,
                                                       source_params=None, image_plot = True, corner_plot=True, flux_ratio_plot=False,
                                                       deep_seed = False, fixcenter= fixcenter)
result = transfer_to_result(data=QSO_im[cut:-cut,cut:-cut],
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=data_C_D,
                            cut=cut, filt=filt, fixcenter=fixcenter,ID=ID)
fit_result.write("#fit with plan b, ave with 'CI', relax center: \n")
fit_result.write(repr(result) + "\n")
fit_result.close()
'''