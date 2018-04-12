#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 20:17:41 2018

@author: Dartoon

On the analysis of CID1174
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

ID = 'CID1174'
filt = 'F140w'
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
fig = QSO_psfs_compare(QSO=QSO_im[cut:-cut,cut:-cut], psfs=psf_list,
                 plt_which_PSF=(1,),
                 mask_list=mask_list,
                 include_QSO=False, grids=30, gridspace= 'log')
'''
psf_ave_dirt, psf_std_dirt=psf_ave(psf_list,mode = 'direct', not_count=(5,6),
                  mask_list=mask_list)

psf_ave_wght, psf_std_wght=psf_ave(psf_list,mode = 'CI', not_count=(5,6),
                  mask_list=mask_list)

#pyfits.PrimaryHDU(psf_ave_wght).writeto('../../PSF_legacy/{0}_PSF.fits'.format(ID),overwrite=True)
#pyfits.PrimaryHDU(psf_std_wght).writeto('../../PSF_legacy/{0}_PSF_std.fits'.format(ID),overwrite=True)

prf_list = [QSO_im,psf_ave_dirt, psf_ave_wght]
scal_list = [1,1,1]
prf_name_list = ['QSO', 'PSF_ave_direct', 'PSF_ave_by_wght']
profiles_compare(prf_list, scal_list, prf_name_list=prf_name_list, gridspace = 'log')

from fit_qso import fit_qso
from transfer_to_result import transfer_to_result

fit_result = open('fit_result.txt','w') 

##############################Fit
print "by psf_ave_dirt"
fixcenter = True
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_ave_dirt, psf_std = psf_std_dirt,
                                                       source_params=None, image_plot = False, corner_plot=False, flux_ratio_plot=True,
                                                       fixcenter=fixcenter)
result = transfer_to_result(data=QSO_im[cut:-cut,cut:-cut],
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=data_C_D,
                            cut=cut, filt=filt, fixcenter=fixcenter,ID=ID)
fit_result.write("#fit with averaged PSF by 'direct': \n")
fit_result.write(repr(result) + "\n")
##############################Fit
print "by psf_ave_wght"
fixcenter = True
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_ave_wght, psf_std = psf_std_wght,
                                                       source_params=None, image_plot = True, corner_plot=False, flux_ratio_plot=True,
                                                       deep_seed = False, fixcenter= fixcenter)
result = transfer_to_result(data=QSO_im[cut:-cut,cut:-cut],
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=data_C_D,
                            cut=cut, filt=filt, fixcenter=fixcenter,ID=ID, savepng=True, plot_compare= True)
fit_result.write("#fit with averaged PSF by 'CI': \n")
fit_result.write(repr(result) + "\n")
##############################Fit
print "by psf_ave_wght, relax center"
fixcenter = False
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_ave_wght, psf_std = psf_std_wght,
                                                       source_params=None, image_plot = False, corner_plot=False, flux_ratio_plot=True,
                                                       deep_seed = True, fixcenter= fixcenter)
result = transfer_to_result(data=QSO_im[cut:-cut,cut:-cut],
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=data_C_D,
                            cut=cut, filt=filt, fixcenter=fixcenter,ID=ID, plot_compare= False)
fit_result.write("#fit with averaged PSF by 'CI', relax center: \n")
fit_result.write(repr(result)+ "\n")

fit_result.close()
'''