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
for i in range(len(psf_name_list)+0):
    if 'PSF{0}.fits'.format(i) in psf_name_list:
        psf_get = pyfits.getdata('PSF{0}.fits'.format(i))
        psf_list.append(psf_get)
    else:
        psf_list.append(None)
frame_size = 61
frame = '{0}'.format(frame_size)
ct = (121-frame_size)/2     # If want to cut to 61, i.e. (121-61)/2=30
mask_list = glob.glob("PSF*.reg")   # Read *.reg files in a list.
QSO_im = pyfits.getdata('{0}_cutout.fits'.format(ID))
#==============================================================================
# Compare the profile and derive the Average image
#==============================================================================
if_QSO_l = [False, True]
gridsp_l = ['log', None]
if_annuli_l = [False, True] 
for i in range(2):
    for j in range(2):
        for k in range(2):
            plt_which_PSF = None
            plt_QSO = False
            if i+k+j == 0:
                plt_which_PSF = (0,1,2,3,4,5,6,7,8)
            if i==1 and j+k ==0:
                plt_QSO = True
            fig_psf_com = QSO_psfs_compare(QSO=QSO_im, psfs=psf_list,
                                               plt_which_PSF=plt_which_PSF,
                                               mask_list=mask_list, grids=30,
                                               include_QSO=if_QSO_l[i], 
                                               plt_QSO = plt_QSO, norm_pix = 3.0,
                                               gridspace= gridsp_l[j], if_annuli=if_annuli_l[k])
            fig_psf_com.savefig('PSFvsQSO{0}_{1}_{2}.pdf'.format(i,['xlog','xlin'][j],['circ','annu'][k]))
            if i==1 and k==1:
                plt.show()
            else:
                plt.close()

psf_ave_pa, psf_std_pa=psf_ave(psf_list,mode = 'CI', not_count=(0,1,3,5,6,7),
                  mask_list=mask_list)

psf_ave_pb, psf_std_pb=psf_ave(psf_list,mode = 'CI', not_count=(0,2,5,6,7,8),
                  mask_list=mask_list)

psf_ave_pa, psf_std_pa = psf_ave_pa[ct:-ct,ct:-ct], psf_std_pa[ct:-ct,ct:-ct]
psf_ave_pb, psf_std_pb = psf_ave_pb[ct:-ct,ct:-ct], psf_std_pb[ct:-ct,ct:-ct]

#pyfits.PrimaryHDU(psf_ave_pa).writeto('psf_ave_pa.fits',overwrite=True)
#pyfits.PrimaryHDU(psf_std_pa).writeto('psf_std_pa.fits',overwrite=True)
#pyfits.PrimaryHDU(psf_ave_pb).writeto('psf_ave_pb.fits',overwrite=True)
#pyfits.PrimaryHDU(psf_std_pb).writeto('psf_std_pb.fits',overwrite=True)

prf_list = [QSO_im,psf_ave_pa, psf_ave_pb]
scal_list = [1,1,1]
prf_name_list = ['QSO','Plan a','Plan b']
fig_pro_compare = profiles_compare(prf_list, scal_list, prf_name_list=prf_name_list,norm_pix = 5.0,
                                   gridspace = 'log',if_annuli=True,astrodrz=True)
plt.show()
# =============================================================================
# Doing the fitting
# =============================================================================
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result
from flux_profile import cr_mask_img
background_rms = 0.0076
#mask_list = glob.glob("QSO_msk*.reg")   # Read *.reg files in a list.
#QSO_msk = cr_mask_img(QSO_im, mask_list)
#QSO_msk = QSO_msk[ct:-ct,ct:-ct]
QSO_im = QSO_im[ct:-ct,ct:-ct]
QSO_msk = None
QSO_std = pyfits.getdata('wht_err.fits')[ct:-ct,ct:-ct]
fit_result = open('fit_result_{0}_reg.txt'.format(frame),'w') 
##############################Fit
print "by plan a"
fixcenter = True
tag = 'pa_fc{0}'.format(frame)
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im, psf_ave=psf_ave_pa, psf_std = psf_std_pa,background_rms=background_rms,
                                                       source_params=None, QSO_msk = QSO_msk, fixcenter=fixcenter, pix_sz = 'drz06', no_MCMC =True,
                                                       QSO_std =QSO_std, tag=tag)
result = transfer_to_result(data=QSO_im, pix_sz = 'drz06',
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=QSO_std**2,
                            filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk = '', tag=tag)
fit_result.write("#fit with plan a: \n")
fit_result.write(repr(result) + "\n")

##############################Fit
print "by plan a, relax center"
fixcenter = False
tag = 'pa_rc{0}'.format(frame)
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im, psf_ave=psf_ave_pa, psf_std = psf_std_pa,background_rms=background_rms,
                                                       source_params=None, deep_seed = False, fixcenter= fixcenter, pix_sz = 'drz06',  no_MCMC =True,
                                                       QSO_msk = QSO_msk, QSO_std =QSO_std, tag=tag)
result = transfer_to_result(data=QSO_im, pix_sz = 'drz06',
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=QSO_std**2,
                            filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk = '',tag=tag)
fit_result.write("#fit with plan b, relax center: \n")
fit_result.write(repr(result) + "\n")
##############################Fit
print "by plan b"
fixcenter = True
tag = 'pb_fc{0}'.format(frame)
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im, psf_ave=psf_ave_pb, psf_std = psf_std_pb,background_rms=background_rms,
                                                       source_params=None, deep_seed = False, fixcenter= fixcenter, pix_sz = 'drz06',
                                                       QSO_msk = QSO_msk, QSO_std =QSO_std, tag=tag)
result = transfer_to_result(data=QSO_im, pix_sz = 'drz06',
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=QSO_std**2,
                            filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk = '',tag=tag)
fit_result.write("#fit with plan b: \n")
fit_result.write(repr(result) + "\n")

###############################Fit
print "by plan b, relax center"
fixcenter = False
tag = 'pb_rc{0}'.format(frame)
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im, psf_ave=psf_ave_pb, psf_std = psf_std_pb,background_rms=background_rms,
                                                       source_params=None, deep_seed = False, fixcenter= fixcenter, pix_sz = 'drz06',
                                                       QSO_msk = QSO_msk, QSO_std =QSO_std, tag=tag)
result = transfer_to_result(data=QSO_im, pix_sz = 'drz06',
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=QSO_std**2,
                            filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk = '',tag=tag)
fit_result.write("#fit with plan b, relax center: \n")
fit_result.write(repr(result) + "\n")
fit_result.close()
