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
    if i ==2:
        psf_get = psf_get #-0.0001
    if i ==4:
        psf_get = psf_get #+0.00025
    psf_list.append(psf_get)
mask_list = glob.glob("PSF*.reg")   # Read *.reg files in a list.
QSO_im = pyfits.getdata('{0}_cutout.fits'.format(ID))
#==============================================================================
# Compare the profile and derive the Average image
#==============================================================================
#cut = 40      #cut_range
#if_QSO_l = [False, True]
#gridsp_l = ['log', None]
#if_annuli_l = [False, True] 
#for i in range(2):
#    for j in range(2):
#        for k in range(2):
#            plt_which_PSF = None
#            plt_QSO = False
##            if i+k+j == 0:
##                plt_which_PSF = (0,1,2,3,4,5,6)
#            if i==1 and j+k ==0:
#                plt_QSO = True
#            fig_psf_com = QSO_psfs_compare(QSO=QSO_im, psfs=psf_list,
#                                               plt_which_PSF=plt_which_PSF,
#                                               mask_list=mask_list, grids=30,
#                                               include_QSO=if_QSO_l[i], 
#                                               plt_QSO = plt_QSO, norm_pix = 5.0,
#                                               gridspace= gridsp_l[j], if_annuli=if_annuli_l[k])
##            fig_psf_com.savefig('PSFvsQSO{0}_{1}_{2}.pdf'.format(i,['xlog','xlin'][j],['circ','annu'][k]))
#            if i==1 and k==1:
#                plt.show()
#            else:
#                plt.close()

psf_ave_pa, psf_std_pa=psf_ave(psf_list,mode = 'CI', not_count=(0,),
                  mask_list=mask_list)

psf_ave_pb, psf_std_pb=psf_ave(psf_list,mode = 'CI', not_count=(0,1,5,6),
                  mask_list=mask_list)

psf_ave_pa, psf_std_pa = psf_ave_pa[30:-30,30:-30], psf_std_pa[30:-30,30:-30]
psf_ave_pb, psf_std_pb = psf_ave_pb[30:-30,30:-30], psf_std_pb[30:-30,30:-30]

#QSO_im -= 0.00024

prf_list = [QSO_im,psf_ave_pa, psf_ave_pb]
scal_list = [1,1,1]
prf_name_list = ['QSO','Plan a','Plan b']
fig_pro_compare = profiles_compare(prf_list, scal_list, prf_name_list=prf_name_list,norm_pix = 6.0,
                                   gridspace = 'log',if_annuli=True,astrodrz=True)
plt.show()
# =============================================================================
# Doing the fitting
# =============================================================================
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result
from flux_profile import cr_mask_img
background_rms = 0.0076
#QSO_im -= 0.00026
#mask_list = glob.glob("QSO_msk*.reg")   # Read *.reg files in a list.
#QSO_msk = cr_mask_img(QSO_im, mask_list)
mask_list = glob.glob('QSO_msk*.reg')   # Read *.reg files in a list.
QSO_msk = cr_mask_img(QSO_im, mask_list)[30:-30,30:-30]

#QSO_msk = QSO_msk[30:-30,30:-30]

QSO_std = pyfits.getdata('wht_err.fits')[30:-30,30:-30]
QSO_im  = QSO_im[30:-30,30:-30]
fit_result = open('fit_result_PSF_std_QSOmsk.txt','w') 
##############################Fit
#print "by plan a"
#fixcenter = True
#tag = 'plan_a_fc'
#source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im, psf_ave=psf_ave_pa, psf_std = psf_std_pa,background_rms=background_rms,
#                                                       source_params=None, QSO_msk = QSO_msk, fixcenter=fixcenter, pix_sz = 'drz06',
#                                                       QSO_std =QSO_std, tag = tag)
#result = transfer_to_result(data=QSO_im, pix_sz = 'drz06',
#                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=QSO_std**2,
#                            cut=0, filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk = '', tag = tag)
#fit_result.write("#fit with plan a, ave with 'CI': \n")
#fit_result.write(repr(result) + "\n")

##############################Fit
print "by plan a, relax center"
fixcenter = False
tag = 'PSF_std_Plan_a_QSOmsk'
source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_im, psf_ave=psf_ave_pa, psf_std = psf_std_pa,background_rms=background_rms,
                                                       source_params=None, deep_seed = False, fixcenter= fixcenter, pix_sz = 'drz06', no_MCMC=True,
                                                       QSO_msk = QSO_msk, QSO_std =QSO_std, tag=tag)
result = transfer_to_result(data=QSO_im, pix_sz = 'drz06',
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, error_map=error_map,
                            filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk = QSO_msk,tag=tag)
fit_result.write("#fit with plan a, ave with 'CI', relax center: \n")
fit_result.write(repr(result) + "\n")
##############################Fit
#print "by plan b"
#fixcenter = True
#tag = 'plan_b_fc'
#source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im, psf_ave=psf_ave_pb, psf_std = psf_std_pb,background_rms=background_rms,
#                                                       source_params=None, deep_seed = False, fixcenter= fixcenter, pix_sz = 'drz06',
#                                                       QSO_msk = QSO_msk, QSO_std =QSO_std)
#result = transfer_to_result(data=QSO_im, pix_sz = 'drz06',
#                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=QSO_std**2,
#                            cut=0, filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk = '',tag=tag)
#fit_result.write("#fit with plan b, ave with 'CI': \n")
#fit_result.write(repr(result) + "\n")

###############################Fit
print "by plan b, relax center"
fixcenter = False
tag = 'PSF_std_Plan_b_QSOmsk'
source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_im, psf_ave=psf_ave_pb, psf_std = psf_std_pb,background_rms=background_rms,
                                                       source_params=None, deep_seed = False, fixcenter= fixcenter, pix_sz = 'drz06', no_MCMC=True,
                                                       QSO_msk = QSO_msk, QSO_std =QSO_std, tag=tag)
result = transfer_to_result(data=QSO_im, pix_sz = 'drz06',
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, error_map=error_map,
                            filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk = QSO_msk,tag=tag)
fit_result.write("#fit with plan b, ave with 'CI', relax center: \n")
fit_result.write(repr(result) + "\n")
fit_result.close()
