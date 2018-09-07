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

ID = 'XID2396'
filt = 'F140w'

# =============================================================================
# Read PSF and QSO image
# =============================================================================
QSO_bkg_value= 0.
psf_name_list = glob.glob("*PSF?.fits") + glob.glob("*PSF??.fits") # Read *.reg files in a list.
PSF_bkg_values = np.zeros(len(psf_name_list))
psf_list = []
#if not count PSF?, just name the file to not_count_PSF?.fits and +1 in the following line.
for i in range(len(psf_name_list)+0):
    if 'PSF{0}.fits'.format(i) in psf_name_list:
        psf_get = pyfits.getdata('PSF{0}.fits'.format(i)) - PSF_bkg_values[i]
        psf_list.append(psf_get)
    else:
        psf_list.append(None)
frame_size = 61
frame = '{0}'.format(frame_size)
ct = (121-frame_size)/2     # If want to cut to 61, i.e. (121-61)/2=30
PSF_mask_name_list = glob.glob("*PSF?_msk.fits") + glob.glob("*PSF??_msk.fits")   # Read *_msk.fits.
PSF_mask_name_list = sorted(PSF_mask_name_list,key=lambda x:x.split()[-1])
PSF_mask_img_list = [pyfits.getdata(PSF_mask_name_list[i]) for i in range(len(psf_list))]
QSO_im = pyfits.getdata('{0}_cutout.fits'.format(ID)) - QSO_bkg_value
QSO_msk = pyfits.getdata('{0}_msk.fits'.format(ID))

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
                plt_which_PSF = (0,1,2,3,4,5,6,7)
            if i==1 and j+k ==0:
                plt_QSO = True
            fig_psf_com = QSO_psfs_compare(QSO=QSO_im, QSO_msk=QSO_msk, psfs=psf_list,
                                               plt_which_PSF=plt_which_PSF,
                                               PSF_mask_img=PSF_mask_img_list, grids=30,
                                               include_QSO=if_QSO_l[i], 
                                               plt_QSO = plt_QSO, norm_pix = 5.0,
                                               gridspace= gridsp_l[j], if_annuli=if_annuli_l[k])
            fig_psf_com.savefig('PSFvsQSO{0}_{1}_{2}.pdf'.format(i,['xlog','xlin'][j],['circ','annu'][k]))
            if i==1 and k==1:
                plt.show()
            else:
                plt.close()
# =============================================================================
# Doing the fitting
# =============================================================================
from fit_qso import fit_qso, fit_ps
from transfer_to_result import transfer_to_result
from flux_profile import cr_mask_img
background_rms = 0.007
QSO_msk = QSO_msk[ct:-ct,ct:-ct]
QSO_im = QSO_im[ct:-ct,ct:-ct]
#QSO_msk = None
QSO_std = pyfits.getdata('wht_err.fits')[ct:-ct,ct:-ct]
#############################Fit
filename = 'fit_result_each/each_PSF_fit_ps.txt'
if_file = glob.glob(filename)   
if if_file == []:
    fit_result =  open(filename,'w') 
elif if_file is not []:
    fit_result = open(filename,"r+")
    fit_result.read()
count = 0
for i in np.array([0,1,2,3,4,5,6,7]):
    print "by PSF{0}".format(i)
    tag = 'fit_result_each/ps_fit_PSF{0}'.format(i)
    mask_list = glob.glob("PSF{0}_*.reg".format(i))
    print "mask_list", mask_list
    psf_i = psf_list[i] * PSF_mask_img_list[i]
    psf_i = psf_i[ct:-ct,ct:-ct]
    source_result, ps_result, image_ps, image_host, error_map=fit_ps(QSO_im, psf_ave=psf_i, psf_std = None,
                                                                     background_rms=background_rms,
                                                                     source_params=None, QSO_msk = QSO_msk, fixcenter=False,
                                                                     pix_sz = 'drz06', no_MCMC =True,
                                                                     QSO_std =QSO_std, tag=tag)
    if count == 0:
        fit_result.write("#QSO_img intensity: {0} \n".format(round(np.sum(QSO_im*QSO_msk),2)))
    fit_result.write("#fit by PSF{0}: \n".format(i))
    fit_result.write('PSF_intensity: {0} '.format(round(np.sum(psf_i),2)))
    ps_amp = round(ps_result[0]['point_amp'][0],2)
    host_amp = round((image_host*QSO_msk).sum(),2)
    ratio = round(host_amp/(host_amp+ps_amp)*100, 2)
    fit_result.write('Point_source_flux: '+repr(ps_amp) + ' host_flux: ' + repr(host_amp)+ ' host_ratio: '+ repr(ratio) + "%\n")
    count += 1
fit_result.close()    

fixcenter = False
filename = 'fit_result_each/each_PSF_fit_qso.txt'
if_file = glob.glob(filename)   
if if_file == []:
    fit_result =  open(filename,'w') 
elif if_file is not []:
    fit_result = open(filename,"r+")
    fit_result.read()
count = 0
for i in np.array([0,1,2,3,4,5,6,7]):
    print "by PSF{0}".format(i)
    tag = 'fit_result_each/qso_fit_PSF{0}'.format(i)
    mask_list = glob.glob("PSF{0}_*.reg".format(i))
    print mask_list
    psf_i = psf_list[i] * PSF_mask_img_list[i]
    psf_i = psf_i[ct:-ct,ct:-ct]
    source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_im, psf_ave=psf_i, psf_std = None,
                                                                     background_rms=background_rms,
                                                                     source_params=None, QSO_msk = QSO_msk, fixcenter=fixcenter,
                                                                     pix_sz = 'drz06', no_MCMC =True,
                                                                     QSO_std =QSO_std, tag=tag)
    result = transfer_to_result(data=QSO_im, pix_sz = 'drz06',
            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, error_map=error_map,
            filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk =QSO_msk, tag=tag)
    if count == 0:
        fit_result.write("#QSO_img intensity: {0} \n".format(round(np.sum(QSO_im*QSO_msk),2)))
    fit_result.write("#fit by PSF{0}: \n".format(i))
    fit_result.write('PSF_intensity:{0} \n'.format(round(np.sum(psf_i),2)))
    fit_result.write(repr(result) + "\n")
    count += 1
fit_result.close()

#from fit_qso import fit_single_host
#from roundme import roundme
#filename = 'fit_result_each/each_PSF_pure_galaxy.txt'
#if_file = glob.glob(filename)   
#if if_file == []:
#    fit_result =  open(filename,'w') 
#elif if_file is not []:
#    fit_result = open(filename,"r+")
#    fit_result.read()
#count = 0
#for i in np.array([???]):
#    print "by PSF{0}".format(i)
#    tag = 'fit_result_each/pure_galaxy_PSF{0}'.format(i)
#    mask_list = glob.glob("PSF{0}_*.reg".format(i))
#    print mask_list
#    psf_i = psf_list[i] * cr_mask_img(image=psf_list[i], mask_list=mask_list)
#    psf_i = psf_i[ct:-ct,ct:-ct]
#    lens_light_result, image_host, error_map = fit_single_host(QSO_im, psf_ave=psf_i,
#                                                              pix_sz = 'drz06',QSO_msk=QSO_msk,
#                                                              QSO_std=QSO_std, tag=tag)
#    if count == 0:
#        fit_result.write("#QSO_img intensity: {0} \n".format(round(np.sum(QSO_im*QSO_msk),2)))
#    fit_result.write("#fit by PSF{0}: \n".format(i))
#    fit_result.write('PSF_intensity:{0} \n'.format(round(np.sum(psf_i),2)))
#    fit_result.write(repr(roundme(lens_light_result[0])) + "\n")
#    fit_result.write('Image_total_flux: '+repr(round(QSO_im.sum(),2))+' Galaxy_fitted_totflux: '+repr(round(image_host.sum(),2)))
#    count += 1
#fit_result.close()

##==============================================================================
## Combining fitting
##==============================================================================
#PSF_mask_list = glob.glob("PSF*.reg")   # Read *.reg files in a list.
#count_ele = range(len(psf_name_list))
#import re
#count_list= [[xxx], [xxx],[xxx]]
#psf_com = [re.sub('\W+', '',repr(count_list[i])) for i in range(len(count_list))] 
##psf_com = ['0157', '57','01579', '579']
#
#not_count_list= [] #[(2,3,4,6,8,9), (0,1,2,3,4,6,8,9),(2,3,4,6,8) ,(0,1,2,3,4,6,8)]
#
#for i in range(len(count_list)):
#    count_i = [x for x in count_ele if x not in count_list[i]]
#    not_count_list.append(count_i)    
    
#PSF_aves, PSF_stds = [], []
#for i in range(len(not_count_list)):
#    psf_ave_i, psf_std_i=psf_ave(psf_list,mode = 'CI', not_count=not_count_list[i],   #0157
#                  mask_img_list = PSF_mask_img_list)
#    psf_ave_i = psf_ave_i[ct:-ct,ct:-ct]
#    psf_std_i = psf_std_i[ct:-ct,ct:-ct]
#    PSF_aves.append(psf_ave_i)
#    PSF_stds.append(psf_std_i)
#    
#prf_list = [QSO_im]+ PSF_aves
#scal_list = np.ones(len(PSF_aves)+1)
#prf_name_list = ['QSO']+ psf_com
#fig_pro_compare = profiles_compare(prf_list, scal_list, prf_name_list=prf_name_list,norm_pix = 5.0,
#                                   gridspace = 'log',if_annuli=True,astrodrz=True)
#plt.show()
#fig_pro_compare.savefig('pro_compare.pdf')
#
#
#filename = 'fit_result/comb_PSF.txt'
#if_file = glob.glob(filename)   
#if if_file == []:
#    fit_result =  open(filename,'w') 
#elif if_file is not []:
#    fit_result = open(filename,"r+")
#    fit_result.read()
#for i in range(len(PSF_aves)):
#    print "NO.",i," For fit PSF", psf_com[i]
#    pyfits.PrimaryHDU(PSF_aves[i]).writeto('psf_ave_{0}.fits'.format(psf_com[i]),overwrite=True)
#    pyfits.PrimaryHDU(PSF_stds[i]).writeto('psf_std_{0}.fits'.format(psf_com[i]),overwrite=True)
#    
#    tag = 'fit_result/comb_PSF_{0}'.format(psf_com[i])
#    source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_im, psf_ave=PSF_aves[i], psf_std = PSF_stds[i]**2,
#                                                                      background_rms=background_rms, pix_sz = 'drz06',  no_MCMC =True,
#                                                                      QSO_msk = QSO_msk, QSO_std =QSO_std, tag=tag)
#    result = transfer_to_result(data=QSO_im, pix_sz = 'drz06', source_result=source_result, ps_result=ps_result,
#                                image_ps=image_ps, image_host=image_host, error_map=error_map,
#                                filt=filt, fixcenter=False,ID=ID,QSO_msk = QSO_msk,tag=tag)
#    fit_result.write("#fit selected PSF{0}: \n".format(psf_com[i]))
#    fit_result.write(repr(result) + "\n")
#fit_result.close()

import os
os.system('say "your program has finished"')
