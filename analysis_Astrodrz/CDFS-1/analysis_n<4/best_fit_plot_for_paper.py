#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 16:02:19 2018

@author: Dartoon
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

import os
path = os.getcwd()
ID = path.split('/')[-2]

from filter_info import filt_info
filt = filt_info[ID]

import re
f = open("fit_result_each/each_PSF_fit_qso.txt","r")
string = f.read()
Chisq = re.findall(r"redu_Chisq':(.*?),",string)
Chisq = [float(i) for i in Chisq]
sort_Chisq = np.argsort(np.asarray(Chisq))
print "best fit PSF ID:", sort_Chisq[0]
psf_id = sort_Chisq[0]

# =============================================================================
# Read PSF and QSO image
# =============================================================================
QSO_bkg_value= 0.
QSO_im = pyfits.getdata('{0}_cutout.fits'.format(ID)) - QSO_bkg_value
QSO_msk = pyfits.getdata('{0}_msk.fits'.format(ID))
frame_size = 61
#frame = '{0}'.format(frame_size)
QSO_fm = len(QSO_im)
ct = (QSO_fm-frame_size)/2     # If want to cut to 61, i.e. (121-61)/2=30
        
import pickle
PSFs_dict = {}
QSOs_dict = {}
for key in filt_info.keys():
    PSFs, QSOs=pickle.load(open('../../{0}/analysis/{0}_PSFs_QSO'.format(key),'rb'))
    PSFs_dict.update({'{0}'.format(key):PSFs})
    QSOs_dict.update({'{0}'.format(key):QSOs})

PSF_list = []
PSF_id = []
filter_list = []
for key in PSFs_dict.keys():
    if filt_info[key] == filt:
        psfs_dict = PSFs_dict[key]
        psfs = [psfs_dict[i] for i in range(len(psfs_dict))]
        PSF_list += psfs
        name_id = [key+"_"+str(i) for i in range(len(psfs_dict))]
        PSF_id = PSF_id + name_id
        filt_ = [filt_info[key]]
        filter_list += filt_ * len(PSFs_dict[key])
        if len(PSF_list) != len(PSF_id):
            raise ValueError("The PSF_list is not consistent with PSF_id")
psf_list = [PSF_list[i][0] for i in range(len(PSF_list))]
PSF_mask_img_list = [PSF_list[i][3] for i in range(len(PSF_list))]
psf_name_list = PSF_id
# =============================================================================
# Doing the fitting
# =============================================================================
from fit_qso import fit_qso, fit_ps
from transfer_to_result import transfer_to_result
#from flux_profile import cr_mask_img
QSO_outer = pyfits.getdata('{0}_cutout_outer.fits'.format(ID))
from photutils import make_source_mask
mask = make_source_mask(QSO_outer, snr=2, npixels=5, dilate_size=11)
plt.imshow(QSO_outer* (1-mask*1), origin='low')
plt.close()
background_rms = np.std(QSO_outer* (1-mask*1))
print "background_rms: ", background_rms
QSO_msk = QSO_msk[ct:-ct,ct:-ct]
QSO_im = QSO_im[ct:-ct,ct:-ct]
QSO_msk = QSO_msk*0 +1    # This means no mask is added
QSO_std = pyfits.getdata('wht_err.fits')[ct:-ct,ct:-ct]
##############################Fit
source_params = None
import time
t1 = time.time()
fixcenter = False
#filename = 'fit_result_each/each_PSF_fit_qso.txt'
#if_file = glob.glob(filename)   
#if if_file == []:
#    fit_result =  open(filename,'w') 
#elif if_file is not []:
#    fit_result = open(filename,"r+")
#    fit_result.read()
count = 0
for i in np.array(range(psf_id, psf_id+1)):
    print "by PSF: {0}".format(psf_name_list[i])
    tag = '../../bestfit_plot/best_fit_{0}'.format(ID)
    psf_i = psf_list[i] * PSF_mask_img_list[i]
    psf_i = psf_i[ct:-ct,ct:-ct]
    source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_im, psf_ave=psf_i, psf_std = None,
                                                                     background_rms=background_rms,
                                                                     source_params=source_params, QSO_msk = QSO_msk, fixcenter=fixcenter,
                                                                     pix_sz = 'drz06', no_MCMC =True,
                                                                     QSO_std =QSO_std, tag=tag, image_plot=False)
    result = transfer_to_result(data=QSO_im, pix_sz = 'drz06',
            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, error_map=error_map,
            filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk =QSO_msk, tag=tag)
#    if count == 0:
#        fit_result.write("#QSO_img intensity: {0} \n".format(round(np.sum(QSO_im*QSO_msk),2)))
#    fit_result.write("#fit by PSF{0}: \n".format(psf_name_list[i]))
#    fit_result.write('PSF_intensity:{0} \n'.format(round(np.sum(psf_i),2)))
#    fit_result.write(repr(result) + "\n")
#    count += 1
#fit_result.close()
#t2 = time.time()
#t_tot= (t2-t1)/60
#print "total time:", t_tot, "mins"

#import os
#os.system('say "your program has finished"')