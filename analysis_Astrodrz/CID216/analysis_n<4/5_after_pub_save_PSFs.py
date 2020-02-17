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

import os
path = os.getcwd()
ID = path.split('/')[-2]

from filter_info import filt_info
filt = filt_info[ID]
del filt_info['CID255']
# =============================================================================
# Read PSF and QSO image
# =============================================================================
QSO_bkg_value= 0.
QSO_im = pyfits.getdata('{0}_cutout.fits'.format(ID)) - QSO_bkg_value
QSO_msk = pyfits.getdata('{0}_msk.fits'.format(ID))
frame_size = 81
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
PSF_keys = PSFs_dict.keys()
index1 = PSF_keys.index('CID1281')
del PSF_keys[index1]
index2 = PSF_keys.index('CID597')
del PSF_keys[index2]
PSF_keys.append('CID1281')
PSF_keys.append('CID597')

for key in PSF_keys:
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
pix_s = 0.0642 
for i in np.array(range(len(psf_name_list))):
    psf_i = psf_list[i] * PSF_mask_img_list[i]
#    psf_i = psf_i[ct:-ct,ct:-ct]
    psf_i = psf_i/psf_i.sum()
    pyfits.PrimaryHDU(psf_i).writeto('PSFs_used_F140w/psf_{0}.fits'.format(i),overwrite=True) 
    
