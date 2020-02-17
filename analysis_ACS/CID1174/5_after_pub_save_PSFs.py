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
sys.path.insert(0,'../../py_tools')
from psfs_average import psf_ave
from flux_profile import QSO_psfs_compare, profiles_compare
from matplotlib.colors import LogNorm

import os
path = os.getcwd()
ID = path.split('/')[-1]

IDs = ['CID1174','CID255','CID50','CID70','XID2138','CID1281','CID3242','CID526',\
'LID1273','XID2202','CID206','CID3570','CID543','LID1538','XID2396','CID216','CID452',\
'CID597','LID360','CID237','CID454','CID607']

filt = 'acs'
# =============================================================================
# Read PSF and QSO image
# =============================================================================
QSO_bkg_value= 0.
QSO_im = pyfits.getdata('{0}_cutout.fits'.format(ID)) - QSO_bkg_value
QSO_msk = pyfits.getdata('{0}_msk.fits'.format(ID))
QSO_fm = len(QSO_im)
       
import pickle
PSFs_dict = {}
QSOs_dict = {}
for key in IDs:
    PSFs, _=pickle.load(open('../{0}/first_analysis/{0}_PSFs_QSO'.format(key),'rb'))
    PSFs_dict.update({'{0}'.format(key):PSFs})
#    QSOs_dict.update({'{0}'.format(key):QSOs})

PSF_list = []
PSF_id = []
for key in IDs:
    psfs_dict = PSFs_dict[key]
    psfs = [psfs_dict[i] for i in range(len(psfs_dict))]
    PSF_list += psfs
    name_id = [key+"_"+str(i) for i in range(len(psfs_dict))]
    PSF_id = PSF_id + name_id
    if len(PSF_list) != len(PSF_id):
        raise ValueError("The PSF_list is not consistent with PSF_id")
psf_list = [PSF_list[i][0] for i in range(len(PSF_list))]
PSF_mask_img_list = [PSF_list[i][3] for i in range(len(PSF_list))]
psf_name_list = PSF_id

count = 0
for i in np.array(range(len(psf_name_list))):
    print i, psf_name_list[i]
    psf_i = psf_list[i] * PSF_mask_img_list[i]
#    psf_i = psf_i[ct:-ct,ct:-ct]
    psf_i = psf_i/psf_i.sum()
    pyfits.PrimaryHDU(psf_i).writeto('../PSFs_used_ACS/psf_{0}.fits'.format(i),overwrite=True) 
   
