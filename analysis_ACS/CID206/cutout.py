#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 14:04:02 2018

@author: Dartoon

Cut PSF and QSO for CID206
"""
import numpy as np
import sys
sys.path.insert(0,'../../py_tools')
from cut_image import cut_image, cut_center_bright, save_loc_png, grab_pos
import copy
import astropy.io.fits as pyfits
ID = 'CID206'

filename= '{0}.reg'.format(ID)
c_psf_list, QSO_loc = grab_pos(filename,reg_ty = 'acs', QSO_reg_return=True)
#print c_psf_list

fitsFile = pyfits.open('../../Cycle25data/ACS_data/{0}_acs_I_mosaic180mas_sci.fits'.format(ID))
img = fitsFile[0].data  #- (-0.003)  # check the back grounp
center_QSO = c_psf_list[QSO_loc]
QSO = cut_center_bright(image=img, center=center_QSO, radius=100)
pyfits.PrimaryHDU(QSO).writeto('{0}_cutout.fits'.format(ID),overwrite=True)
count=0
psf_list = np.delete(c_psf_list, (QSO_loc), axis=0)
psf_list = psf_list[psf_list[:,0].argsort()]
#psf_list[[3,4]] = psf_list[[4,3]]
for i in range(len(psf_list)):
    PSF = cut_center_bright(image=img, center=psf_list[i], radius=60)
    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
    count += 1
    
extra_psfs = np.array([[351,3092],[1055,3358]])
for i in range(len(extra_psfs)):
    PSF = cut_center_bright(image=img, center=extra_psfs[i], radius=60)
    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
    count += 1
save_loc_png(img,center_QSO,psf_list, extra_psfs, ID=ID,reg_ty = 'acs')

##Check and find that the brightest point of PSF1.fits are not at the center.
#PSF = cut_image(image=img, center=(705, 843), radius=20)
#pyfits.PrimaryHDU(PSF).writeto('PSF1.fits'.format(i),overwrite=True)
