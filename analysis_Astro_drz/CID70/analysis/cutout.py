#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 14:04:02 2018

@author: Dartoon

Cut PSF and QSO for CID70
"""
import numpy as np
import sys
sys.path.insert(0,'../../../py_tools')
from cut_image import cut_image, cut_center_bright, save_loc_png, grab_pos
import copy
import astropy.io.fits as pyfits
ID = 'CID70'

#filename= 'stars_and_QSO.reg'
#c_psf_list = grab_pos(filename,reg_ty = 'astrodrz_06')
#print c_psf_list

fitsFile = pyfits.open('../astrodrz/final_drz.fits')
img = fitsFile[1].data  #- (-0.003)  # check the back grounp
center_QSO = np.array([1172,1805])
QSO = cut_center_bright(image=img, center=center_QSO, radius=100)
pyfits.PrimaryHDU(QSO).writeto('{0}_cutout.fits'.format(ID),overwrite=True)
count=0
#psf_list = copy.deepcopy(c_psf_list[:-1])
#psf_list = psf_list[psf_list[:,0].argsort()]
#psf_list[[3,4]] = psf_list[[4,3]]
psf_list = np.array([[2006,2058],[1334,1766],[1484,1674],[1002,1484],[2388,1784],[2140,1522]])
for i in range(len(psf_list)):
    PSF = cut_center_bright(image=img, center=psf_list[i], radius=60)
    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
    count += 1
    
#extra_psfs = np.array([[xxx,xxx],[xxx,xxx],[xxx,xxx],[xxx,xxx]])
#for i in range(len(extra_psfs)):
#    PSF = cut_center_bright(image=img, center=extra_psfs[i], radius=60)
#    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
#    count += 1
save_loc_png(img,center_QSO,psf_list, ID=ID,reg_ty = 'astrodrz_06')

##Check and find that the brightest point of PSF1.fits are not at the center.
#PSF = cut_image(image=img, center=(705, 843), radius=20)
#pyfits.PrimaryHDU(PSF).writeto('PSF1.fits'.format(i),overwrite=True)
