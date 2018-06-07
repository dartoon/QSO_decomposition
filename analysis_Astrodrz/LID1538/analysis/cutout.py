#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 14:04:02 2018

@author: Dartoon

Cut PSF and QSO for LID1538
"""
import numpy as np
import sys
sys.path.insert(0,'../../../py_tools')
from cut_image import cut_image, cut_center_bright, save_other_target, grab_pos
import copy
import astropy.io.fits as pyfits
ID = 'LID1538'

filename= 'stars_and_QSO.reg'
c_psf_list = grab_pos(filename,reg_ty = 'astrodrz_06')
#print c_psf_list

fitsFile = pyfits.open('../astrodrz/final_drz.fits')
img = fitsFile[1].data  #- (-0.003)  # check the back grounp
center_QSO = np.array([1539,1210])
QSO = cut_center_bright(image=img, center=center_QSO, radius=100)
#pyfits.PrimaryHDU(QSO).writeto('{0}_cutout.fits'.format(ID),overwrite=True)
center_star_brust = np.array([1811,1726])
count=0
psf_list = copy.deepcopy(c_psf_list)
psf_list = psf_list[psf_list[:,0].argsort()]
for i in range(len(psf_list)):
    PSF = cut_center_bright(image=img, center=psf_list[i], radius=60)
#    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
    count += 1

extra_psfs = np.array([[737,1641],[1159,1110],[2054,1959],[2327,981],[1818,370]])
for i in range(len(extra_psfs)):
    PSF = cut_center_bright(image=img, center=extra_psfs[i], radius=60)
#    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
    count += 1
save_other_target(img,center_QSO,other_target=center_star_brust,target_name = "ALMA starbusts",
                  c_psf_list=psf_list,extra_psfs=extra_psfs, ID=ID,reg_ty = 'astrodrz_06')
##Check and find that the brightest point of PSF1.fits are not at the center.
#PSF = cut_image(image=img, center=(705, 843), radius=20)
#pyfits.PrimaryHDU(PSF).writeto('PSF1.fits'.format(i),overwrite=True)
