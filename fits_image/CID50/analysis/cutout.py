#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 14:04:02 2018

@author: Dartoon

Cut PSF and QSO for CID50
"""
import numpy as np
import sys
sys.path.insert(0,'../../../py_tools')
from cut_image import cut_image, cut_center_bright, save_loc_png, grab_pos
import astropy.io.fits as pyfits
ID = 'CID50'

filename= 'stars_and_QSO.reg'
c_psf_list = grab_pos(filename)
print c_psf_list

fitsFile = pyfits.open('../swarp/coadd.fits')
img = fitsFile[0].data  #- (-0.003)  # check the back grounp
center_QSO = c_psf_list[-1]
QSO = cut_center_bright(image=img, center=center_QSO, radius=50)
pyfits.PrimaryHDU(QSO).writeto('{0}_cutout.fits'.format(ID),overwrite=True)

count=0
for i in range(len(c_psf_list[:-1])):
    PSF = cut_center_bright(image=img, center=c_psf_list[i], radius=30)
    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
    count += 1

extra_psfs = np.array([[641,381],[1158,398],[1195,726]])
for i in range(len(extra_psfs)):
    PSF = cut_center_bright(image=img, center=extra_psfs[i], radius=30)
    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
    count += 1
save_loc_png(img,center_QSO,c_psf_list[:-1],extra_psfs, ID=ID)

##Check and find that the brightest point of PSF1.fits are not at the center.
#PSF = cut_image(image=img, center=(705, 843), radius=20)
#pyfits.PrimaryHDU(PSF).writeto('PSF1.fits'.format(i),overwrite=True)
