#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 14:04:02 2018

@author: Dartoon

Cut PSF and QSO for CID1174
"""
#import numpy as np
import sys
sys.path.insert(0,'../../../py_tools')
from cut_image import cut_image, cut_center_bright, save_loc_png
import astropy.io.fits as pyfits
ID = 'CID1174'

fitsFile = pyfits.open('../swarp/coadd.fits')
img = fitsFile[0].data - (-0.003)  # check the back grounp
center_QSO = (997, 915)
QSO = cut_center_bright(image=img, center=center_QSO, radius=50)
pyfits.PrimaryHDU(QSO).writeto('{0}_cutout.fits'.format(ID),overwrite=True)

c_psf_list = [(495, 892), (634, 791), (825, 919), (943, 676), (1048, 877), (1052, 341), (1301, 631)]

for i in range(len(c_psf_list)):
    PSF = cut_center_bright(image=img, center=c_psf_list[i], radius=20)
    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(i),overwrite=True)
save_loc_png(img,center_QSO,c_psf_list,ID=ID)

##Check and find that the brightest point of PSF1.fits are not at the center.
#PSF = cut_image(image=img, center=(705, 843), radius=20)
#pyfits.PrimaryHDU(PSF).writeto('PSF1.fits'.format(i),overwrite=True)
