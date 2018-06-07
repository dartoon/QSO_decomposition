#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 14:04:02 2018

@author: Dartoon

Cut PSF and QSO for LID1538
"""
import numpy as np
import sys
sys.path.insert(0,'../../py_tools')
from cut_image import cut_image, cut_center_bright, save_other_target, grab_pos
import copy
import astropy.io.fits as pyfits
ID = 'LID1538'

filename= '{0}.reg'.format(ID)
c_psf_list = grab_pos(filename,reg_ty = 'acs')
#print c_psf_list

fitsFile = pyfits.open('../../Cycle25data/ACS_data/{0}_acs_I_mosaic_180mas_sci.fits'.format(ID))
img = fitsFile[0].data  #- (-0.003)  # check the back grounp
c_psf_list[[0,-1]] = c_psf_list[[-1,0]]
center_QSO = c_psf_list[-1]
QSO = cut_center_bright(image=img, center=center_QSO, radius=100)
#pyfits.PrimaryHDU(QSO).writeto('{0}_cutout.fits'.format(ID),overwrite=True)
count=0
psf_list = copy.deepcopy(c_psf_list[:-1])
psf_list = psf_list[psf_list[:,0].argsort()]
for i in range(len(psf_list)):
    PSF = cut_center_bright(image=img, center=psf_list[i], radius=60)
#    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
    count += 1
    
extra_psfs = np.array([[1277,3936],[2187,2797],[4106,4619],[4682,2523],[3585,1218],[4423,4082]])
for i in range(len(extra_psfs)):
    PSF = cut_center_bright(image=img, center=extra_psfs[i], radius=60)
#    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
    count += 1
    
center_star_brust = np.array([3607,4096])
save_other_target(img,center_QSO,other_target=center_star_brust,target_name = "ALMA starbusts",
                  c_psf_list=psf_list,extra_psfs=extra_psfs, ID=ID,reg_ty = 'acs')

##Check and find that the brightest point of PSF1.fits are not at the center.
#PSF = cut_image(image=img, center=(705, 843), radius=20)
#pyfits.PrimaryHDU(PSF).writeto('PSF1.fits'.format(i),overwrite=True)
