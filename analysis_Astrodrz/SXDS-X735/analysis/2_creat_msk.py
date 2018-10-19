#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 14:39:05 2018

@author: dartoon

Test create mask auto
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'../../../py_tools')
from mask_objects import mask_obj

import os
path = os.getcwd()
ID = path.split('/')[-2]

QSO_name = ID + "_cutout.fits"
QSO_img  = pyfits.getdata(QSO_name)
_,QSO_obj_mask  = mask_obj(img=QSO_img, exp_sz=1)
QSO_obj_mask = np.sum(QSO_obj_mask,axis=0)
QSO_obj_mask = (1 - (QSO_obj_mask != 0)*1.)
QSO_fit_mask0,_  = mask_obj(img=QSO_img, exp_sz=1.4)
QSO_msk0 = QSO_obj_mask*QSO_fit_mask0
QSO_fit_mask1,_  = mask_obj(img=QSO_img, exp_sz=2.4)
#QSO_msk1 = QSO_obj_mask*QSO_fit_mask1
QSO_msk1 = QSO_fit_mask1

print "QSO image:"
plt.imshow((QSO_img), norm=LogNorm(),origin='lower')
plt.show()
plt.imshow((QSO_msk0), origin='low') 
plt.show()
plt.imshow((QSO_msk1), origin='low') 
plt.show()
#pyfits.PrimaryHDU(QSO_msk0).writeto('{0}_msk0.fits'.format(ID),overwrite=True)
pyfits.PrimaryHDU(QSO_msk1).writeto('{0}_msk.fits'.format(ID),overwrite=True)

#import glob
#PSFs = glob.glob("*PSF?.fits") + glob.glob("*PSF??.fits")
#PSFs = sorted(PSFs,key=lambda x:x.split()[-1])
#for PSF_i in PSFs:
#    PSF_img = pyfits.getdata(PSF_i)
#    _, obj_mask  =  mask_obj(PSF_img) 
#    PSF_msk = np.sum(obj_mask,axis=0)
#    print "{0} image:".format(PSF_i)
##    plt.imshow((PSF_img), norm=LogNorm(),origin='lower')
##    plt.show()
#    wt_msk = (1 - (PSF_msk != 0)*1.)
##    plt.imshow((wt_msk), origin='low') 
##    plt.show()
#    pyfits.PrimaryHDU(wt_msk).writeto('{0}_msk.fits'.format(PSF_i[:-5]),overwrite=True)

#import glob
#from flux_profile import cr_mask_img
#psf_id = (6,7)
#for i in psf_id:
#    PSF_msk = pyfits.getdata('PSF{0}_msk.fits'.format(i))
#    ex_name = glob.glob("PSF{0}_?.reg".format(i))
#    PSF_ex_msk = cr_mask_img(PSF_msk, ex_name)
#    plt.imshow((PSF_ex_msk * PSF_msk), origin='low') 
#    plt.show()
#    pyfits.PrimaryHDU(PSF_ex_msk * PSF_msk).writeto('PSF{0}_msk.fits'.format(i),overwrite=True)
