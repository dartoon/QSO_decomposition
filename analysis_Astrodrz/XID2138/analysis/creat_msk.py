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
from cut_image import make_side_msk

ID='XID2138'
QSO_name = ID + "_cutout.fits"
QSO_img  = pyfits.getdata(QSO_name)
QSO_msk = make_side_msk(QSO_img,snr=2.5, npixels=10, dilate_size=5)

print "QSO image:"
plt.imshow((QSO_img), norm=LogNorm(),origin='lower')
plt.show()
plt.imshow((QSO_msk), origin='low') 
plt.show()
pyfits.PrimaryHDU(QSO_msk*1).writeto('{0}_msk.fits'.format(ID),overwrite=True)

import glob
PSFs = glob.glob("*PSF?.fits") + glob.glob("*PSF??.fits")
PSFs = sorted(PSFs,key=lambda x:x.split()[-1])
for PSF_i in PSFs:
    PSF_img = pyfits.getdata(PSF_i)
    PSF_msk =  make_side_msk(PSF_img,snr=2.5, npixels=10, dilate_size=5) 
    print "{0} image:".format(PSF_i)
    plt.imshow((PSF_img), norm=LogNorm(),origin='lower')
    plt.show()
    plt.imshow((PSF_msk), origin='low') 
    plt.show()
    pyfits.PrimaryHDU(PSF_msk*1).writeto('{0}_msk.fits'.format(PSF_i[:-5]),overwrite=True)

#
#QSO_mask = mask