#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 13:46:11 2018

@author: Dartoon

Measure the total flux within 3 arcsec aperture.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import sys
sys.path.insert(0,'../py_tools')
from flux_profile import pix_region, flux_in_region

filt_info = [ 'CDFS-1', 'F140W',
'CDFS-229', 'F125W',
'CDFS-724', 'F125W',
'COSMOS-CID1174', 'F140W',
'COSMOS-CID206', 'F140W',
'COSMOS-CID216', 'F140W',
'COSMOS-CID3242', 'F140W',
'COSMOS-CID3570', 'F125W',
'COSMOS-CID452', 'F125W',
'COSMOS-CID454', 'F140W',
'COSMOS-CID50', 'F125W',
'COSMOS-CID607', 'F125W',
'COSMOS-CID70', 'F140W',
'COSMOS-LID1273', 'F140W',
'COSMOS-LID360', 'F140W',
'COSMOS-XID2138', 'F140W',
'COSMOS-XID2202', 'F140W',
'COSMOS-XID2396', 'F140W',
'ECDFS-358', 'F140W',
'SXDS-X1136', 'F125W',
'SXDS-X50', 'F125W',
'SXDS-X735', 'F140W',
'COSMOS-CID543', 'F125W',
'COSMOS-LID1538', 'F140W']

def return_flt(name, filt_info=filt_info):
    for i in range(len(filt_info)):
        if name in filt_info[i]:
            return filt_info[i+1]

QSO_list = glob.glob("CID*")
filename =  open('total_flux.txt','w')
filename.write("Target name,\t tot_flux, \t filter,\t magnitude \n")
for i in range(len(QSO_list)):
    QSO_image = pyfits.getdata('{0}/analysis/{0}_cutout.fits'.format(QSO_list[i]))
    c_p = len(QSO_image)/2
    deltapix = 0.0642
    radius = 3./2/0.0642
    region = pix_region(center=([c_p,c_p]), radius=radius)
    tot_flux =  flux_in_region(QSO_image, region)
    ax=plt.subplot(1,1,1)
    cb = ax.imshow(np.log10(QSO_image), origin='low')
    ax.add_patch(region.as_patch(facecolor='none', edgecolor='red'))
    plt.colorbar(cb)
    print QSO_list[i]
    plt.close()
    flt = return_flt(QSO_list[i])
    if flt == 'F140W':
        zp = 26.4524
    elif flt== 'F125W':
        zp = 26.2303
    mag = - 2.5 * np.log10(tot_flux) + zp 
    filename.write(QSO_list[i] +'\t' + repr(round(tot_flux,2)) +'\t' + flt +'\t'+ repr(round(mag,2))+'\n')
filename.close()