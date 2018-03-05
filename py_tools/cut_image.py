#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 11:45:45 2018

@author: Dartoon

The scripts for cut out the PSF candidate.
"""

import numpy as np
from regions import PixCoord, CirclePixelRegion 
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from flux_profile import pix_region

def cut_image(image, center, radius):
    region = pix_region(center, radius=radius)
    cut = region.to_mask(mode='exact')
    cut_image = cut.cutout(image)
    return cut_image

def cut_center_bright(image, center, radius):
    temp_center = np.asarray(center)
    radius = radius
    img_test = cut_image(image=image, center=temp_center, radius=radius)
    test_center =  np.asarray(np.where(img_test == img_test.max()))[:,0]
    center_shift = np.array((test_center- radius))[::-1]
    center = (temp_center + center_shift)
    cut_c_b = cut_image(image=image, center=center, radius=radius)
    return cut_c_b
    
def save_loc_png(img, center_QSO, c_psf_list):
    fig = plt.figure(figsize=(15,15))
    ax=fig.add_subplot(1,1,1)
    cax=ax.imshow(np.log10(img+0.15),origin='lower', cmap='gist_heat', interpolation='nearest', vmax=2.1)
    QSO_box_size = 30
    QSO_reg = pix_region(center_QSO, radius= QSO_box_size)
    QSO_mask = QSO_reg.to_mask(mode='center')
    ax.text(center_QSO[0]-2*QSO_box_size, center_QSO[1]+QSO_box_size, 'QSO', fontsize=20)
    ax.add_patch(QSO_mask.bbox.as_patch(facecolor='none', edgecolor='white'))
    PSF_box_size = 20
    for i in range(len(c_psf_list)):
        PSF_reg = pix_region(c_psf_list[i], radius= PSF_box_size)
        PSF_mask = PSF_reg.to_mask(mode='center')
        ax.add_patch(PSF_mask.bbox.as_patch(facecolor='none', edgecolor='blue'))
        ax.text(c_psf_list[i][0]-2*PSF_box_size, c_psf_list[i][1]+PSF_box_size, 'PSF{0}'.format(i), fontsize=15)
        fig.savefig('QSO_PSF_loc.png')
