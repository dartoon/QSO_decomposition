#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 14:04:02 2018

@author: Dartoon

Cut PSF and QSO for CID1174
"""
import numpy as np
import sys
sys.path.insert(0,'../../py_tools')
from cut_image import cut_image, cut_center_bright, save_loc_png, grab_pos
import copy
import astropy.io.fits as pyfits
from flux_profile import pix_region, flux_in_region
import matplotlib.pylab as plt

ID = 'CID1174'

filename= 'CID1174.reg'
c_psf_list = grab_pos(filename,reg_ty = 'astrodrz_06')
#print c_psf_list

fitsFile = pyfits.open('CID1174_sub_coadd.fits')
img = fitsFile[0].data # check the back grounp

center_QSO = c_psf_list[-1]

with open('CID1174_RaDec.reg', 'r') as myfile:
#    loc=myfile.read().replace('\n', '')
    content = myfile.readlines()
content.pop(0)
content.pop(0)
content.pop(0)

filename =  open('total_flux.txt','w')
filename.write("Target name, tot_flux, filter, magnitude, Ra Dec \n")
psf_list = copy.deepcopy(c_psf_list[:-1])
for i in range(len(psf_list)):
    PSF, PSF_center = cut_center_bright(image=img, center=psf_list[i], radius=60, return_center=True,  plot=False)
    c_p = len(PSF)/2
    deltapix = 0.0642
    radius = 3./2/0.0642
    region = pix_region(center=([c_p,c_p]), radius=radius)
    tot_flux =  flux_in_region(PSF, region)
    ax=plt.subplot(1,1,1)
    cb = ax.imshow(np.log10(PSF), origin='low')
    ax.add_patch(region.as_patch(facecolor='none', edgecolor='red'))
    plt.colorbar(cb)
    print 'PSF',i
    plt.show()
    zp = 26.4524
    mag = - 2.5 * np.log10(tot_flux) + zp 
    filename.write('PSF{0}'.format(i) +'\t' + repr(round(tot_flux,2))
    +'\t' + 'F140w' +'\t'+ repr(round(mag,2))+'\t'+content[i][7:28] +'\n')
filename.close()
                          

save_loc_png(img,center_QSO,psf_list, ID=ID,reg_ty = 'astrodrz_06')
