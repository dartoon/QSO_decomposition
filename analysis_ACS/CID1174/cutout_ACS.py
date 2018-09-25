#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 14:04:02 2018

@author: Dartoon

Cut PSF and QSO for xxx
"""
import numpy as np
import sys
sys.path.insert(0,'../../py_tools')
from cut_image import cut_image, cut_center_bright, save_loc_png, grab_pos
import copy
import astropy.io.fits as pyfits
ID = 'CID1174'

#fitsFile = pyfits.open('{0}_acs_I_mosaic_180mas_sci.fits'.format(ID))
#img = fitsFile[0].data # check the back grounp
#from astropy.visualization import SqrtStretch
#from astropy.stats import SigmaClip
#from photutils import Background2D, SExtractorBackground  
#from astropy.visualization.mpl_normalize import ImageNormalize
#import matplotlib.pyplot as plt
#norm = ImageNormalize(stretch=SqrtStretch())         
#sigma_clip = SigmaClip(sigma=3., iters=10)
#bkg_estimator = SExtractorBackground()
#from photutils import make_source_mask
#mask_0 = make_source_mask(img, snr=2, npixels=5, dilate_size=11)
#mask_1 = (np.isnan(img))
#mask = mask_0 + mask_1
#bkg = Background2D(img, (50, 50), filter_size=(3, 3),
#                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
#                   mask=mask)
#from matplotlib.colors import LogNorm
#fig=plt.figure(figsize=(15,15))
#ax=fig.add_subplot(1,1,1)
#ax.imshow(img, norm=LogNorm(), origin='lower') 
##bkg.plot_meshes(outlines=True, color='#1f77b4')
#ax.xaxis.set_visible(False)
#ax.yaxis.set_visible(False)
#plt.show()  
#fig=plt.figure(figsize=(15,15))
#ax=fig.add_subplot(1,1,1)
#ax.imshow(mask, origin='lower') 
##bkg.plot_meshes(outlines=True, color='#1f77b4')
#ax.xaxis.set_visible(False)
#ax.yaxis.set_visible(False)
#plt.show()  
#
#back = bkg.background* ~mask_1
#fig=plt.figure(figsize=(15,15))
#ax=fig.add_subplot(1,1,1)
#ax.imshow(back, origin='lower', cmap='Greys_r')
#ax.xaxis.set_visible(False)
#ax.yaxis.set_visible(False)
#plt.show()
#
#img -= back              
#pyfits.PrimaryHDU(img).writeto('sub_coadd.fits',overwrite=True)
img = pyfits.getdata('sub_coadd.fits')

filename= '{0}.reg'.format(ID)
c_psf_list, QSO_loc = grab_pos(filename,reg_ty = 'acs', QSO_reg_return=True)
center_QSO = c_psf_list[QSO_loc]

QSO, cut_center = cut_center_bright(image=img, center=center_QSO, radius=100, return_center=True, plot=True)
QSO_outer = cut_image(image=img, center=cut_center, radius=200)
pyfits.PrimaryHDU(QSO).writeto('{0}_cutout.fits'.format(ID),overwrite=True)
pyfits.PrimaryHDU(QSO_outer).writeto('{0}_cutout_outer.fits'.format(ID),overwrite=True)

count=0
psf_list = None
psf_list = np.delete(c_psf_list, (QSO_loc), axis=0)
dist = (psf_list-center_QSO)[:,0]**2+(psf_list-center_QSO)[:,1]**2
psf_list = psf_list[dist.argsort()]
#psf_list[[3,4]] = psf_list[[4,3]]
for i in range(len(psf_list)):
    print 'PSF',i
    PSF, PSF_center = cut_center_bright(image=img, center=psf_list[i], radius=60, return_center=True, plot=True)
    PSF_outer = cut_image(image = img, center = PSF_center, radius=200)
    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
    pyfits.PrimaryHDU(PSF_outer).writeto('PSF_outer_{0}.fits'.format(count),overwrite=True)
    count += 1

extra_psfs=None
#extra_psfs = np.array([[xxx,xxx],[xxx,xxx],[xxx,xxx],[xxx,xxx]])
#extra_psfs = extra_psfs[extra_psfs[:,0].argsort()]
#for i in range(len(extra_psfs)):
#    print 'PSF',count
#    PSF, PSF_center = cut_center_bright(image=img, center=extra_psfs[i], radius=60, return_center=True, plot=True)
#    PSF_outer = cut_image(image = img, center = PSF_center, radius=200)
#    pyfits.PrimaryHDU(PSF).writeto('PSF{0}.fits'.format(count),overwrite=True)
#    pyfits.PrimaryHDU(PSF_outer).writeto('PSF_outer_{0}.fits'.format(count),overwrite=True)
#    count += 1

if extra_psfs is None:
    save_loc_png(img,center_QSO,psf_list, ID=ID,reg_ty = 'acs')
else:
    save_loc_png(img,center_QSO,psf_list,extra_psfs, ID=ID,reg_ty = 'acs')
