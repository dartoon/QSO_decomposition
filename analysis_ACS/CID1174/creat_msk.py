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
sys.path.insert(0,'../../py_tools')
from cut_image import make_side_msk, make_circle_msk
from photutils import make_source_mask
ID='CID1174'
QSO_name = ID + "_cutout.fits"
QSO_img  = pyfits.getdata(QSO_name)
#QSO_msk = make_side_msk(QSO_img,snr=2.5, npixels=6, dilate_size=6, ct_QSO_mask = True)
#rad = np.sqrt(QSO_msk.sum()/np.pi)
#QSO_msk_c = make_circle_msk(QSO_img,x=len(QSO_img)/2, y=len(QSO_img)/2, radius=rad*1.2)
#
#print "QSO image:"
#plt.imshow((QSO_img), norm=LogNorm(),origin='lower')
#plt.show()
#plt.imshow((QSO_msk), origin='low') 
#plt.show()
#plt.imshow((QSO_msk_c), origin='low') 
#plt.show()
#pyfits.PrimaryHDU(QSO_msk*1).writeto('{0}_msk.fits'.format(ID),overwrite=True)

from photutils import detect_threshold
threshold = detect_threshold(QSO_img, snr=2.)
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources,deblend_sources
sigma = 3.0 * gaussian_fwhm_to_sigma    # FWHM = 3.
kernel = Gaussian2DKernel(sigma, x_size=5, y_size=5)
kernel.normalize()
segm = detect_sources(QSO_img, threshold, npixels=10, filter_kernel=kernel)
npixels = 20
segm_deblend = deblend_sources(QSO_img, segm, npixels=npixels,
                                filter_kernel=kernel, nlevels=32,
                                contrast=0.001)
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
norm = ImageNormalize(stretch=SqrtStretch())
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
ax1.imshow(QSO_img, origin='lower', norm=LogNorm())
ax1.set_title('Data')
ax2.imshow(segm_deblend, origin='lower', cmap=segm_deblend.cmap(random_state=12345))
ax2.set_title('Segmentation Image')
plt.show()

from photutils import source_properties
columns = ['id', 'xcentroid', 'ycentroid', 'source_sum', 'area']
cat = source_properties(QSO_img, segm_deblend)
tbl = cat.to_table(columns=columns)
tbl['xcentroid'].info.format = '.2f'  # optional format
tbl['ycentroid'].info.format = '.2f'
print(tbl)

from photutils import source_properties, EllipticalAperture
cat = source_properties(QSO_img, segm_deblend)
segm_deblend_size = segm_deblend.areas
apertures = []
for obj in cat:
    size = segm_deblend_size[obj.id]
    position = (obj.xcentroid.value, obj.ycentroid.value)
    a_o = obj.semimajor_axis_sigma.value
    b_o = obj.semiminor_axis_sigma.value
    size_o = np.pi * a_o * b_o
    r = np.sqrt(size/size_o)*1.2
    a, b = a_o*r, b_o*r
    theta = obj.orientation.value
    apertures.append(EllipticalAperture(position, a, b, theta=theta))

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
norm = ImageNormalize(stretch=SqrtStretch())
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
ax1.imshow(QSO_img, origin='lower', cmap='Greys_r', norm=norm)
ax1.set_title('Data')
ax2.imshow(segm_deblend, origin='lower',
           cmap=segm_deblend.cmap(random_state=12345))
ax2.set_title('Segmentation Image')
for aperture in apertures:
    aperture.plot(color='white', lw=1.5, ax=ax1)
    aperture.plot(color='white', lw=1.5, ax=ax2)

from regions import PixCoord, EllipsePixelRegion
from astropy.coordinates import Angle
center = PixCoord(x=37.14706524, y=130)
theta = Angle(aperture.theta/np.pi*180.,'deg')
reg = EllipsePixelRegion(center=center, width=aperture.a*2, height=aperture.b*2, angle=theta)
patch = reg.as_artist(facecolor='none', edgecolor='red', lw=2)
ax2.add_patch(patch)
plt.show()
mask_set = reg.to_mask(mode='center')
mask = mask_set.to_image((200,200))
plt.imshow(mask, origin='lower')
plt.show()

#labels = [1,2,3,4]
#cat = source_properties(QSO_img, segm_deblend, labels=labels)
#columns = ['id', 'xcentroid', 'ycentroid', 'source_sum', 'area']
#tbl3 = cat.to_table(columns=columns)
#tbl3['xcentroid'].info.format = '.4f'  # optional format
#tbl3['ycentroid'].info.format = '.4f'
#tbl3['source_sum'].info.format = '.4f'
#print(tbl3)

'''
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
#    pyfits.PrimaryHDU(PSF_msk*1).writeto('{0}_msk.fits'.format(PSF_i[:-5]),overwrite=True)
'''
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
