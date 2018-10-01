#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 15:01:17 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from photutils import detect_threshold
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources,deblend_sources
from matplotlib.colors import LogNorm


def mask_obj(img, snr=2.5, exp_sz= 1.2):
    threshold = detect_threshold(img, snr=snr)
    center_img = len(img)/2
    sigma = 3.0 * gaussian_fwhm_to_sigma# FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=5, y_size=5)
    kernel.normalize()
    segm = detect_sources(img, threshold, npixels=10, filter_kernel=kernel)
    npixels = 20
    segm_deblend = deblend_sources(img, segm, npixels=npixels,
                                    filter_kernel=kernel, nlevels=32,
                                    contrast=0.001)
    #Number of objects segm_deblend.data.max()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 10))
    import copy, matplotlib
    my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
    my_cmap.set_bad('black')
    vmin = 1.e-3
    vmax = 2.1 
    ax1.imshow(img, origin='lower', cmap=my_cmap, norm=LogNorm(), vmin=vmin, vmax=vmax)
    ax1.set_title('Data')
    ax2.imshow(segm_deblend, origin='lower', cmap=segm_deblend.cmap(random_state=12345))
    ax2.set_title('Segmentation Image')
    plt.close()
    
    from photutils import source_properties
    columns = ['id', 'xcentroid', 'ycentroid', 'source_sum', 'area']
    cat = source_properties(img, segm_deblend)
    tbl = cat.to_table(columns=columns)
    tbl['xcentroid'].info.format = '.2f'  # optional format
    tbl['ycentroid'].info.format = '.2f'
    print(tbl)
    
    from photutils import source_properties, EllipticalAperture
    cat = source_properties(img, segm_deblend)
    segm_deblend_size = segm_deblend.areas
    apertures = []
    for obj in cat:
        size = segm_deblend_size[obj.id]
        position = (obj.xcentroid.value, obj.ycentroid.value)
        a_o = obj.semimajor_axis_sigma.value
        b_o = obj.semiminor_axis_sigma.value
        size_o = np.pi * a_o * b_o
        r = np.sqrt(size/size_o)*exp_sz
        a, b = a_o*r, b_o*r
        theta = obj.orientation.value
        apertures.append(EllipticalAperture(position, a, b, theta=theta))
    
    dis_sq = [np.sqrt((apertures[i].positions[0][0]-center_img)**2+(apertures[i].positions[0][1]-center_img)**2) for i in range(len(apertures))]
    dis_sq = np.asarray(dis_sq)
    c_index= np.where(dis_sq == dis_sq.min())[0][0]
    #from astropy.visualization.mpl_normalize import ImageNormalize
    #norm = ImageNormalize(stretch=SqrtStretch())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 10))
    ax1.imshow(img, origin='lower', cmap=my_cmap, norm=LogNorm(), vmin=vmin, vmax=vmax)
    ax1.set_title('Data')
    ax2.imshow(segm_deblend, origin='lower',
               cmap=segm_deblend.cmap(random_state=12345))
    ax2.set_title('Segmentation Image')
    for i in range(len(apertures)):
        aperture = apertures[i]
        aperture.plot(color='white', lw=1.5, ax=ax1)
        aperture.plot(color='white', lw=1.5, ax=ax2)
    plt.show()
    
    from regions import PixCoord, EllipsePixelRegion
    from astropy.coordinates import Angle
    
    obj_masks = []
    for i in range(len(apertures)):
        aperture = apertures[i]
        x, y = aperture.positions[0]
        center = PixCoord(x=x, y=y)
        theta = Angle(aperture.theta/np.pi*180.,'deg')
        reg = EllipsePixelRegion(center=center, width=aperture.a*2, height=aperture.b*2, angle=theta)
        patch = reg.as_artist(facecolor='none', edgecolor='red', lw=2)
        fig, axi = plt.subplots(1, 1, figsize=(10, 12.5))
        axi.add_patch(patch)
        mask_set = reg.to_mask(mode='center')
        mask = mask_set.to_image((len(img),len(img)))
        axi.imshow(mask, origin='lower')
        plt.close()
        if i != c_index:
            obj_masks.append(mask)
        elif i == c_index:
            target_mask = mask
    return target_mask, obj_masks