#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 11:20:18 2018

@author: dartoon

Propose:
    derive the SB/tot_flux profile when give a fit.
"""

import numpy as np
from regions import PixCoord, CirclePixelRegion 
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

def pix_region(center=([49,49]), radius=5):
    '''
    Creat a region file, with pixel units
    
    Parameter
    --------
        center: The center of the region, with ([reg_x, reg_y]);
        radius: The radius of the region.
        
    Return
    --------
        A region which is ds9-like.
    '''
    center= PixCoord(x=center[0],y=center[1])
    region = CirclePixelRegion(center, radius)
    return region

def flux_in_region(image,region,mode='exact'):
    '''
    Calculate the total flux inside a given region.
    
    Parameter
    --------
        image: 2-D array image;
        region: The region generated by pix_region;
        mode: mode type of the mask, 'exact', 'center', default is 'exact'.
        
    Returns
    --------
        Total flux
    '''
    mask = region.to_mask(mode=mode)
    data = mask.cutout(img)
    tot_flux= np.sum(mask.data * data)
    return tot_flux

def flux_profile(image, center, radius=35, grids=20, ifplot=True, fits_plot=True):
    '''
    Derive the flux profile of one image start at the center.
    
    Parameters
    --------
        image: A 2-D array image;
        center: The center point of the profile;
        radius: The radius of the profile favourable with default equals to 35;
        grids: The number of points to sample the flux with default equals to 20;
        ifplot: if plot the profile
        fits_plot: if plot the fits file with the regions.
    Returns
    --------
        A 1-D array of the tot_flux value of each 'grids' in the profile sampled radius.    
    '''
    r_grids=(np.linspace(0,1,grids+1)*radius)[1:]
    r_flux = np.empty(grids)
    regions = []
    for i in range(len(r_grids)):
        region = pix_region(center, r_grids[i])
        r_flux[i] =flux_in_region(image, region)
        regions.append(region)
    if fits_plot == True:
        ax=plt.subplot(1,1,1)
        cax=ax.imshow((image),origin='lower')#,vmin=0,vmax=1)
        #ax.add_patch(mask.bbox.as_patch(facecolor='none', edgecolor='white'))
        for i in range(grids):
            ax.add_patch(regions[i].as_patch(facecolor='none', edgecolor='orange'))
        plt.colorbar(cax)
    if ifplot == True:
        from matplotlib.ticker import AutoMinorLocator
        minorLocator = AutoMinorLocator()
        fig, ax = plt.subplots()
        plt.plot(r_grids, r_flux, 'x-')
        ax.xaxis.set_minor_locator(minorLocator)
        plt.tick_params(which='both', width=2)
        plt.tick_params(which='major', length=7)
        plt.tick_params(which='minor', length=4, color='r')
        plt.grid()
        ax.set_ylabel("Total Flux")
        ax.set_xlabel("Pixels")
        plt.grid(which="minor")
        plt.show()
    return r_flux, r_grids, regions

def SB_profile(image, center, radius=35, grids=20, ifplot=True, fits_plot=True):
    '''
    Derive the SB profile of one image start at the center.
    
    Parameters
    --------
        image: A 2-D array image;
        center: The center point of the profile;
        radius: The radius of the profile favourable with default equals to 35;
        grids: The number of points to sample the flux with default equals to 20;
        ifplot: if plot the profile
        fits_plot: if plot the fits file with the regions.
    Returns
    --------
        A 1-D array of the SB value of each 'grids' in the profile with the sampled radius.
    '''
    r_flux, r_grids, regions=flux_profile(image, center, radius=radius, grids=grids, ifplot=False, fits_plot=fits_plot)
    region_area = np.zeros([len(r_flux)])
    for i in range(len(r_flux)):
        circle=regions[i].to_mask(mode='exact')
        region_area[i]=circle.data.sum()
    r_SB= r_flux/region_area
    if ifplot == True:
        from matplotlib.ticker import AutoMinorLocator
        minorLocator = AutoMinorLocator()
        fig, ax = plt.subplots()
        plt.plot(r_grids, r_SB, 'x-')
        ax.xaxis.set_minor_locator(minorLocator)
        plt.tick_params(which='both', width=2)
        plt.tick_params(which='major', length=7)
        plt.tick_params(which='minor', length=4, color='r')
        plt.grid()
        ax.set_ylabel("Total Flux")
        ax.set_xlabel("Pixels")
        plt.grid(which="minor")
        plt.show()
    return r_SB, r_grids

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

def mask(filename='test_circle.reg', image = np.ones([99,99])):
    '''
    The creat a mask with a .reg file. The pixels in the region is 0, otherwise 1.
    
    Parameter
    --------
        filename: filename of the .reg
        image:
        
    Return
    --------
        A image.shape array. Pixels in the region is 0, otherwise 1.
    '''
    ##### Note the numpy starts from 0, especially in the center,
    ####Need to check the center, shape of the boxtype, the size?
    with open(filename, 'r') as input_file:
        reg_string=input_file.read().replace('\n', '')
    print reg_string
    if "physicalcircle" in reg_string:
        abc=find_between(reg_string, "(", ")")
        reg_info=np.fromstring(abc, dtype=float, sep=',')
        center, radius = reg_info[:2]-1, reg_info[2]
        region = pix_region(center, radius)
        box = 1-region.to_mask(mode='center').data
    elif "physicalbox" in reg_string:
        abc=find_between(reg_string, "(", ")")
        reg_info=np.fromstring(abc, dtype=float, sep=',')
        center = reg_info[:2] - 1
        x_r, y_r = reg_info[2:4]
        box = np.zeros([np.int(y_r), np.int(x_r)])  #!!! need to test x and y position
    else:
        raise ValueError("The input reg is un-defined yet")
    frame_size = image.shape
    box_size = box.shape
    a = np.int(center[0]-box_size[0]/2)
    b = a + box_size[0]
    c = np.int(center[1]-box_size[1]/2)
    d = c + box_size[1]
    mask = np.ones(frame_size)
    mask_box_part = mask[a:b,c:d]
    mask_box_part *= box
    return mask
        
    
# =============================================================================
# Example:
# =============================================================================
fitsFile = pyfits.open('psf.fits')
img = fitsFile[0].data 
#region = pix_region(center=([49,49]), radius=5)
#flux_in_region(img, region)
#SB_profile(img, center=([49,49]),ifplot=False, fits_plot=False)

mask = mask(filename='test_circle.reg', image=img)
plt.imshow((mask),origin='lower')
plt.show()