#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 14:28:33 2018

@author: dartoon

A test of the region file
"""
import numpy as np
from regions import PixCoord, CirclePixelRegion 
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
reg_x, reg_y=59, 49         #One must note that python starts from 0!!!
center= PixCoord(x=reg_x,y=reg_y)
radius=1
region = CirclePixelRegion(center, radius)

fitsFile = pyfits.open('psf.fits')
img = fitsFile[0].data 
#img[np.isnan(img)] = 0

#cutout=img[(reg_x-radius):(reg_x+radius+1), (reg_y-radius):(reg_y+radius+1)]
mask = region.to_mask(mode='exact')
data = mask.cutout(img)

#from regions import DS9Parser
#parser = DS9Parser(region)
#print parser

#from regions import write_ds9, read_ds9
#filename = 'test_circle.reg'
#write_ds9(region, filename)
#region_test=read_ds9(filename)
#region_=ds9_objects_to_string(region_test)
#print region_

#ax=plt.subplot(1,1,1)
#cax=ax.imshow((img),origin='lower')#,vmin=0,vmax=1)
##ax.add_patch(mask.bbox.as_patch(facecolor='none', edgecolor='white'))
##ax.add_patch(region.as_patch(facecolor='none', edgecolor='orange'))
##ax.add_patch(mask.bbox.as_patch(facecolor='none', edgecolor='white'))
#ax.add_patch(region_test.as_patch(facecolor='none', edgecolor='orange'))
#plt.colorbar(cax)
#
