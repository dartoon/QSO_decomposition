#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 12:03:26 2018

@author: dartoon
"""

import numpy as np
import os
filt = 'F140w'
pix_sz = 'drz06'

if filt == 'F140w':
    zp = 26.4524
elif filt == 'F125w':
    zp = 26.2303
elif filt == 'acs':
    zp = 25.94333  # The AB zp for F814w, after 2006/7/4 

if pix_sz == 'drz06':
    deltaPix = 0.0642
elif pix_sz == 'acs':
    deltaPix = 0.03

filename = 'galfit.01' 
fit_out = open('{0}'.format(filename),'r')
lines = fit_out.readlines()
sersic_re = float(lines[42][4:12]) * deltaPix
sersic_n = float(lines[43][4:11])

import astropy.io.fits as pyfits
gal_data = pyfits.open('imgblock_QSO.fits')[1].data.copy()
gal_bestfit = pyfits.open('imgblock_QSO.fits')[2].data.copy()
gal_residual = pyfits.open('imgblock_QSO.fits')[3].data.copy()

sersic_mag = -2.5 * np.log10(gal_bestfit.sum()) + zp

noise_map = pyfits.getdata('noise_level.fits')
chiq_map = (gal_residual/noise_map)**2
pixels=len(noise_map)**2
reduced_Chisq = chiq_map.sum()/pixels 

lenstronomy_redisual = pyfits.open('plan_b_residual.fits')[0].data.copy()
print  "{0}\t{1}\t{2}\t{3}".format(sersic_mag, sersic_n, sersic_re,
        reduced_Chisq)
import copy
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import matplotlib
f = plt.figure(0, figsize=(22.333,5.333))
ax1 = plt.subplot2grid((3,3), (0,0), rowspan=3)
ax2 = plt.subplot2grid((3,3), (0,1), rowspan=3)
ax3 = plt.subplot2grid((3,3), (0,2), rowspan=3)
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')
c_ax1 = ax1.imshow(gal_data,origin='lower',cmap=my_cmap, norm=LogNorm(), vmax = gal_data.max() )
clim=c_ax1.properties()['clim']
frame_size = len(gal_data[0])
ax1.text(frame_size*0.05, frame_size*0.9, 'Input to Galfit', fontsize=15,color='w',backgroundcolor='k')
ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)
ax2.imshow(gal_bestfit,origin='lower',cmap=my_cmap, norm=LogNorm(),  clim=clim)
pos2_o = ax2.get_position() # get the original position
pos2 = [pos2_o.x0 -0.07, pos2_o.y0, pos2_o.width, pos2_o.height]
ax2.set_position(pos2) # set a new position
ax2.text(frame_size*0.05, frame_size*0.9, 'Galfit best-fit', fontsize=15,color='w',backgroundcolor='k')
ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)
ax3.imshow((gal_data-gal_bestfit),origin='lower',cmap=my_cmap, norm=LogNorm(),  clim=clim)
ax3.text(frame_size*0.05, frame_size*0.9, 'Residual', fontsize=15,color='w',backgroundcolor='k')
pos3_o = ax3.get_position() # get the original position
pos3 = [pos3_o.x0 -0.14, pos3_o.y0, pos3_o.width, pos3_o.height]
ax3.set_position(pos3) # set a new position
ax3.get_xaxis().set_visible(False)
ax3.get_yaxis().set_visible(False)
#plt.colorbar(c_ax1)