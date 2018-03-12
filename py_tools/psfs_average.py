#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 10:57:14 2018

@author: Dartoon

Doing the PSFs average given a list of PSF
"""

import numpy as np
from flux_profile import text_in_string_list,cr_mask
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

def psf_ave(psfs_list, not_count=(), mode = 'direct',  mask_list=['default.reg'],scale=3):
    '''
    Produce the average for a list of psfs.
    
    
    Parameter
    --------
        psfs_list:
            A raw list of psfs array. 
        mode: the way to do the average.
            'direcity':
                directly do the average of the scaled PSF 
            'CI':
                Consider Intensity of the PSF. (weighted by the root-square of the intensity),
                ---as the noise of PSF is most related to the Poission noise. The SNR of the image is related to root-square of the image.
        no_count:
            The serial No. of psf which not considered.
        mask_list: 
            A list include all the mask.reg names
    Return
    --------
        A averaged and scaled PSF.
    '''
    ### masked PSF give the region.
    psf_NO = len(psfs_list)
    
    psfs_masked_list = np.ones_like(psfs_list)  # To load the PSF plus masked area 
    mask_4ave = np.ones_like(psfs_list)
    for i in range(psf_NO):
        if i in not_count:
            print "The PSF{0} is not count".format(i)
            psfs_masked_list[i] = np.zeros_like(psfs_list[i])
            mask_4ave[i] *= 0
        else:
            msk_counts, mask_lists = text_in_string_list("PSF{0}".format(i), mask_list)
            mask = np.ones(psfs_list[i].shape)
            if msk_counts != 0:
                for j in range(msk_counts):
                    mask *= cr_mask(image=psfs_list[i], filename=mask_lists[j])
            psfs_masked_list[i] = psfs_list[i] * mask
            mask_4ave[i] = mask_4ave[i] * mask
#    for i in range(psf_NO):
#            print "plot psfs_list", i
#            plt.matshow(psfs_masked_list[i], origin= 'low', norm=LogNorm())
#            plt.colorbar()
#            plt.show()  
    ### Doing the average.
    if mode =='direct':
        for i in range(psf_NO):
            if psfs_masked_list[i].sum() != 0:
                psfs_masked_list[i] /= psfs_masked_list[i].sum()  # scale the image to a same level
        psf_total = np.sum(psfs_masked_list, axis=0)
        sum_4ave = np.sum(mask_4ave, axis=0)
        psf_final = psf_total/sum_4ave
    elif mode == 'CI':
        for i in range(psf_NO):
            if psfs_masked_list[i].sum() != 0:
                psfs_masked_list[i] /= np.sqrt(np.sum(psfs_masked_list[i]))  # scale the image based on their intensity (SNR)
        psf_total = np.sum(psfs_masked_list, axis=0)
        sum_4ave = np.sum(mask_4ave, axis=0)
        print sum_4ave
        psf_final = psf_total/sum_4ave
    #### The PSF are found not very necessary to be shiftted. !!!! Note the high_CI is not ready --- high_res. mask is not OK.
#    if mode == 'high_CI':
#        psfs_high_list = np.empty([psf_NO, psfs_list[0].shape[0]*scale, psfs_list[0].shape[1]*scale])
#        #Creat a mask_high_list:
#        mask_high_list = np.ones_like(psfs_high_list)
#        for i in range(psf_NO):
#            psfs_high_list[i] = im_2_high_res(psfs_list[i], scale=scale)  # scale the image to a same level
#            psfs_high_list[i] *= mask_high_list[i]
#            psfs_high_list[i] /= np.sqrt(np.sum(psfs_high_list[0]))
#        psf_high_total = np.sum(psfs_high_list, axis=0)
#        sum_4ave = np.sum(mask_high_list, axis=0)
#        print sum_4ave
#        psf_high_final = psf_high_total/sum_4ave
#        psf_final = rebin(psf_high_final, scale = scale)

    psf_final /= psf_final.sum()
    return psf_final
        
def im_2_high_res(img, scale=3, shift_center = True):
    '''
    To derive a higher image resolution
    
    Parameter
    --------
        img: 
            The array of the inpt image.
        scale:
            An int value, the factor for incresing the resolution.
        shift_center:
            If True, the brightest center will be at the center.
    Return
    --------
        A higher resolution image (array) of the input image.
    '''
    y,x = coords((int(len(img)*scale),int(len(img)*scale))) 
    center=np.where(img==img.max())
    psf_info=get_high_res(img=img,scale=scale)
    psf_info.nx, psf_info.ny= center[0][0]*scale, center[1][0]*scale  #need to be improved if want to get the center automaticlly???
    psf_high=psf_info.evaluateSource(x,y)
    if shift_center == True:
        center_o_high = np.reshape(np.asarray(np.where(psf_high== psf_high.max())),(2))[::-1]
        center_should = np.asarray(np.shape(psf_high))/2
        shift = center_should - center_o_high
        print "center shift:",shift
        psf_info.nx += shift[0]
        psf_info.ny += shift[1]
        psf_high=psf_info.evaluateSource(x,y)
        center_o_high = np.reshape(np.asarray(np.where(psf_high== psf_high.max())),(2))[::-1]
        if abs(center_should - center_o_high).sum() !=0 :
            raise ValueError("The high res. image is not shiftted to the center:",  "center_should:",center_should,"However:", center_o_high)
            print  "if",center_should,"?:", center_o_high
#        print  "if",center_should,"?:", center_o_high
#    if shift_center is True:
#        print True
    return psf_high

from scipy import ndimage
class get_high_res():
  def __init__(self,img=None,scale=None,order=1):  #Matt said the order is 1
      self.img = img #/np.sum(img)
      self.scale = scale
      self.dy,self.dx = img.shape
      self.x0 = self.dx/2
      self.y0 = self.dy/2
      self.order = order
      self.nx, self.ny = None, None   #For coordinate the pixel position
      if  self.order==1:
          self.model = self.img
      else:
          self.model = ndimage.spline_filter(self.img,output=np.float64,order=order)
  def getPix(self,sx,sy):
      x = (sx-self.nx)/self.scale+self.x0
      y = (sy-self.ny)/self.scale+self.y0
      return x,y
  def evaluateSource(self,sx,sy):
      ix,iy = self.getPix(sx,sy)
      return ndimage.map_coordinates(self.model,[[iy],[ix]],order=self.order)[0,:,:]/self.scale**2
      # Other function names, for convenience

def coords(shape):
    return np.indices(shape).astype(np.float64)
    
def rebin(image, scale=3):         
    '''
    Rebin a image to lower resolution 
    
    Parameter
    --------
        image:
        scale:
    Return
    --------
    '''
    shape = (len(image)/scale, len(image)/scale)
    sh = shape[0], scale, shape[1], scale
    return image.reshape(sh).mean(-1).mean(1)*scale**2

