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
from matplotlib.colors import LogNorm

def cut_image(image, center, radius):
    region = pix_region(center, radius=radius)
    cut = region.to_mask(mode='exact')
    cut_image = cut.cutout(image)
    return cut_image

def cut_center_bright(image, center, radius, return_center=False):
    """
    Auto cut the image, with with brightest center in the center.
    """
    temp_center = np.asarray(center)
    radius = radius
    img_test = cut_image(image=image, center=temp_center, radius=radius)
    frm_qrt = len(img_test)/4
    test_center =  np.asarray(np.where(img_test == img_test[frm_qrt:-frm_qrt,frm_qrt:-frm_qrt].max()))[:,0]
    center_shift = np.array((test_center- radius))[::-1]
    center = (temp_center + center_shift)
    cut_c_b = cut_image(image=image, center=center, radius=radius)
    if return_center==False:
        return cut_c_b
    if return_center==True:
        return cut_c_b, center
 
def save_loc_png(img, center_QSO, c_psf_list=None,extra_psfs=None,ID=None,
                 reg_ty=None,ifsave=True, label_shift_NO=(), shift_where=None):
    '''
    label shift_where: 1,2,3,4 --- up, right, down, left
    '''
    fig = plt.figure(figsize=(15,15))
    ax=fig.add_subplot(1,1,1)
    import copy, matplotlib
    my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
    my_cmap.set_bad('black')
    if reg_ty == None:
        vmax = 2.2
        vmin = 1.e-2
        QSO_box_size = 30
        PSF_box_size = 20
    elif reg_ty == 'astrodrz_06':
        vmax = 2.1 
        vmin = 1.e-3
        QSO_box_size = 60
        PSF_box_size = 40
    elif reg_ty == 'acs':
        vmax = 2.1 
        vmin = 1.e-3
        QSO_box_size = 100
        PSF_box_size = 60 
    cax=ax.imshow(img,origin='lower', cmap=my_cmap, norm=LogNorm(), vmin=vmin, vmax=vmax)
    QSO_reg = pix_region(center_QSO, radius= QSO_box_size)
    QSO_mask = QSO_reg.to_mask(mode='center')
    ax.text(center_QSO[0]-2*QSO_box_size, center_QSO[1]+1.5*QSO_box_size, 'QSO',color='white', fontsize=20)
    ax.add_patch(QSO_mask.bbox.as_patch(facecolor='none', edgecolor='white', linewidth=2))
    count=0
    count_shift = 0
    if c_psf_list is not None:
        for i in range(len(c_psf_list)):
            PSF_reg = pix_region(c_psf_list[i], radius= PSF_box_size)
            PSF_mask = PSF_reg.to_mask(mode='center')
            ax.add_patch(PSF_mask.bbox.as_patch(facecolor='none', edgecolor='blue', linewidth=2))
            if count not in label_shift_NO:
                ax.text(c_psf_list[i][0]-2*PSF_box_size, c_psf_list[i][1]+2*PSF_box_size, 'PSF{0}'.format(count),color='white', fontsize=15)
            else:
                if count in label_shift_NO:
                    shift = shift_label_index(shift_where[count_shift])
                    ax.text(c_psf_list[i][0]+shift[0]*PSF_box_size, c_psf_list[i][1]+shift[1]*PSF_box_size,
                            'PSF{0}'.format(count),color='white', fontsize=15)
                    count_shift += 1
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            count += 1
#    plt.colorbar(cax)
    if extra_psfs is not None:
        for i in range(len(extra_psfs)):
            PSF_reg = pix_region(extra_psfs[i], radius= PSF_box_size)
            PSF_mask = PSF_reg.to_mask(mode='center')
            ax.add_patch(PSF_mask.bbox.as_patch(facecolor='none', edgecolor='yellow', linewidth=2))
            if count not in label_shift_NO:
                ax.text(extra_psfs[i][0]-2*PSF_box_size, extra_psfs[i][1]+2*PSF_box_size, 'PSF{0}?'.format(count),color='white', fontsize=15)
            else:
                if count in label_shift_NO:
                    shift = shift_label_index(shift_where[count_shift])
                    print extra_psfs[i][0], shift[0]*PSF_box_size
                    ax.text(extra_psfs[i][0]+shift[0]*PSF_box_size, extra_psfs[i][1]+shift[1]*PSF_box_size,
                            'PSF{0}?'.format(count),color='white', fontsize=15)
                    count_shift += 1
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            count += 1
#    plt.colorbar(cax)
    if not ID == None:
        ax.text(len(img)*0.05, len(img)*0.8, ID,color='white', fontsize=30)
    if ifsave == True:
        fig.savefig('QSO_PSF_loc.pdf')

def shift_label_index(s):
    if s==1:
        return np.array([-2, 2])
    elif s==2:
        return np.array([1, -1])
    elif s==3:
        return np.array([-2, -3])
    elif s==4:
        return np.array([-6, -1])    
    
def save_other_target(img, center_QSO, other_target,target_name,c_psf_list,extra_psfs,ID=None,reg_ty=None,ifsave=True):
    fig = plt.figure(figsize=(15,15))
    ax=fig.add_subplot(1,1,1)
    import copy, matplotlib
    my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
    my_cmap.set_bad('black')
    if reg_ty == None:
        vmax = 2.2
        vmin = 1.e-2
        QSO_box_size = 30
        PSF_box_size = 20
    elif reg_ty == 'astrodrz_06':
        vmax = 2.1 
        vmin = 1.e-3
        QSO_box_size = 60
        PSF_box_size = 40
    elif reg_ty == 'acs':
        vmax = 2.1 
        vmin = 1.e-3
        QSO_box_size = 100
        PSF_box_size = 60 
    cax=ax.imshow(img,origin='lower', cmap=my_cmap, norm=LogNorm(), vmin=vmin, vmax=vmax)
    QSO_reg = pix_region(center_QSO, radius= QSO_box_size)
    QSO_mask = QSO_reg.to_mask(mode='center')
    ax.text(center_QSO[0]-2*QSO_box_size, center_QSO[1]+1.5*QSO_box_size, 'QSO',color='white', fontsize=20)
    ax.add_patch(QSO_mask.bbox.as_patch(facecolor='none', edgecolor='white', linewidth=2))
#   setting for other target
    target_reg = pix_region(other_target, radius= PSF_box_size)
    target_mask = target_reg.to_mask(mode='center')
    ax.add_patch(target_mask.bbox.as_patch(facecolor='none', edgecolor='white', linewidth=2))
    ax.text(other_target[0]-5*PSF_box_size, other_target[1]+2*PSF_box_size, '{0}'.format(target_name),color='white', fontsize=15)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
#    plt.colorbar(cax)
    count = 0
    if c_psf_list is not None:
        for i in range(len(c_psf_list)):
            PSF_reg = pix_region(c_psf_list[i], radius= PSF_box_size)
            PSF_mask = PSF_reg.to_mask(mode='center')
            ax.add_patch(PSF_mask.bbox.as_patch(facecolor='none', edgecolor='blue', linewidth=2))
            ax.text(c_psf_list[i][0]-2*PSF_box_size, c_psf_list[i][1]+2*PSF_box_size, 'PSF{0}'.format(count),color='white', fontsize=15)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            count += 1
    if extra_psfs is not None:
        for i in range(len(extra_psfs)):
            PSF_reg = pix_region(extra_psfs[i], radius= PSF_box_size)
            PSF_mask = PSF_reg.to_mask(mode='center')
            ax.add_patch(PSF_mask.bbox.as_patch(facecolor='none', edgecolor='yellow', linewidth=2))
            ax.text(extra_psfs[i][0]-2*PSF_box_size, extra_psfs[i][1]+2*PSF_box_size, 'PSF{0}?'.format(count),color='white', fontsize=15)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            count += 1
#    plt.colorbar(cax)
    if not ID == None:
        ax.text(len(img)*0.05, len(img)*0.8, ID,color='white', fontsize=30)
    if ifsave == True:
        fig.savefig('QSO_PSF_loc.pdf')
    
def QSO_star_mag(img, center_QSO, QSO_mag, psf_list,mag, ID=None, reg_ty='astrodrz_06', ifsave=False):
    fig = plt.figure(figsize=(15,15))
    ax=fig.add_subplot(1,1,1)
    import copy, matplotlib
    my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
    my_cmap.set_bad('black')
    if reg_ty == None:
        vmax = 2.2
        vmin = 1.e-2
        QSO_box_size = 30
        PSF_box_size = 20
    elif reg_ty == 'astrodrz_06':
        vmax = 2.1 
        vmin = 1.e-3
        QSO_box_size = 60
        PSF_box_size = 40
    elif reg_ty == 'acs':
        vmax = 2.1 
        vmin = 1.e-3
        QSO_box_size = 100
        PSF_box_size = 60 
    ax.imshow(img,origin='lower', cmap=my_cmap, norm=LogNorm(), vmin=vmin, vmax=vmax)
    QSO_reg = pix_region(center_QSO, radius= QSO_box_size)
    QSO_mask = QSO_reg.to_mask(mode='center')
    ax.text(center_QSO[0]-3*QSO_box_size, center_QSO[1]+1.5*QSO_box_size, 'QSO:',color='white', fontsize=15)
    ax.text(center_QSO[0]+0*QSO_box_size, center_QSO[1]+1.5*QSO_box_size, '{0}'.format(round(QSO_mag,2)),color='white', fontsize=15)
    ax.add_patch(QSO_mask.bbox.as_patch(facecolor='none', edgecolor='white', linewidth=2))
    count=0
    for i in range(len(psf_list)):
        PSF_reg = pix_region(psf_list[i], radius= PSF_box_size)
        PSF_mask = PSF_reg.to_mask(mode='center')
        ax.add_patch(PSF_mask.bbox.as_patch(facecolor='none', edgecolor='blue', linewidth=2))
        ax.text(psf_list[i][0]-5*PSF_box_size, psf_list[i][1]+2*PSF_box_size, 'PSF{0}:'.format(count),color='white', fontsize=13)
        ax.text(psf_list[i][0]-1*PSF_box_size, psf_list[i][1]+2*PSF_box_size, '{0}'.format(round(mag[i],2)),color='c', fontsize=13)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        count += 1
    if not ID == None:
        ax.text(len(img)*0.05, len(img)*0.8, ID,color='white', fontsize=30)
    if ifsave == True:
        fig.savefig('PSF_mag.pdf')
        
def QSO_star_color(img, QSO_pos, QSO_mags, psf_list, mag_0, mag_1,
                   mag_diff, ID, reg_ty='astrodrz_06', ifsave=False,
                   label_shift_NO=(),shift_where=None):
    fig = plt.figure(figsize=(15,15))
    ax=fig.add_subplot(1,1,1)
    import copy, matplotlib
    my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
    my_cmap.set_bad('black')
    if reg_ty == 'astrodrz_06':
        vmax = 2.1 
        vmin = 1.e-3
        QSO_box_size = 40
        PSF_box_size = 40
    elif reg_ty == 'acs':
        vmax = 2.1 
        vmin = 1.e-3
        QSO_box_size = 60
        PSF_box_size = 60 
    ax.imshow(img,origin='lower', cmap=my_cmap, norm=LogNorm(), vmin=vmin, vmax=vmax)
    QSO_reg = pix_region(QSO_pos, radius= QSO_box_size)
    QSO_mask = QSO_reg.to_mask(mode='center')
    ax.text(QSO_pos[0]-4*QSO_box_size, QSO_pos[1]+2.6*QSO_box_size, 'QSO',color='white', fontsize=13)
    ax.text(QSO_pos[0]-1*QSO_box_size, QSO_pos[1]+3.6*QSO_box_size, '{0}'.format(round(QSO_mags[0],2)),color='yellow', fontsize=13)
    ax.text(QSO_pos[0]-1*QSO_box_size, QSO_pos[1]+2.6*QSO_box_size, '{0}'.format(round(QSO_mags[1],2)),color='c', fontsize=13)
    ax.text(QSO_pos[0]-1*QSO_box_size, QSO_pos[1]+1.6*QSO_box_size, '{0}'.format(round(QSO_mags[0]-QSO_mags[1],2)),color='white', fontsize=13)
    ax.add_patch(QSO_mask.bbox.as_patch(facecolor='none', edgecolor='white', linewidth=2))
    count=0
    count_shift = 0
    for i in range(len(psf_list)):
        PSF_reg = pix_region(psf_list[i], radius= PSF_box_size)
        PSF_mask = PSF_reg.to_mask(mode='center')
        if count not in label_shift_NO:
            ax.add_patch(PSF_mask.bbox.as_patch(facecolor='none', edgecolor='blue', linewidth=2))
            ax.text(psf_list[i][0]-4*PSF_box_size, psf_list[i][1]+2*PSF_box_size, 'PSF{0}'.format(count),color='white', fontsize=13)
            ax.text(psf_list[i][0]-1*PSF_box_size, psf_list[i][1]+3*PSF_box_size, '{0}'.format(round(mag_0[i],2)),color='yellow', fontsize=13)
            ax.text(psf_list[i][0]-1*PSF_box_size, psf_list[i][1]+2*PSF_box_size, '{0}'.format(round(mag_1[i],2)),color='c', fontsize=13)
            ax.text(psf_list[i][0]-1*PSF_box_size, psf_list[i][1]+1*PSF_box_size, '{0}'.format(round(mag_diff[i],2)),color='white', fontsize=13)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
        elif count in label_shift_NO:
            shift = shift_label_index(shift_where[count_shift])        
            ax.add_patch(PSF_mask.bbox.as_patch(facecolor='none', edgecolor='blue', linewidth=2))
            ax.text(psf_list[i][0]+shift[0]*PSF_box_size, psf_list[i][1]+shift[1]*PSF_box_size, 'PSF{0}'.format(count),color='white', fontsize=13)
            ax.text(psf_list[i][0]+(shift[0]+3)*PSF_box_size, psf_list[i][1]+(shift[1]+1)*PSF_box_size, '{0}'.format(round(mag_0[i],2)),color='yellow', fontsize=13)
            ax.text(psf_list[i][0]+(shift[0]+3)*PSF_box_size, psf_list[i][1]+(shift[1]+0)*PSF_box_size, '{0}'.format(round(mag_1[i],2)),color='c', fontsize=13)
            ax.text(psf_list[i][0]+(shift[0]+3)*PSF_box_size, psf_list[i][1]+(shift[1]-1)*PSF_box_size, '{0}'.format(round(mag_diff[i],2)),color='white', fontsize=13)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            count_shift += 1
        count += 1
    if not ID == None:
        ax.text(len(img)*0.05, len(img)*0.8, ID,color='white', fontsize=30)
    if ifsave == True:
        fig.savefig('PSF_color.pdf')
    
def grab_pos(filename, reg_ty='swarp', QSO_reg_return = False):
    '''
    Grab all the positions of the QSO and stars from a region name. The last one are always the QSO.
    '''
    with open(filename, 'r') as input_file:
        reg_string=input_file.read()
    string_list=reg_string.split('\n')
    pos_string = []
    for i in range(len(string_list)):
        if reg_ty == 'swarp':
            string = string_find_between(string_list[i],"(", ",78.134039")
        elif reg_ty == 'org':
            string = string_find_between(string_list[i],"(", ",77.972708")
        elif reg_ty == 'astrodrz_06':
            string = string_find_between(string_list[i],"(", ",155.76324")
        elif reg_ty == 'acs':
            string = string_find_between(string_list[i],"(", ",333.33335")        
        if string.split(',')[0] != '':
            pos_list = [float(j) for j in string.split(',')]
            pos_string.append(pos_list)
    if QSO_reg_return == True:
        for i in range(1,len(string_list)):
            if 'red' in string_list[i]:
                count_QSO = i
            if '(' in string_list[i] and '(' not in string_list[i-1]:
                count_reg = i
        QSO_loc = count_QSO - count_reg
        return np.asarray(pos_string), QSO_loc
    elif QSO_reg_return == False:
        return np.asarray(pos_string)
    
def string_find_between(s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""
