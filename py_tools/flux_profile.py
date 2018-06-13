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
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import LogNorm
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

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
    data = mask.cutout(image)
    tot_flux= np.sum(mask.data * data)
    return tot_flux

def flux_profile(image, center, radius=35, grids=20, gridspace=None, ifplot=False, fits_plot=False):
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
        1. A 1-D array of the tot_flux value of each 'grids' in the profile sampled radius. 
        2. The grids of each pixel radius.
        3. The region file for each radius.
    '''
    if gridspace == None:
        r_grids=(np.linspace(0,1,grids+1)*radius)[1:]
        diff = 0.5 - r_grids[0]
        r_grids += diff             #starts from pixel 0.5
    elif gridspace == 'log':
        r_grids=(np.logspace(-2,0,grids+1)*radius)[1:]
        diff = 0.5 - r_grids[0]
        r_grids += diff             #starts from pixel 0.5
    r_flux = np.empty(grids)
    regions = []
    for i in range(len(r_grids)):
        region = pix_region(center, r_grids[i])
        r_flux[i] =flux_in_region(image, region)
        regions.append(region)
    if fits_plot == True:
        ax=plt.subplot(1,1,1)
        cax=ax.imshow((image),origin='lower',cmap='gist_heat')
        #ax.add_patch(mask.bbox.as_patch(facecolor='none', edgecolor='white'))
        for i in range(grids):
            ax.add_patch(regions[i].as_patch(facecolor='none', edgecolor='orange'))
        plt.colorbar(cax)
        plt.show()
    if ifplot == True:
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
        if gridspace == 'log':
            ax.set_xscale('log')
            plt.xlim(0.5, ) 
        plt.grid(which="minor")
        plt.show()
    return r_flux, r_grids, regions

def SB_profile(image, center, radius=35, grids=20, gridspace = None, 
               ifplot=False, fits_plot=False,
               mask_list=None, mask_plot = False,
               mask_cut = 0., if_annuli= False, msk_image=None):
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
        msk_image: if not is not one, will use this image as mask.
    Returns
    --------
        A 1-D array of the SB value of each 'grids' in the profile with the sampled radius.
    '''
    if mask_list == None:
        r_flux, r_grids, regions=flux_profile(image, center, radius=radius, grids=grids, gridspace=gridspace, ifplot=False, fits_plot=False)
        region_area = np.zeros([len(r_flux)])
        for i in range(len(r_flux)):
            circle=regions[i].to_mask(mode='exact')
            edge_mask = circle.cutout(np.ones(image.shape))
            region_area[i]=(circle.data * edge_mask).sum()
    elif mask_list != None:
        mask = np.ones(image.shape)
        if msk_image is None:
            for i in range(len(mask_list)):
                mask *= cr_mask(image=image, filename=mask_list[i], mask_reg_cut = mask_cut)
        else:
            mask = msk_image
#        plt.imshow(mask, origin='low')
#        plt.show()
        r_flux, r_grids, regions=flux_profile(image*mask, center, radius=radius, grids=grids, gridspace=gridspace, ifplot=False, fits_plot=False)
        region_area = np.zeros([len(r_flux)])
        for i in range(len(r_flux)):
            circle=regions[i].to_mask(mode='exact')
            circle_mask =  circle.cutout(mask)
            if i ==len(r_flux)-1 and mask_plot == True:
#                print "plt circle_mask"
                plt.imshow((circle_mask),origin='lower')
                plt.show()
#                plt.imshow((circle.data),origin='lower') #The circle data is zero at outer circle area is zero
#                plt.show()
            region_area[i]=(circle.data * circle_mask).sum()
    if if_annuli ==False:
        r_SB= r_flux/region_area
    elif if_annuli == True:
        r_SB = np.zeros_like(r_flux)
        r_SB[0] = r_flux[0]/region_area[0]
        r_SB[1:] = (r_flux[1:]-r_flux[:-1]) / (region_area[1:]-region_area[:-1])
    if fits_plot == True:
        ax=plt.subplot(1,1,1)
        if mask_list != None:
            cax=ax.imshow(image*mask, norm=LogNorm(),origin='lower')
        elif mask_list == None:
            cax=ax.imshow(image,  norm=LogNorm(),origin='lower')
        #ax.add_patch(mask.bbox.as_patch(facecolor='none', edgecolor='white'))
        for i in range(grids):
            ax.add_patch(regions[i].as_patch(facecolor='none', edgecolor='orange'))
        plt.colorbar(cax)
        plt.show()
    if ifplot == True:
        minorLocator = AutoMinorLocator()
        fig, ax = plt.subplots()
        plt.plot(r_grids, r_SB, 'x-')
        ax.xaxis.set_minor_locator(minorLocator)
        plt.tick_params(which='both', width=2)
        plt.tick_params(which='major', length=7)
        plt.tick_params(which='minor', length=4, color='r')
        plt.grid()
        ax.set_ylabel("Surface Brightness")
        ax.set_xlabel("Pixels")
        if gridspace == 'log':
            ax.set_xscale('log')
            plt.xlim(0.5, ) 
        plt.grid(which="minor")
        plt.show()
    return r_SB, r_grids

def text_in_string_list(text, string_list):
    '''
    Parameter
    --------
        text: a string that could existed in the string_list
        string_list: The list of name
    Return
    --------
        counts:
            how many items have the text,
            The remainds in the string_list which have text.
    '''
    counts = 0
    text_string=[]
    for i in range(len(string_list)):
        if text in string_list[i]:
            counts += 1
            text_string.append(string_list[i])
    return counts, text_string
            

def QSO_psfs_compare(QSO, psfs, mask_list=None, plt_which_PSF=None,
                     include_QSO = True, gridspace = None , grids=30, norm_pix = 3.0,
                     if_annuli=False,plt_QSO=False,astrodrz=False, not_plt=()):
    """
    Plot the QSO and PSFs SB compare.
    ------
    norm_pix : normalized position. If norm_pix = 'no', means no normalizaion.
    """
    psfs_not_none = [x for x in psfs if x is not None][0]
    if gridspace == None and astrodrz==False:
        radius = 6
    elif gridspace == None and astrodrz==True:
        radius = 8
    elif gridspace == 'log':
        radius = len(psfs_not_none)/2
    frame_mid = len(QSO)/2
    if include_QSO == True:
        if plt_QSO ==True:
            print "Plot for QSO:"
        center_QSO = np.reshape(np.asarray(np.where(QSO== QSO[frame_mid-20:frame_mid+20,frame_mid-20:frame_mid+20].max())),(2))[::-1]
#        print "center_QSO:", center_QSO
        r_SB_QSO, r_grids_QSO = SB_profile(QSO, center=center_QSO, radius=radius, grids=grids, fits_plot=plt_QSO, gridspace=gridspace, if_annuli=if_annuli)
        if isinstance(norm_pix,int) or isinstance(norm_pix,float):
            count = r_grids_QSO <= norm_pix
            idx = count.sum()-1
#            print "idx:",idx
            r_SB_QSO /= r_SB_QSO[idx]      #normalize the curves
    psfs_NO = len(psfs)
    center = np.reshape(np.asarray(np.where(psfs_not_none== psfs_not_none[frame_mid-20:frame_mid+20,frame_mid-20:frame_mid+20].max())),(2))[::-1]
#    print "center_PSF:", center
    if plt_which_PSF != None:
        for i in range(len(plt_which_PSF)):
            j = plt_which_PSF[i]
            msk_counts, mask_lists = text_in_string_list("PSF{0}_".format(j), mask_list)
            print "Plot for fits: PSF{0}.fits".format(j)
            if msk_counts == 0:
                SB_profile(psfs[j], center, radius=radius, grids=grids, fits_plot=True, gridspace=gridspace)
            elif msk_counts >0:
                print mask_lists
                SB_profile(psfs[j], center, radius=radius, grids=grids, fits_plot=True, gridspace=gridspace,
                                       mask_plot = False, mask_list=mask_lists)
    minorLocator = AutoMinorLocator()
    fig, ax = plt.subplots(figsize=(10,7))
    for i in range(psfs_NO):
        if psfs[i] is not None:
            msk_counts, mask_lists = text_in_string_list("PSF{0}_".format(i), mask_list)
            if i ==0 and include_QSO == True:
                    plt.plot(r_grids_QSO, r_SB_QSO, '-', color = 'red', label="QSO", linewidth=5)
                    plt.legend()
            if msk_counts == 0:
                r_SB, r_grids = SB_profile(psfs[i], center, radius=radius, grids=grids, gridspace=gridspace, if_annuli=if_annuli)
            elif msk_counts >0:
                r_SB, r_grids = SB_profile(psfs[i], center, radius=radius, grids=grids, gridspace=gridspace, if_annuli=if_annuli,
                                           mask_list=mask_lists)
            if isinstance(norm_pix,int) or isinstance(norm_pix,float):
                count = r_grids <= norm_pix
                idx = count.sum() -1
#                print "idx:",idx
                r_SB /= r_SB[idx]      #normalize the curves
            if gridspace == None:
                plt.text(0.25,r_SB[0],i)
            elif gridspace == 'log':
                plt.text(0.40,r_SB[0],i)
#            print r_grids[idx]
            if i not in not_plt:
                plt.plot(r_grids, r_SB, 'x-', label="PSF{0}".format(i))
                plt.legend()
    ax.xaxis.set_minor_locator(minorLocator)
    plt.tick_params(which='both', width=2)
    plt.tick_params(which='major', length=7)
    plt.tick_params(which='minor', length=4, color='r')
    plt.grid()
    if isinstance(norm_pix,int):
        ax.set_ylabel("Scaled Surface Brightness")
    else:
        ax.set_ylabel("Surface Brightness")
    ax.set_xlabel("Pixels")
    if gridspace == 'log':
        ax.set_xscale('log')
        plt.xlim(0.4, ) 
    plt.grid(which="minor")
#    plt.close() 
    return fig

def profiles_compare(prf_list, scal_list, prf_name_list = None,gridspace = None ,
                     grids = 20,  norm_pix = 3, if_annuli=False, astrodrz=False):
    '''
    Compare the profile between different profile (prf?). One can set the scal to uniformize the resolution.
    Note that the SB center from the center of the image; Not allow mask yet.
    Parameter
    --------
        prf_list: profile list of the arraies
        scal_list: a list for the scaled value for the resultion, default as zero, i.e. same resolution.
    Return
    --------
        A plot of SB comparison.
    '''
    if gridspace == None and astrodrz==False:
        radius = 6
    elif gridspace == None and astrodrz==True:
        radius = 8
    elif gridspace == 'log':
        radius = len(prf_list[0])/2
    
    from matplotlib.ticker import AutoMinorLocator
    minorLocator = AutoMinorLocator()
    fig, ax = plt.subplots(figsize=(10,7))
    prf_NO = len(prf_list)
    for i in range(prf_NO):
        center = np.reshape(np.asarray(np.where(prf_list[i]== prf_list[i].max())),(2))[::-1]
        scale = scal_list[i]
        r_SB, r_grids = SB_profile(prf_list[i], center, radius=radius*scale,
                                   grids=grids, gridspace=gridspace,if_annuli=if_annuli)
        
        if isinstance(norm_pix,int) or isinstance(norm_pix,float):
            count = r_grids <= norm_pix
            idx = count.sum() -1
#            print "idx:",idx
            r_SB /= r_SB[idx]      #normalize the curves
        r_grids /= scale
        if prf_name_list == None:
            plt.plot(r_grids, r_SB, 'x-', label="prf_list{0}".format(i))
        elif len(prf_name_list)==len(prf_list):
            plt.plot(r_grids, r_SB, 'x-', label=prf_name_list[i])
        else:
            raise ValueError("The profile name is not in right length")
        plt.legend()
        
    ax.xaxis.set_minor_locator(minorLocator)
    plt.tick_params(which='both', width=2)
    plt.tick_params(which='major', length=7)
    plt.tick_params(which='minor', length=4, color='r')
    plt.grid()
    ax.set_ylabel("Scaled Surface Brightness")
    ax.set_xlabel("Pixels")
    if gridspace == 'log':
        ax.set_xscale('log')
        plt.xlim(0.5, ) 
    plt.grid(which="minor")
#    plt.close() 
    return fig


def string_find_between(s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

def cr_mask_img(image, mask_list, mask_reg_cut = 0., cut_org=False):
    '''
    Creat a mask image given a mask_list
        if set cut_org=True, the input image don't need to be cutout (e.g. QSO[cut,cut...]).
    '''
    if cut_org ==False:
        mask = np.ones(image.shape)
        for i in range(len(mask_list)):
            mask *= cr_mask(image=image, filename=mask_list[i],mask_reg_cut=mask_reg_cut)
    elif cut_org ==True:
        cut = mask_reg_cut
        mask_org = np.ones(image.shape)
        for i in range(len(mask_list)):
            mask_org *= cr_mask(image=image, filename=mask_list[i])
        mask = np.ones(image[cut:-cut,cut:-cut].shape)
        mask = mask_org[cut:-cut,cut:-cut]
    return mask

def cr_mask(image, filename='test_circle.reg', mask_reg_cut = 0.):
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
    ####!!!Need to check the center, shape of the boxtype, the size?
    with open(filename, 'r') as input_file:
        reg_string=input_file.read().replace('\n', '')
    if "physicalcircle" in reg_string:
        abc=string_find_between(reg_string, "(", ")")
        reg_info=np.fromstring(abc, dtype=float, sep=',')
        center, radius = reg_info[:2]-1 - mask_reg_cut, reg_info[2]
        region = pix_region(center, radius)
        box = 1-region.to_mask(mode='center').data
    elif "physicalbox" in reg_string:
        abc=string_find_between(reg_string, "(", ")")
        reg_info=np.fromstring(abc, dtype=float, sep=',')
        center = reg_info[:2] - 1 - mask_reg_cut
        x_r, y_r = reg_info[2:4]  # x_r is the length of the x, y_r is the length of the y
        box = np.zeros([np.int(x_r)+1, np.int(y_r)+1]).T
    else:
        print reg_string
        raise ValueError("The input reg is un-defined yet")
    frame_size = image.shape
    box_size = box.shape
    x_edge = np.int(center[1]-box_size[0]/2) #The position of the center is x-y switched.
    y_edge = np.int(center[0]-box_size[1]/2)
    mask = np.ones(frame_size)
    mask_box_part = mask[x_edge:x_edge+box_size[0],y_edge: y_edge + box_size[1]]
    mask_box_part *= box
    return mask

def total_compare(label_list, flux_list,
                  facility = 'F140w' , plot_type= 4, target_ID = 'target_ID',
                  add_background=0.0, data_mask_list = None, data_cut = 0.,plot_compare=False,
                  pix_sz = 'swarp', msk_image=None):
    if facility == 'F140w':
        zp = 26.4524
    elif facility == 'F125w':
        zp = 26.2303
    elif facility == 'acs':
        zp = 25.94333
    
    if pix_sz == 'swarp':
        delatPixel = 0.127985
    elif pix_sz == 'drz06':
        delatPixel = 0.0642
    elif pix_sz == 'acs':
        delatPixel = 0.03
        
    norm = ImageNormalize(stretch=SqrtStretch())
    f = plt.figure(0, figsize=(16.75,4))
    ax1 = plt.subplot2grid((6,4), (0,0), rowspan=6)
    ax2 = plt.subplot2grid((6,4), (0,1), rowspan=6)
    ax3 = plt.subplot2grid((6,4), (0,2), rowspan=6)
    ax4 = plt.subplot2grid((6,4), (0,3), rowspan=5)
    ax5 = plt.subplot2grid((6,4), (5,3), sharex=ax4)
    c_ax1 = ax1.imshow(flux_list[0] + add_background,origin='lower',cmap='Greys', norm=norm, vmax = flux_list[0].max()/5 )
    clim=c_ax1.properties()['clim']
    frame_size = len(flux_list[0])
    ax1.set_ylabel(target_ID, fontsize=15)
    ax1.text(frame_size*0.05, frame_size*0.9, label_list[0],
         fontsize=20)
    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)
    scale_bar(ax1, frame_size, dist=1/delatPixel, text='1"')
    coordinate_arrows(ax1, frame_size, arrow_size=0.03)
    
    ax2.imshow(flux_list[1] + flux_list[2] + add_background,origin='lower',cmap='Greys', norm=norm, clim=clim)
    pos2_o = ax2.get_position() # get the original position
    pos2 = [pos2_o.x0 -0.03, pos2_o.y0, pos2_o.width, pos2_o.height]
    ax2.set_position(pos2) # set a new position
    ax2.text(frame_size*0.05, frame_size*0.9, label_list[-2],
         fontsize=20)
    scale_bar(ax2, frame_size, dist=1/delatPixel, text='1"')
    coordinate_arrows(ax2, frame_size, arrow_size=0.03)
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    
    ax3.imshow(flux_list[0]-(flux_list[1]+flux_list[2]),origin='lower',cmap='Greys', norm=norm, clim=clim)
    ax3.text(frame_size*0.05, frame_size*0.9, label_list[-1],
         fontsize=20)
    pos3_o = ax3.get_position() # get the original position
    pos3 = [pos3_o.x0 -0.06, pos3_o.y0, pos3_o.width, pos3_o.height]
    scale_bar(ax3, frame_size, dist=1/delatPixel, text='1"')
    coordinate_arrows(ax3, frame_size, arrow_size=0.03)
    ax3.set_position(pos3) # set a new position
    ax3.get_xaxis().set_visible(False)
    ax3.get_yaxis().set_visible(False)
#    f.colorbar(ax1_c, ax=ax1.ravel().tolist())
#    plt.colorbar(c_ax1)
    make_ticklabels_invisible(plt.gcf())
    for i in range(len(flux_list)-1):
        if i == 0:
            model_flux = flux_list[i+1] +0  # Don't share a same space
        else:
            model_flux += flux_list[2]
    model_flux = flux_list[1] + flux_list[2]

    label_SB_list = [label_list[0], label_list[-2], label_list[1], label_list[2]] 
    flux_SB_list = [flux_list[0], model_flux, flux_list[1], flux_list[2]]
    for i in range(len(label_SB_list)):
        center = len(flux_SB_list[i])/2, len(flux_SB_list[i])/2
        if label_SB_list[i] == 'data':
            print "data_mask_lists:\t", data_mask_list
            r_SB, r_grids = SB_profile(flux_SB_list[i], center, gridspace = 'log',
                                       radius= 20, grids = 40, mask_list=data_mask_list,
                                       mask_cut = data_cut, msk_image=msk_image)
        else:
            r_SB, r_grids = SB_profile(flux_SB_list[i], center, gridspace = 'log', radius= 20, mask_list=None)
        r_mag = - 2.5 * np.log10(r_SB) + zp 
        if label_SB_list[i] == 'data':
            ax4.plot(r_grids * delatPixel, r_mag, 'o', color = 'whitesmoke',markeredgecolor="black", label=label_SB_list[i])
        else:
            ax4.plot(r_grids * delatPixel, r_mag, '-', label=label_SB_list[i])
    ax4.set_xscale('log')
    ax4.invert_yaxis()
    ax4.set_ylabel('$\mu$(mag, pixel$^{-2}$)', fontsize=12)
    plt.gca().invert_yaxis()
    ax4.legend()
    pos4_o = ax4.get_position() # get the original position
    pos4 = [pos4_o.x0 -0.04, pos4_o.y0 + 0.08, pos4_o.width, pos4_o.height*0.9]
    ax4.set_position(pos4) # set a new position

    x = np.linspace(1.e-4, 100, 2)
    y = x * 0
    r_mag_0 = 2.5 * np.log10(SB_profile(flux_SB_list[0], center, gridspace = 'log',
                                        radius= 20, mask_list=data_mask_list, mask_cut = data_cut,
                                        msk_image=msk_image)[0])
    r_mag_1 = 2.5 * np.log10(SB_profile(flux_SB_list[1], center, gridspace = 'log', radius= 20)[0])
    ax5.plot(r_grids*delatPixel, r_mag_0-r_mag_1, 'ro')   
    ax5.set_ylabel('$\Delta\mu$', fontsize=15)
    ax5.set_xlabel('arcsec', fontsize=15)
    ax5.set_xscale('log')
    ax5.set_xticks([0.1, 0.2, 0.5, 1, 3])
    ax5.set_yticks([-1.0, -0.5, 0., 0.5, 1.0])
    import matplotlib
    ax5.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax5.plot(x, y, 'k--')  
    plt.xlim([(r_grids*delatPixel).min()* 0.7,(r_grids*delatPixel).max() + 1])
    plt.ylim([-1,1])
    pos5_o = ax5.get_position() # get the original position
    pos5 = [pos5_o.x0 -0.04, pos5_o.y0, pos5_o.width, pos5_o.height*2]
    ax5.set_position(pos5) # set a new position
    if plot_compare == True:
        plt.show()
    return f
#
#def data_model_compare(label_list, flux_list, fig = None):
#    center = len(flux_list[0])/2, len(flux_list[0])/2
#    for i in range(len(label_list)):
#        r_SB, r_grids = SB_profile(flux_list, center)
#        fig.plot(r_grids, r_SB, 'x-', label="label_list{0}".format(i))
#        fig.legend()
    
def make_ticklabels_invisible(fig):
    for i, ax in enumerate(fig.axes):
#        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        if i !=3 and i !=4:
            for tl in ax.get_xticklabels() + ax.get_yticklabels():
                tl.set_visible(False)

def coordinate_arrows(ax, d, color='black', arrow_size=0.02):
    d0 = d / 12.
    p0 = d / 12.
    pt = d / 7.
    deltaPix = 1
    ra0, dec0 = (d - d0) / deltaPix, d0 / deltaPix
    xx_, yy_ = ra0, dec0
    xx_ra, yy_ra = (ra0 - p0, dec0)
    xx_dec, yy_dec = (ra0, dec0 + p0)
    xx_ra_t, yy_ra_t = (ra0 - pt, dec0)
    xx_dec_t, yy_dec_t = (ra0, dec0 + pt)

    ax.arrow(xx_ * deltaPix, yy_ * deltaPix, (xx_ra - xx_) * deltaPix, (yy_ra - yy_) * deltaPix,
             head_width=arrow_size * d, head_length=arrow_size * d, fc=color, ec=color, linewidth=1.2)
    ax.text(xx_ra_t * deltaPix, yy_ra_t * deltaPix, "E", color=color, fontsize=12, ha='center')
    ax.arrow(xx_ * deltaPix, yy_ * deltaPix, (xx_dec - xx_) * deltaPix, (yy_dec - yy_) * deltaPix,
             head_width=arrow_size * d, head_length=arrow_size * d, fc
             =color, ec=color, linewidth=1.2)
    ax.text(xx_dec_t * deltaPix, yy_dec_t * deltaPix, "N", color=color, fontsize=12, ha='center')
    
def scale_bar(ax, d, dist=1/0.13, text='1"', color='black', flipped=False):
    if flipped:
        p0 = d - d / 15. - dist
        p1 = d / 15.
        ax.plot([p0, p0 + dist], [p1, p1], linewidth=2, color=color)
        ax.text(p0 + dist / 2., p1 + 0.02 * d, text, fontsize=15, color=color, ha='center')
    else:
        p0 = d / 15.
        ax.plot([p0, p0 + dist], [p0, p0], linewidth=2, color=color)
        ax.text(p0 + dist / 2., p0 + 0.02 * d, text, fontsize=15, color=color, ha='center')