#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 16:45:01 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import copy

import matplotlib as matt
import matplotlib.lines as mlines
from matplotlib import colors
matt.rcParams['font.family'] = 'STIXGeneral'

filt_info = {'CDFS-1': 'F140w', 'CDFS-229': 'F125w', 'CDFS-724': 'F125w',\
'CID1174': 'F140w', 'CID206': 'F140w', 'CID216': 'F140w', 'CID3242': 'F140w',\
'CID3570': 'F125w', 'CID452': 'F125w', 'CID454': 'F140w', 'CID50': 'F125w',\
'CID607': 'F125w', 'CID70': 'F140w', 'LID1273': 'F140w', 'LID360': 'F140w',\
'XID2138': 'F140w', 'XID2202': 'F140w', 'XID2396': 'F140w', 'ECDFS-358': 'F140w',\
'SXDS-X1136': 'F125w', 'SXDS-X50': 'F125w', 'SXDS-X735': 'F140w',\
'CID543': 'F125w', 'LID1538': 'F140w', 'CID237': 'F140w', 'SXDS-X717': 'F125w',\
'SXDS-X763': 'F125w', 'SXDS-X969': 'F140w', 'CDFS-321': 'F140w'}


PSFs_dict = {}
QSOs_dict = {}
for key in filt_info.keys():
    ID = key
    PSFs, QSOs=pickle.load(open('{0}/analysis/{0}_PSFs_QSO'.format(ID),'rb'))
    PSFs_dict.update({'{0}'.format(ID):PSFs})
    QSOs_dict.update({'{0}'.format(ID):QSOs})

PSF_list = []
PSF_id = []
QSO_list = []
filter_list = []

IDs = PSFs_dict.keys()

for ID in IDs:
    psfs_dict = PSFs_dict[ID]
    psfs = [psfs_dict[i] for i in range(len(psfs_dict))]
    PSF_list += psfs
    name_id = [ID+"_"+str(i) for i in range(len(psfs_dict))]
    PSF_id = PSF_id + name_id
    filt = [filt_info[ID]]
    filter_list += filt * len(PSFs_dict[ID])
    qso = QSOs_dict[ID]
    QSO_list.append(qso)
    if len(PSF_list) != len(PSF_id):
        raise ValueError("The PSF_list is not consistent with PSF_id")
QSO_id = IDs
#==============================================================================
#Plot the flux distribution 
#==============================================================================
fluxs = [] 
locs = []
id_star_s = []
for i in range(len(QSO_id)):
    PSFs = [PSF_list[j] for j in range(len(PSF_list)) if QSO_id[i] in PSF_id[j]]
#    print len(PSFs)
    flux = [np.sum(PSFs[j][0] * PSFs[j][3]) for j in range(len(PSFs))]
    fluxs += flux
    loc = [PSFs[j][2] for j in range(len(PSFs))]
    locs += loc
    id_star = [PSFs[j][1] for j in range(len(PSFs))]
    id_star_s += id_star
plt.figure(figsize=(10,6))
common_params = dict(bins=30)
#                         label=QSO_id)
plt.hist((fluxs), **common_params)
plt.legend()
plt.xticks(np.arange(0, 4500, step=500))
plt.yticks(np.arange(0,65,step=5))
plt.xlabel('Total Flux (counts)',fontsize=15)
plt.ylabel('Number of PSFs', fontsize=15)
plt.tick_params(labelsize=15)
plt.show()

#==============================================================================
#Location 
#==============================================================================
#temp_frame = pyfits.getdata('CDFS-1/astrodrz/final_drz.fits')
frame_y, frame_x = (2183, 2467) #temp_frame.shape
x_len = 15.
ratio = x_len/frame_x
y_len = frame_y * ratio
#fig, ax= plt.subplots(1)
plt.figure(figsize=(x_len,y_len))
color_label = copy.deepcopy(QSO_id)
for i in range(len(QSO_id)):
    loc = QSO_list[i][1]
    plt.plot(loc[0]*ratio, loc[1]*ratio, marker='X', label = color_label[i], ms = 20,  linestyle = 'None')
    color_label[i] = "_nolegend_"
    PSFs = [PSF_list[j] for j in range(len(PSF_list)) if QSO_id[i] in PSF_id[j]]
    for j in range(len(PSFs)):
        loc = PSFs[j][2]
        if PSFs[j][1] == 1:
            plt.plot(loc[0]*ratio, loc[1]*ratio, marker='*', ms = 10, linestyle = 'None')
        elif PSFs[j][1] == 0:
            plt.plot(loc[0]*ratio, loc[1]*ratio, marker='o', ms = 7, linestyle = 'None')

dith_fx, dith_fy = (2128,1916)
dith_fx *= ratio
dith_fy *= ratio
box_wx, box_wy = dith_fx, dith_fy
dx = (x_len-dith_fx)/5
dy = (y_len-dith_fy)/5
for i in range(6):
    rectangle = plt.Rectangle((dx*i, dy*i), box_wx, box_wy, fill=None, alpha=1)
    plt.gca().add_patch(rectangle)

plt.tick_params(labelsize=20)
plt.legend(bbox_to_anchor=(0.97, 1.14),prop={'size': 12},ncol=7)
plt.xticks([]),plt.yticks([])
plt.xlim(0, x_len)
plt.ylim(0, y_len)
plt.show()

##==============================================================================
## Plot mask
##==============================================================================
#from matplotlib.colors import LogNorm
#ncols = 3
#nrows = 3
#count = 0
#import math
#ceils = math.ceil(len(PSF_list)/9.)
#for k in range(int(ceils)):
#    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(9,9))
#    for i in range(ncols):
#        for j in range(nrows):
#            if count < len(PSF_list) :
#                axs[i,j].imshow(PSF_list[count][0]*PSF_list[count][3], origin='low', norm=LogNorm())
#                t = axs[i,j].text(5,110,PSF_id[count], {'color': 'b', 'fontsize': 20})
#                t.set_bbox(dict(facecolor='yellow', alpha=0.5, edgecolor='red'))
#                if PSF_list[count][1] ==1:
#                    t1 = axs[i,j].text(100,5,'star', {'color': 'b', 'fontsize': 10})
#                    t1.set_bbox(dict(facecolor='yellow', alpha=0.5, edgecolor='red'))
#            axs[i,j].set_xticks([])
#            axs[i,j].set_yticks([])
#            count += 1
#    plt.tight_layout()
#    plt.show()
    
#==============================================================================
# Gaussian Fit
#==============================================================================
from scipy.optimize import curve_fit
test_data = PSF_list[0][0]
center = len(test_data)/2
#def func(x, a, x0, sigma):
#    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def measure_FWHM(image,count=0 ,line_range= (50,70)):
    seed = range(line_range[0],line_range[1])
    frm = len(image)
    q_frm = frm/4
    center = np.where(image == image[q_frm:-q_frm,q_frm:-q_frm].max())[0][0]
    x_n = np.asarray([image[i][center] for i in seed]) # The x value, vertcial 
    y_n = np.asarray([image[center][i] for i in seed]) # The y value, horizontal 
    #popt, pcov = curve_fit(func, x, yn)
    from astropy.modeling import models, fitting
    g_init = models.Gaussian1D(amplitude=1., mean=60, stddev=1.)
    fit_g = fitting.LevMarLSQFitter()
    g_x = fit_g(g_init, seed, x_n)
    g_y = fit_g(g_init, seed, y_n)
    FWHM_ver = g_x.stddev.value * 2.355  # The FWHM = 2*np.sqrt(2*np.log(2)) * stdd = 2.355*stdd
    FWHM_hor = g_y.stddev.value * 2.355
    if (FWHM_hor-FWHM_ver)/FWHM_ver > 0.20:
        print "Warning, the {0} have inconsistent FWHM".format(count), FWHM_ver, FWHM_hor
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    seed = np.linspace(seed.min(),seed.max(),50)
#    ax.plot(seed, g_x(seed), c='r', label='Gaussian')
#    ax.legend()
#    ax.scatter(x, sample)
#    plt.show()
    return (FWHM_ver+FWHM_hor)/2., FWHM_ver, FWHM_hor
FWHM = []
for i in range(len(PSF_list)):
    FWHM_i = measure_FWHM(PSF_list[i][0], count = i)[0]
    FWHM.append(FWHM_i)
FWHM = np.asarray(FWHM)

fig, ax = plt.subplots(figsize=(10,6))
color_dict = {"F140w": 'b', "F125w": 'g'}
label = {"F140w1":'F140w star', "F140w0":'F140w selected', "F125w1":'F125w star', "F125w0":'F125w selected'}
marker =  {1:'*', 0:'o'}
for i in range(len(FWHM)):
    filt = filter_list[i]
    label_key = filt + str(PSF_list[i][1])
    ax.scatter(FWHM[i], fluxs[i], color =color_dict[filt], label = label[label_key], marker=marker[PSF_list[i][1]])
    label[label_key] = "_nolegend_"
import matplotlib
ax.set_yscale('log')
ax.set_yticks([50,100, 200, 300, 500, 1000, 2000, 5000])
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.xlabel('FWHM (pixel)',fontsize=15)
plt.ylabel('Total flux', fontsize=15)
plt.legend(prop={'size': 12})
plt.tick_params(labelsize=15)
plt.ylim(30,6000)
plt.show()

flux_dict, FWHM_dict, locs_dict, filter_dict, id_stars_dict = {}, {}, {}, {}, {}
for i in range(len(PSF_list)):
    if filter_list[i] == 'F140w':
        print PSF_id[i], round(fluxs[i],3), round(FWHM[i],3), locs[i], filter_list[i], id_star_s[i]
        flux_dict.update({PSF_id[i]:fluxs[i]})
        FWHM_dict.update({PSF_id[i]:FWHM[i]})
        locs_dict.update({PSF_id[i]:locs[i]})
        filter_dict.update({PSF_id[i]:filter_list[i]})
        id_stars_dict.update({PSF_id[i]:id_star_s[i]})
print '\n'

for i in range(len(PSF_list)):
    if filter_list[i] == 'F125w':
        print PSF_id[i], round(fluxs[i],3), round(FWHM[i],3), locs[i], filter_list[i], id_star_s[i]
        flux_dict.update({PSF_id[i]:fluxs[i]})
        FWHM_dict.update({PSF_id[i]:FWHM[i]})
        locs_dict.update({PSF_id[i]:locs[i]})
        filter_dict.update({PSF_id[i]:filter_list[i]})
        id_stars_dict.update({PSF_id[i]:id_star_s[i]})

#import pickle
#filename='PSFs_lib_dict'.format(ID)
#datafile = open(filename, 'wb')
#PSFs_lib = [flux_dict, FWHM_dict, locs_dict, filter_dict, id_stars_dict]
#pickle.dump(PSFs_lib, open(filename, 'wb'))
#datafile.close()