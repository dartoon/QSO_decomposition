#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 20:28:30 2018

@author: Dartoon

Plot the best 8 PSFs location and compare with the QSO location.
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import re
from adjustText import adjust_text   # avoid the overlapping while ploting
import sys
sys.path.insert(0,'../py_tools')

file_list = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358',\
]

for ID in file_list:
    f = open("{0}/analysis/fit_result_each/each_PSF_fit_qso.txt".format(ID),"r")
    string = f.read()
    
    PSF_id = re.findall(r"by PSF(.*?):",string)
    Chisq = re.findall(r"redu_Chisq':(.*?),",string)
    Chisq = [float(i) for i in Chisq]
    flux_dict, FWHM_dict, locs_dict, filter_dict, id_stars_dict=pickle.load(open('PSFs_lib_dict','rb'))
    sort_Chisq = np.argsort(np.asarray(Chisq))
    count_n = 8
    
    #==============================================================================
    #Location 
    #==============================================================================
    _, QSOs=pickle.load(open('{0}/analysis/{0}_PSFs_QSO'.format(ID),'rb'))
    #temp_frame = pyfits.getdata('CDFS-1/astrodrz/final_drz.fits')
    frame_y, frame_x = (2183, 2467) #temp_frame.shape
    x_len = 15.
    ratio = x_len/frame_x
    y_len = frame_y * ratio
    #fig, ax= plt.subplots(1)
    plt.figure(figsize=(x_len,y_len))
    
    texts = []
    loc = QSOs[1]
    plt.plot(loc[0]*ratio, loc[1]*ratio, marker='X', label = ID, ms = 20,  linestyle = 'None')
    texts.append(plt.text(loc[0]*ratio, loc[1]*ratio, ID, fontsize=12))
    for i in range(count_n):
        PSF_name = PSF_id[sort_Chisq[i]]
        loc = locs_dict[PSF_name]
        if id_stars_dict[PSF_name] == 1:
            plt.plot(loc[0]*ratio, loc[1]*ratio, marker='*', ms = 10, linestyle = 'None')
            texts.append(plt.text(loc[0]*ratio, loc[1]*ratio, "NO.{0} ".format(i+1) + PSF_name, fontsize=12))
        elif id_stars_dict[PSF_name] == 0:
            plt.plot(loc[0]*ratio, loc[1]*ratio, marker='o', ms = 7, linestyle = 'None')
            texts.append(plt.text(loc[0]*ratio, loc[1]*ratio, "NO.{0} ".format(i+1) + PSF_name, fontsize=12))
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))    
    
    
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
    plt.xticks([]),plt.yticks([])
    plt.xlim(0, x_len)
    plt.ylim(0, y_len)
    plt.savefig('top8_PSFs_loc/{0}.pdf'.format(ID))
    