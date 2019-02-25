#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 18:25:21 2019

@author: Dartoon

Read the overall fitting result
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import re
import matplotlib
from math import log10
from adjustText import adjust_text   # avoid the overlapping while ploting
import sys
sys.path.insert(0,'../py_tools')
from filter_info import filt_info, redshift_info
import pickle

def load_result(ID, flt, count_rank=8, outputBHmag = 0):
    '''
    Read the fitting result by given the ID, 
    
    Parameter
    --------
        ID: The blash of blash
        b: The blash of blash
        sort_PSF_by: 'fit_result_each' or 'fit_ps_each'
    Return
    --------
        A sth sth
    '''
    ratio_results, Re_results, n_results, total_flux_results, host_amp_results, mag_result, mag_lh_result = [], [], [], [], [], [], []
    chisq_list, inf_list, best_PSF_id = [],[], []
    flux_dict, FWHM_dict, locs_dict, filter_dict, id_stars_dict=pickle.load(open('../analysis_Astrodrz/PSFs_lib_dict','rb'))
    if flt == 'WFC3':
        filt = filt_info[ID]
        f = open("../analysis_Astrodrz/{0}/deep_analysis/fit_result_each/each_PSF_fit_qso.txt".format(ID),"r")
    elif flt == 'ACS':
        filt = 'F814w'
        f = open("../analysis_ACS/{0}/fit_result_each_fix/each_PSF_fit_qso.txt".format(ID),"r")
    string = f.read()
#    PSF_id = re.findall(r"by PSF(.*?):",string)
    S_n_list = re.findall(r"n_sersic':(.*?),",string)
    Re = re.findall(r"R_sersic':(.*?),",string)
    host_flux_ratio = re.findall(r"host_flux_ratio_percent':(.*?)}",string)
    Chisq = re.findall(r"redu_Chisq':(.*?),",string)
    re.findall(r"redu_Chisq':(.*?),",string)
    QSO_amp = re.findall(r"QSO_amp':(.*?),",string)
    host_amp = re.findall(r"host_amp':(.*?),",string)
    
    S_n_list = [float(value) for value in S_n_list]
    Re = [float(i) for i in Re]
    host_flux_ratio = [float(i) for i in host_flux_ratio]
    Chisq = [float(i) for i in Chisq]
    QSO_amp = [float(i) for i in QSO_amp]
    host_amp =  [float(i) for i in host_amp]
    bh_amp = [host_amp[i]*(1/host_flux_ratio[i]*100. -1) for i in range(len(host_amp))]
    total_flux = np.asarray(QSO_amp) + np.asarray(host_amp)
    sort_Chisq = np.argsort(np.asarray(Chisq))
    count_n = count_rank
    Chisq_best = Chisq[sort_Chisq[0]]
    Chisq_last= Chisq[sort_Chisq[count_n-1]]
    inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
    weight = np.zeros(len(Chisq))
    for i in sort_Chisq[:count_n]:
        weight[i] = np.exp(-1/2. * (Chisq[i]-Chisq_best)/(Chisq_best* inf_alp))
    # =============================================================================
    # Weighting result
    # =============================================================================
    weighted_host_ratio = np.sum(np.array(host_flux_ratio)*weight) / np.sum(weight)
    rms_host_ratio = np.sqrt(np.sum((np.array(host_flux_ratio)-weighted_host_ratio)**2*weight) / np.sum(weight))
    weighted_Re = np.sum(np.array(Re)*weight) / np.sum(weight)
    rms_Re = np.sqrt(np.sum((np.array(Re)-weighted_Re)**2*weight) / np.sum(weight))
    weighted_index = np.sum(np.array(S_n_list)*weight) / np.sum(weight)
    rms_index = np.sqrt(np.sum((np.array(S_n_list)-weighted_index)**2*weight) / np.sum(weight))
    weighted_total_flux = np.sum(total_flux*weight) / np.sum(weight)
    rms_total_flux = np.sqrt(np.sum((total_flux-weighted_total_flux)**2*weight) / np.sum(weight))
    weighted_host_flux = np.sum(np.asarray(host_amp)*weight) / np.sum(weight)
    rms_host_flux = np.sqrt(np.sum((np.asarray(host_amp)-weighted_host_flux)**2*weight) / np.sum(weight))
    weighted_bh_flux = np.sum(np.asarray(bh_amp)*weight) / np.sum(weight)
    rms_bh_flux = np.sqrt(np.sum((np.asarray(bh_amp)-weighted_bh_flux)**2*weight) / np.sum(weight))
    if filt == "F140w":
        zp = 26.4524
    elif filt == "F125w":
        zp = 26.2303
    elif filt == "F814w":
        zp = 25.94333
        
    amp_min, amp, amp_max = weighted_host_flux - rms_host_flux, weighted_host_flux, weighted_host_flux + rms_host_flux
    mag_min =-2.5*np.log10(amp_min) + zp
    mag =-2.5*np.log10(amp) + zp
    mag_max =-2.5*np.log10(amp_max) + zp
    mag_result = mag 
    mag_lh_result = [mag_min-mag, mag-mag_max]
    
    bh_min, bh_amp, bh_max = weighted_bh_flux - rms_bh_flux, weighted_bh_flux, weighted_bh_flux + rms_bh_flux
    bh_mag_min =-2.5*np.log10(bh_min) + zp
    bh_mag =-2.5*np.log10(bh_amp) + zp
    bh_mag_max =-2.5*np.log10(bh_max) + zp
    bh_mag_result =bh_mag
    bh_mag_lh_result = [bh_mag_min-bh_mag, bh_mag-bh_mag_max]

    best_PSF_id.append(sort_Chisq[:8])
    ratio_results = [weighted_host_ratio, rms_host_ratio]
    Re_results = [weighted_Re, rms_Re]
    n_results = [weighted_index,rms_index]
    host_amp_results = [weighted_host_flux, rms_host_flux]
    total_flux_results = [weighted_total_flux, rms_total_flux]
    chisq_list.append(Chisq_best)
    inf_list.append(inf_alp)  
    if outputBHmag != 1 :
        return ratio_results, Re_results, n_results, host_amp_results, total_flux_results, mag_result, mag_lh_result, Chisq_best
    elif outputBHmag == 1 :
        return ratio_results, Re_results, n_results, host_amp_results, total_flux_results, mag_result, mag_lh_result, Chisq_best, bh_mag_result, bh_mag_lh_result
    
    
WFC3_list = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174', 'CID216', 'LID360',
'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597', 'CID1281'\
]
WFC3_list.sort()

tab_list = ['CID1174', 'CID1281', 'CID206', 'CID216',
            'CID237', 'CID255', 'CID3242', 'CID3570', 'CID452', 'CID454',
            'CID50', 'CID543', 'CID597', 'CID607', 'CID70', 'LID1273', 'LID1538',
            'LID360', 'XID2138', 'XID2202', 'XID2396', 'CDFS-1', 'CDFS-229', 'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']
tab_list.sort()

#filt_list, ACS_list = [],[]
#z_list = []
#IR_result, UV_result = [], []
#for i in range(len(tab_list)):
#    if tab_list[i] in WFC3_list:   #Input the WFC3 filter
#        filt_list.append(filt_info[tab_list[i]])
#        z_list.append(redshift_info[tab_list[i]])
#        IR_result.append(load_result(tab_list[i],'WFC3', outputBHmag=1))
#    else:
#        filt_list.append("xxx")
#        z_list.append('xxx')
#        IR_result.append('xxx')
        
#for i in range(len(tab_list)):
#    if IR_result[i]=='xxx':
#        IR_result[i] = [np.asarray(IR_result[0][j])*0-99 for j in range(len(IR_result[0]))]
#    print(tab_list[i], round(IR_result[i][-3],3),
#    "{0}%pm{1}%".format(round(IR_result[i][0][0],1), round(IR_result[i][0][1],1)),  #host flux ratio
#    "{0}pm{1}".format(round(IR_result[i][1][0],3), round(IR_result[i][1][1],3)),  #host Re
#    "{0}pm{1}".format(round(IR_result[i][2][0],3), round(IR_result[i][2][1],3)),  #host n
#    "{0}+{1}-{2}".format(round(IR_result[i][5],3), round(IR_result[i][6][0],3),round(IR_result[i][6][1],3)))
#
# =============================================================================
# Load the result:
# =============================================================================
n_results, ratio_results, Re_results = [], [], []
ID = []
for i in range(len(tab_list)):
    if tab_list[i] in WFC3_list:
        ID.append(tab_list[i])
        result = load_result(tab_list[i],'WFC3', outputBHmag=1)
        ratio_results.append([result[0][0], result[0][1]])
        Re_results.append([result[1][0], result[1][1]])
        n_results.append([result[2][0], result[2][1]])

import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
cmap = matplotlib.cm.get_cmap('viridis')
normalize = matplotlib.colors.Normalize(vmin=0.3, vmax=7.0)
colors = [cmap(normalize(value[0])) for value in n_results]
#==============================================================================
# Sersic_n, R_eff relation, host_flux_ratio corner relation
#==============================================================================
num_boxs = 2
gridshape = (num_boxs, num_boxs)
print "Our multivariate grid will therefore be of shape", gridshape
fig = plt.figure(figsize=(20, 15))
axes = [[False for i in range(num_boxs)] for j in range(num_boxs)]
axis_lab = [ "Effective Radius(arcsec)", "S\'ersic $n$", "host flux ratio percent (%)"]
value = [[Re_results[i], n_results[i], ratio_results[i]] for i in range(len(ratio_results))]
n=1
for j in range(num_boxs):
    for i in range(num_boxs):
        if i <= j :
            y_j = j+1
            ax = fig.add_subplot(num_boxs, num_boxs, n)
            plt.setp(ax.spines.values(),linewidth=2)
            ax.tick_params(labelsize=12)
            ax.get_xaxis().set_tick_params(direction='in', width=1.5, length=6)
            ax.get_yaxis().set_tick_params(direction='in', width=1.5, length=6)
            ax.yaxis.set_ticks_position('both')
            ax.xaxis.set_ticks_position('both')
            texts = []
            for k in range(len(value)):
                filt = filt_info[ID[i]]
                if filt == "F140w":
                    ma = 'o'
                elif filt == "F125w":
                    ma = 's'
                ax.errorbar(value[k][i][0], value[k][y_j][0], xerr= value[k][i][1], yerr= value[k][y_j][1], marker = ma, color = colors[k])
                texts.append(ax.text(value[k][i][0], value[k][y_j][0], ID[k], fontsize=15))
            adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red')) 
            if i == 0:
#                plt.ylabel('j+1={0}={1}'.format(axis_lab[y_j],y_j), fontsize=15)
                plt.ylabel('{0}'.format(axis_lab[y_j]), fontsize=25)
#                plt.xlim(0,5)   #plot the limit for goodness for x axis
#            elif i == 1:
#                plt.xlim(0,20)  #plot the limit for precision for x axis
            if y_j == 2:
#                plt.xlabel('i={0}={1}'.format(axis_lab[i],i), fontsize =15)
                plt.xlabel('{0}'.format(axis_lab[i]), fontsize =25)
#                plt.ylim(-20,20)   #plot the limit for accuracy for y axis
#            elif y_j ==1:
#                plt.ylim(0,20)   #plot the limit for precision for y axis
            if i>=1:
                ax.yaxis.set_ticklabels([])
            if y_j <2:
                ax.xaxis.set_ticklabels([])
            plt.tick_params(labelsize=25)     
        fig.tight_layout(h_pad=-1.,w_pad=-0.3)
        if i== 0 and j ==0:
            cax, _ = matplotlib.colorbar.make_axes(ax)
            cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
            cbar.set_label('Sersic n',  fontsize=25)
            cax.tick_params(labelsize=30)
            pos_o = cax.get_position()
            pos = [pos_o.x0+0.09, pos_o.y0, pos_o.width, pos_o.height]
            cax.set_position(pos)
        n += 1
#plt.savefig('flux_r_n_corner.pdf')
plt.show()