#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 14:47:47 2018

@author: dxh
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import re
import matplotlib
from adjustText import adjust_text   # avoid the overlapping while ploting
import sys
sys.path.insert(0,'../py_tools')
from filter_info import filt_info, redshift_info

ID = ['CID50','CID70','XID2138','CID3242',\
'LID1273','XID2202','CID206','CID543','LID1538','XID2396','CID452',\
'LID360','CID237','CID454','CID607','CID3570']

#CID1281, CID255, CID526, CID597 is not in our targets

import pickle
ratio_results, Re_results, n_results, total_flux_results, host_amp_results = [], [], [], [], []
chisq_list, inf_list, best_PSF_id = [],[], []
flux_dict, FWHM_dict, locs_dict, filter_dict, id_stars_dict=pickle.load(open('PSFs_lib_dict','rb'))
for j in range(len(ID)):
    f = open("{0}/fit_result_each/each_PSF_fit_qso.txt".format(ID[j]),"r")
    string = f.read()
    PSF_id = re.findall(r"by PSF(.*?):",string)
    S_n_list = re.findall(r"n_sersic':(.*?),",string)
    Re = re.findall(r"R_sersic':(.*?),",string)
    host_flux_ratio = re.findall(r"host_flux_ratio_percent':(.*?)}",string)
    Chisq = re.findall(r"redu_Chisq':(.*?),",string)
    QSO_amp = re.findall(r"QSO_amp':(.*?),",string)
    host_amp = re.findall(r"host_amp':(.*?),",string)
    
    S_n_list = [float(value) for value in S_n_list]
    Re = [float(i) for i in Re]
    host_flux_ratio = [float(i) for i in host_flux_ratio]
    Chisq = [float(i) for i in Chisq]
    QSO_amp = [float(i) for i in QSO_amp]
    host_amp =  [float(i) for i in host_amp]
    
#    PSFs_dict = {}
#    for key in filt_info.keys():
#        if filt_info[key] == filt:
#            PSFs, _=pickle.load(open('{0}/{0}_PSFs_QSO'.format(key),'rb'))
#            PSFs_dict.update({'{0}'.format(key):PSFs})
    total_flux = np.asarray(QSO_amp) + np.asarray(host_amp)
    
    sort_Chisq = np.argsort(np.asarray(Chisq))
    count_n = 8
    Chisq_best = Chisq[sort_Chisq[0]]
    Chisq_last= Chisq[sort_Chisq[count_n-1]]
    inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
    weight = np.zeros(len(Chisq))
    for i in sort_Chisq[:count_n]:
        weight[i] = np.exp(-1/2. * (Chisq[i]-Chisq_best)/(Chisq_best* inf_alp))
#    for i in sort_Chisq:
#        print PSF_id[i], Chisq[i], weight[i], host_flux_ratio[i], Re[i], S_n_list[i], round(total_flux[i],3),  round(host_amp[i],3),\
#    round(flux_dict[PSF_id[i]],3), round(FWHM_dict[PSF_id[i]],3), "[{0},{1}]".format(int(round(locs_dict[PSF_id[i]][0])) , int(round(locs_dict[PSF_id[i]][1]))), round(id_stars_dict[PSF_id[i]],3)
    
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
    
    print "ID:{0}".format(ID[j])
    print "the best PSFs for {0} is ID ".format(ID[j]), sort_Chisq[:8]
    best_PSF_id.append(sort_Chisq[:8])
    print "Weighted: HOST_Ratio, Reff, Sersic_n, total_flux, host_flux: {0}+-{1} {2}+-{3} {4}+-{5} {6}+-{7} {8}+-{9}".format(\
    round(weighted_host_ratio,3), round(rms_host_ratio,3),\
    round(weighted_Re,3), round(rms_Re,3), round(weighted_index,3),round(rms_index,3),\
    round(weighted_total_flux,3), round(rms_total_flux,3),\
    round(weighted_host_flux,3), round(rms_host_flux,3))
    ratio_results.append([round(weighted_host_ratio,3), round(rms_host_ratio,3)])
    Re_results.append([round(weighted_Re,3), round(rms_Re,3)])
    n_results.append([round(weighted_index,3),round(rms_index,3)])
    host_amp_results.append([round(weighted_host_flux,3), round(rms_host_flux,3)])
    total_flux_results.append([round(weighted_total_flux,3), round(rms_total_flux,3)])
    chisq_list.append(Chisq_best)
    inf_list.append(inf_alp)

import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
cmap = matplotlib.cm.get_cmap('viridis')
normalize = matplotlib.colors.Normalize(vmin=0.3, vmax=7.0)
S_n_list = [value for value in n_results]
colors = [cmap(normalize(value[0])) for value in S_n_list]

fig, ax = plt.subplots(figsize=(15,9))
for i in range(len(ID)):
    filt = filt_info[ID[i]]
#    if filt == "F140w":
#        s = 'o'
#    elif filt == "F125w":
    s = 's'
    ax.errorbar(Re_results[i][0], ratio_results[i][0], xerr= Re_results[i][1], yerr=ratio_results[i][1],
                fmt=s, color=colors[i],ecolor='gray' )
plt.tick_params(labelsize=15)
plt.xlabel('Effective Radius',fontsize=15)
plt.ylabel('Host flux ratio', fontsize=15)
plt.tick_params(labelsize=15)
texts = []
for i in range(len(S_n_list)):
    texts.append(ax.text(float(Re_results[i][0])+0.002, ratio_results[i][0]+0.2, ID[i], fontsize=12))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))    
cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
cbar.set_label('Sersic n',  fontsize=12)
plt.show()

for i in range(len(ID)):
    print "open {0}/fit_result_each/*_PSF{1}_f*".format(ID[i],best_PSF_id[i][0])

#==============================================================================
# Plot Chisq and alpha distribution
#==============================================================================
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(chisq_list, inf_list, 'bo')
plt.xlabel('Chisq best',fontsize=15)
plt.ylabel("$alpha$", fontsize=15)
texts = []
for i in range(len(chisq_list)):
     texts.append(ax.text(float(chisq_list[i])+0.002, inf_list[i]+0.002, ID[i], fontsize=12))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))    
plt.tick_params(labelsize=15)
plt.show()

fit_mag= []
fit_mag_lh = []
from math import log10
for i in range(len(ID)):
    filt = filt_info[ID[i]]
    zp = 25.94333
    amp_min, amp, amp_max = host_amp_results[i][0] - host_amp_results[i][1],host_amp_results[i][0], host_amp_results[i][0] + host_amp_results[i][1]
    mag_min =-2.5*log10(amp_min) + zp
    mag =-2.5*log10(amp) + zp
    mag_max =-2.5*log10(amp_max) + zp
#    print "mag for", ID[i], "{0}~{1}~{2}".format(round(mag_min,3),round(mag,3),round(mag_max,3))
    fit_mag.append(mag)
    fit_mag_lh.append([mag_min-mag, mag-mag_max])

fig, ax = plt.subplots(figsize=(15,9))
for i in range(len(ID)):
    filt = filt_info[ID[i]]
#    if filt == "F140w":
#        s = 'o'
#    elif filt == "F125w":
    s = 's'
    ax.errorbar(Re_results[i][0], fit_mag[i], xerr= Re_results[i][1], yerr=np.array([fit_mag_lh[i]]).T,
                fmt=s, color=colors[i],ecolor='gray' )
plt.tick_params(labelsize=15)
plt.xlabel('Effective Radius',fontsize=15)
plt.ylabel('Host magnitude', fontsize=15)
texts = []
for i in range(len(S_n_list)):
    texts.append(ax.text(float(Re_results[i][0]), fit_mag[i], ID[i], fontsize=12))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))    
cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
cbar.set_label('Sersic n',  fontsize=12)
plt.tick_params(labelsize=15)
plt.show()   
               
        
#==============================================================================
# Sersic_n, R_eff relation, host_flux_ratio corner relation
#==============================================================================
num_boxs = 2
gridshape = (num_boxs, num_boxs)
print "Our multivariate grid will therefore be of shape", gridshape

fig = plt.figure(figsize=(20, 15))
axes = [[False for i in range(num_boxs)] for j in range(num_boxs)]
axis_lab = [ "Effective Radius", "Sersic n", "host_flux_ratio"]
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
        n += 1
#        if i==1 and j ==0:
#            ax.legend(bbox_to_anchor=(1.5, 1), loc=2, borderaxespad=0.,prop={'size': 16})
#            axes[j][i] = ax
fig.tight_layout(h_pad=-1.,w_pad=-0.6)
plt.show()

#for i in range(len(ID)):
#    print ID[i], redshift_info[ID[i]], filt_info[ID[i]], fit_mag[i], host_amp_results[i][0], ratio_results[i][0]
#    
#for i in range(len(ID)):
#    print ID[i], redshift_info[ID[i]], filt_info[ID[i]], fit_mag[i], Re_results[i][0], n_results[i][0]