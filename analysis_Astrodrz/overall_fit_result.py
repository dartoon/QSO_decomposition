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

import sys
sys.path.insert(0,'../py_tools')
from filter_info import filt_info

#'CID206',
ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202', 'XID2396'
]

import pickle
ratio_results, Re_results, n_results, total_flux_results, host_amp_results = [], [], [], [], []
chisq_list, inf_list = [],[]
flux_dict, FWHM_dict, locs_dict, filter_dict, id_stars_dict=pickle.load(open('PSFs_lib_dict','rb'))
for j in range(len(ID)):
    filt = filt_info[ID[j]]
    f = open("{0}/analysis/fit_result_each/each_PSF_fit_qso.txt".format(ID[j]),"r")
    string = f.read()
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
    
    PSFs_dict = {}
    for key in filt_info.keys():
        if filt_info[key] == filt:
            PSFs, _=pickle.load(open('{0}/analysis/{0}_PSFs_QSO'.format(key),'rb'))
            PSFs_dict.update({'{0}'.format(key):PSFs})
    PSF_id = []
    filter_list = []
    for key in PSFs_dict.keys():
        psfs_dict = PSFs_dict[key]
        name_id = [key+"_"+str(i) for i in range(len(psfs_dict))]
        PSF_id = PSF_id + name_id
        filt = [filt_info[key]]
        filter_list += filt * len(PSFs_dict[key])
    psf_name_list = PSF_id
    
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
    
    print "ID:{0}, filt:{1}".format(ID[j], filt[0])
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
import matplotlib.lines as mlines
matt.rcParams['font.family'] = 'STIXGeneral'
cmap = matplotlib.cm.get_cmap('viridis')
normalize = matplotlib.colors.Normalize(vmin=0.3, vmax=7.0)
S_n_list = [value for value in n_results]
colors = [cmap(normalize(value[0])) for value in S_n_list]

fig, ax = plt.subplots(figsize=(15,9))
for i in range(len(ID)):
    ax.errorbar(Re_results[i][0], ratio_results[i][0], xerr= Re_results[i][1], yerr=ratio_results[i][1],
                fmt='o', color=colors[i],ecolor='gray' )
#plt.ylim(0, 100) 

plt.xlabel('Effective Radius',fontsize=15)
plt.ylabel('Host flux ratio', fontsize=15)

cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
cbar.set_label('Sersic n',  fontsize=12)
plt.tick_params(labelsize=15)
for i in range(len(S_n_list)):
    ax.text(float(Re_results[i][0])+0.002, ratio_results[i][0]+0.2, ID[i], fontsize=12)
plt.show()
