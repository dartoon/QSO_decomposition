#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 09:40:52 2018

@author: Dartoon

Compare the histogram of the fitting with different approach.
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

def load_result(ID, count_rank=8, sort_PSF_by='fit_result_each'):
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
    flux_dict, FWHM_dict, locs_dict, filter_dict, id_stars_dict=pickle.load(open('PSFs_lib_dict','rb'))
    for j in range(len(ID)):
        filt = filt_info[ID[j]]
        f = open("{0}/analysis/fit_result_each/each_PSF_fit_qso.txt".format(ID[j]),"r")
        if sort_PSF_by == 'fit_result_each':
            f_ps = open("{0}/analysis/fit_result_each/each_PSF_fit_qso.txt".format(ID[j]),"r")
        elif sort_PSF_by == 'fit_ps_each':
            f_ps = open("{0}/analysis/fit_ps_each/each_PSF_fit_ps.txt".format(ID[j]),"r")
        string = f.read()
        string_ps = f_ps.read()
        PSF_id = re.findall(r"by PSF(.*?):",string)
        S_n_list = re.findall(r"n_sersic':(.*?),",string)
        Re = re.findall(r"R_sersic':(.*?),",string)
        host_flux_ratio = re.findall(r"host_flux_ratio_percent':(.*?)}",string)
        if sort_PSF_by == 'fit_result_each':
            Chisq = re.findall(r"redu_Chisq':(.*?),",string)
        elif sort_PSF_by == 'fit_ps_each':
            Chisq = re.findall(r"redu_Chisq':(.*?)}",string_ps)
        re.findall(r"redu_Chisq':(.*?),",string)
        QSO_amp = re.findall(r"QSO_amp':(.*?),",string)
        host_amp = re.findall(r"host_amp':(.*?),",string)
        
        S_n_list = [float(value) for value in S_n_list]
        Re = [float(i) for i in Re]
        host_flux_ratio = [float(i) for i in host_flux_ratio]
        Chisq = [float(i) for i in Chisq]
        QSO_amp = [float(i) for i in QSO_amp]
        host_amp =  [float(i) for i in host_amp]
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
        if filt == "F140w":
            zp = 26.4524
        elif filt == "F125w":
            zp = 26.2303
        amp_min, amp, amp_max = weighted_host_flux - rms_host_flux, weighted_host_flux, weighted_host_flux + rms_host_flux
        mag_min =-2.5*log10(amp_min) + zp
        mag =-2.5*log10(amp) + zp
        mag_max =-2.5*log10(amp_max) + zp
    #    print "mag for", ID[i], "{0}~{1}~{2}".format(round(mag_min,3),round(mag,3),round(mag_max,3))
        mag_result.append(mag)
        mag_lh_result.append([mag_min-mag, mag-mag_max])
#        print "ID:{0}, filt:{1}".format(ID[j], filt[0])
#        print "the best PSFs for {0} is ID ".format(ID[j]), sort_Chisq[:8]
#        print "Weighted: HOST_Ratio, Reff, Sersic_n, total_flux, host_flux: {0}+-{1} {2}+-{3} {4}+-{5} {6}+-{7} {8}+-{9}".format(\
#        round(weighted_host_ratio,3), round(rms_host_ratio,3),\
#        round(weighted_Re,3), round(rms_Re,3), round(weighted_index,3),round(rms_index,3),\
#        round(weighted_total_flux,3), round(rms_total_flux,3),\
#        round(weighted_host_flux,3), round(rms_host_flux,3))
        best_PSF_id.append(sort_Chisq[:8])
        ratio_results.append([weighted_host_ratio, rms_host_ratio])
        Re_results.append([weighted_Re, rms_Re])
        n_results.append([weighted_index,rms_index])
        host_amp_results.append([weighted_host_flux, rms_host_flux])
        total_flux_results.append([weighted_total_flux, rms_total_flux])
        chisq_list.append(Chisq_best)
        inf_list.append(inf_alp)  
    return ratio_results, Re_results, n_results, host_amp_results, total_flux_results, mag_result, mag_lh_result


ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724'\
]

ratio_results_0 = np.asarray(load_result(ID, count_rank=5)[0])
ratio_results_1 = np.asarray(load_result(ID, count_rank=8)[0])
ratio_results_2 = np.asarray(load_result(ID, count_rank=10)[0])
ratio_results_3 = np.asarray(load_result(ID, count_rank=8, sort_PSF_by = 'fit_ps_each')[0])


plt.figure(figsize=(10,6))
common_params = dict(label=('Rank5','Rank8','Rank10','fit_PS'))
plt.hist((ratio_results_0[:,0], ratio_results_1[:,0], ratio_results_2[:,0], ratio_results_3[:,0]), **common_params)
plt.legend()
plt.show()



#print ratio_results