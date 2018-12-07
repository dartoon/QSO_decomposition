#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 10:29:28 2018

@author: Dartoon

Return the zs and mag as an array based on the input names
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import re
from filter_info import filt_info

redshift_info = {'CID50': 1.239, 'CID206': 1.483, 'CID216': 1.567, 'CID237': 1.618,\
'CID255': 1.664, 'CID452': 1.407, 'CID597': 1.272, 'CID607': 1.294, 'CID1174': 1.552,\
'CID1281': 1.445, 'CID3242': 1.532, 'CID3570': 1.244, 'XID2396': 1.600, 'CID70': 1.667,\
'LID1273': 1.617, 'XID2202': 1.516, 'CID454': 1.478, 'CID543': 1.301, 'XID2138': 1.551,
'LID1538': 1.527, 'CDFS-1': 1.630, 'CDFS-724': 1.337, 'LID360': 1.579, 'CDFS-321': 1.570,\
'CDFS-229': 1.326, 'ECDFS-358': 1.626, 'SXDS-X50': 1.411, 'SXDS-X717': 1.276,\
'SXDS-X735': 1.447, 'SXDS-X763': 1.412, 'SXDS-X969': 1.585, 'SXDS-X1136': 1.325}

def load_result(ID, count_rank=8, sort_PSF_by='fit_result_each', cam = 'WFC3'):
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
#    chisq_list, inf_list, best_PSF_id = [],[], []
    for j in range(len(ID)):
        if cam == 'WFC3':
            filt = filt_info[ID[j]]
            f = open("../analysis_Astrodrz/{0}/analysis/fit_result_each/each_PSF_fit_qso.txt".format(ID[j]),"r")
            if sort_PSF_by == 'fit_result_each':
                f_ps = open("../analysis_Astrodrz/{0}/analysis/fit_result_each/each_PSF_fit_qso.txt".format(ID[j]),"r")
            elif sort_PSF_by == 'fit_ps_each':
                f_ps = open("../analysis_Astrodrz/{0}/analysis/fit_ps_each/each_PSF_fit_ps.txt".format(ID[j]),"r")
        elif cam == 'ACS':
            filt = 'F814w'
            folder = "fit_result_each_fix"
            f = open("../analysis_ACS/{0}/{1}/each_PSF_fit_qso.txt".format(ID[j],folder),"r")
            f_ps = open("../analysis_ACS/{0}/{1}/each_PSF_fit_qso.txt".format(ID[j],folder),"r")
        string = f.read()
        string_ps = f_ps.read()
#        PSF_id = re.findall(r"by PSF(.*?):",string)
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
        elif filt == "F814w":
            zp = 25.94333
        amp_min, amp, amp_max = weighted_host_flux - rms_host_flux, weighted_host_flux, weighted_host_flux + rms_host_flux
        mag_min =-2.5*np.log10(amp_min) + zp
        mag =-2.5*np.log10(amp) + zp
        mag_max =-2.5*np.log10(amp_max) + zp
        mag_result.append(mag)
        mag_lh_result.append([mag_min-mag, mag-mag_max])
        ratio_results.append([weighted_host_ratio, rms_host_ratio])
        total_flux_results.append([weighted_total_flux, rms_total_flux])
        Re_results.append([weighted_Re, rms_Re])
        n_results.append([weighted_index,rms_index])
        host_amp_results.append([weighted_host_flux, rms_host_flux])
#        best_PSF_id.append(sort_Chisq[:count_n])
#        chisq_list.append(Chisq_best)
#        inf_list.append(inf_alp)  
    return mag_result, mag_lh_result, ratio_results, Re_results, n_results, host_amp_results, total_flux_results
#    return np.asarray(mag_result), np.asarray(mag_lh_result)

def load_zs(ID):
    '''
    Get a z list based on the ID list 
    
    Parameter
    --------
        list: The list of the file name
        
    Return
    --------
        A list of the zs
    '''
    zs_list = []
    for i in range(len(ID)):
        zs_list.append(redshift_info[ID[i]])
    return zs_list
        
    
def load_mag(ID, cam = 'WFC3'):
    '''
    Load the name as a list based on the ID list 
    
    Parameter
    --------
        list: The list of the file name
        
    Return
    --------
        A list of the zs
    '''
    mag, mag_lh_result,_ ,_ ,_ ,_ , _ = load_result(ID, cam = cam)
    return mag, mag_lh_result
    
    
#    
##test:
#ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
#'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
#zs_ = load_zs(ID)
#mag, mag_lh_result = load_mag(ID)