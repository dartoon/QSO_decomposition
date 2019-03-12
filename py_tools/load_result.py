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
import pickle

redshift_info = {'CID50': 1.239, 'CID206': 1.483, 'CID216': 1.567, 'CID237': 1.618,\
'CID255': 1.664, 'CID452': 1.407, 'CID597': 1.272, 'CID607': 1.294, 'CID1174': 1.552,\
'CID1281': 1.445, 'CID3242': 1.532, 'CID3570': 1.244, 'XID2396': 1.600, 'CID70': 1.667,\
'LID1273': 1.617, 'XID2202': 1.516, 'CID454': 1.478, 'CID543': 1.301, 'XID2138': 1.551,
'LID1538': 1.527, 'CDFS-1': 1.630, 'CDFS-724': 1.337, 'LID360': 1.579, 'CDFS-321': 1.570,\
'CDFS-229': 1.326, 'ECDFS-358': 1.626, 'SXDS-X50': 1.411, 'SXDS-X717': 1.276,\
'SXDS-X735': 1.447, 'SXDS-X763': 1.412, 'SXDS-X969': 1.585, 'SXDS-X1136': 1.325}

def load_result(IDs, flt, count_rank=8, outputBHmag = 0, folder = '../',result_folder = 'analysis'):
    '''
    Read the fitting result by given the ID, 
    
    Parameter
    --------
        IDs: A list of ID
        b: The blash of blash
        sort_PSF_by: 'fit_result_each' or 'fit_ps_each'
    Return
    --------
        A sth sth
    '''
    ratio_results, Re_results, n_results, total_flux_results, host_amp_results, mag_result, mag_lh_result = [], [], [], [], [], [], []
    chisq_list, inf_list, best_PSF_id = [],[], []
#    flux_dict, FWHM_dict, locs_dict, filter_dict, id_stars_dict=pickle.load(open('{0}analysis_Astrodrz/PSFs_lib_dict'.format(folder),'rb'))
    for j in range(len(IDs)):
        if flt == 'WFC3':
            filt = filt_info[IDs[j]]
            f = open("{1}analysis_Astrodrz/{0}/{2}/fit_result_each/each_PSF_fit_qso.txt".format(IDs[j],folder,result_folder),"r")
        elif flt == 'ACS':
            filt = 'F814w'
            f = open("{1}analysis_ACS/{0}/first_analysis/fit_result_each_fix/each_PSF_fit_qso.txt".format(IDs[j],folder,result_folder),"r")
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
        mag_result.append(mag)
        mag_lh_result.append([mag_min-mag, mag-mag_max])
        
        bh_min, bh_amp, bh_max = weighted_bh_flux - rms_bh_flux, weighted_bh_flux, weighted_bh_flux + rms_bh_flux
        bh_mag_min =-2.5*np.log10(bh_min) + zp
        bh_mag =-2.5*np.log10(bh_amp) + zp
        bh_mag_max =-2.5*np.log10(bh_max) + zp
        bh_mag_result =bh_mag
        bh_mag_lh_result = [bh_mag_min-bh_mag, bh_mag-bh_mag_max]
    
        best_PSF_id.append([sort_Chisq[:8]])
        ratio_results.append([weighted_host_ratio, rms_host_ratio])
        Re_results.append([weighted_Re, rms_Re])
        n_results.append([weighted_index,rms_index])
        host_amp_results.append([weighted_host_flux, rms_host_flux])
        total_flux_results.append([weighted_total_flux, rms_total_flux])
        chisq_list.append([Chisq_best])
        inf_list.append([inf_alp])  
    if outputBHmag != 1 :
        return ratio_results, Re_results, n_results, host_amp_results, total_flux_results, mag_result, mag_lh_result, Chisq_best
    elif outputBHmag == 1 :
        return ratio_results, Re_results, n_results, host_amp_results, total_flux_results, mag_result, mag_lh_result, Chisq_best, bh_mag_result, bh_mag_lh_result
    

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
        
    
def load_mag(ID, flt = 'WFC3', folder='../',result_folder = 'analysis'):
    '''
    Load the name as a list based on the ID list 
    
    Parameter
    --------
        list: The list of the file name
        
    Return
    --------
        A list of the zs
    '''
    _, _, _, _, _, mag_result, mag_lh_result, _ = load_result(ID, flt = flt, folder=folder,result_folder=result_folder)
    return mag_result, mag_lh_result
    
        
def load_host_ratio(ID, flt = 'WFC3', folder='../',result_folder = 'analysis'):
    '''
    Load the name as a list based on the ID list 
    
    Parameter
    --------
        list: The list of the file name
        
    Return
    --------
        A list of the zs
    '''
    ratio_results, _, _, _, _, _, _, _ = load_result(ID, flt = flt, folder=folder,result_folder=result_folder)
    return ratio_results
   
def load_n(ID, folder='../',result_folder = 'analysis'):
    '''
    Load the name as a list based on the ID list 
    
    Parameter
    --------
        list: The list of the file name
        
    Return
    --------
        A list of the zs
    '''
    _, _, n_results, _, _, _, _, _ = load_result(ID, flt = 'WFC3', folder=folder,result_folder=result_folder)
    return n_results

def load_re(ID, folder='../',result_folder = 'analysis'):
    '''
    Load the name as a list based on the ID list 
    
    Parameter
    --------
        list: The list of the file name
        
    Return
    --------
        A list of the zs
    '''
    _, re_results, _, _, _, _, _, _ = load_result(ID, flt = 'WFC3', folder=folder,result_folder=result_folder)
    return re_results
    
    
##test:
#ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
#'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
#zs_ = load_zs(ID)
#mag, mag_lh_result = load_mag(ID)