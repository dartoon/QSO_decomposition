#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 09:25:32 2018

@author: Dartoon
"""

import sys
sys.path.insert(0,'../py_tools')
from filter_info import filt_info, redshift_info
import pickle, re
import numpy as np
from math import log10

def load_result(ID, flt, count_rank=8):
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
        f = open("../analysis_Astrodrz/{0}/analysis/fit_result_each/each_PSF_fit_qso.txt".format(ID),"r")
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
    mag_min =-2.5*log10(amp_min) + zp
    mag =-2.5*log10(amp) + zp
    mag_max =-2.5*log10(amp_max) + zp
    mag_result = mag 
    mag_lh_result = [mag_min-mag, mag-mag_max]
    best_PSF_id.append(sort_Chisq[:8])
    ratio_results = [weighted_host_ratio, rms_host_ratio]
    Re_results = [weighted_Re, rms_Re]
    n_results = [weighted_index,rms_index]
    host_amp_results = [weighted_host_flux, rms_host_flux]
    total_flux_results = [weighted_total_flux, rms_total_flux]
    chisq_list.append(Chisq_best)
    inf_list.append(inf_alp)  
    return ratio_results, Re_results, n_results, host_amp_results, total_flux_results, mag_result, mag_lh_result, Chisq_best

##Test loading
#print load_result('CID1174', 'WFC3')
#print load_result('CID1174', 'ACS')

##import os 
##dir_path = os.path.dirname(os.path.realpath(__file__))
#
WFC3_list = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724'\
]
    
ACS_IDs = ['CID1174','CID216', 'CID50','CID70','XID2138','CID3242',\
'LID1273','XID2202','CID206','CID3570','CID543','LID1538','XID2396','CID452',\
'LID360','CID237','CID454','CID607']     

tab_list = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 'CID3570', 'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 'CID70', 'LID1273', 'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396', 'CDFS-1', 'CDFS-229', 'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']
filt_list, ACS_list = [],[]
z_list = []
IR_result, UV_result = [], []
for i in range(len(tab_list)):
    if tab_list[i] in WFC3_list:   #Input the WFC3 filter
        filt_list.append(filt_info[tab_list[i]])
        z_list.append(redshift_info[tab_list[i]])
        IR_result.append(load_result(tab_list[i],'WFC3'))
    else:
        filt_list.append("xxx")
        z_list.append('xxx')
        IR_result.append('xxx')
    if tab_list[i] in ACS_IDs:
        ACS_list.append('ACS')
        UV_result.append(load_result(tab_list[i],'ACS'))
    else:
        ACS_list.append("xxx")
        UV_result.append('xxx')

for i in range(len(tab_list)):
    if IR_result[i]=='xxx':
        IR_result[i] = [np.asarray(IR_result[0][j])*0-99 for j in range(len(IR_result[0]))]
    if UV_result[i]=='xxx':
        UV_result[i] = [np.asarray(UV_result[0][j])*0-99 for j in range(len(UV_result[0]))]
#    print(tab_list[i], round(IR_result[i][-1],3),
#    "{0}%pm{1}%".format(round(IR_result[i][0][0],1), round(IR_result[i][0][1],1)),  #host flux ratio
#    "{0}pm{1}".format(round(IR_result[i][1][0],3), round(IR_result[i][1][1],3)),  #host Re
#    "{0}pm{1}".format(round(IR_result[i][2][0],3), round(IR_result[i][2][1],3)),  #host n
#    "{0}+{1}-{2}".format(round(IR_result[i][5],3), round(IR_result[i][6][0],3),round(IR_result[i][6][1],3)),  #host mag, l, h
#    UV_result[i][-1],
#    "{0}%pm{1}%".format(round(UV_result[i][0][0],1), round(UV_result[i][0][1],1)),  #host flux ratio
#    "{0}+{1}-{2}".format(round(UV_result[i][5],3), round(UV_result[i][6][0],3),round(UV_result[i][6][1],3))) #host mag, l, h

#XID2202 to LID1622
#XID2138 to LID1820
#XID2396 to LID1878
#CDFS-321 to ECDFS321
ext_ID = {'XID2202':'LID1622', 'XID2138':'LID1820', 'XID2396':'LID1878', 'CDFS-321':'ECDFS-321'}
MBs = []
f_mbh = open("../M_BH_relation/fmos_MBH_table","r")
with f_mbh as g:
    lines = g.readlines()
porp_list = lines[0].replace('#','').split(' ')
samples = [lines[i].split(' ') for i in range(1,len(lines))]
ID_ser_dic =  {}
import copy
tab_sub_list = copy.deepcopy(tab_list)
for i in range(len(tab_sub_list)):
    if tab_sub_list[i] in ext_ID.keys():
        tab_sub_list[i] = ext_ID[tab_sub_list[i]]
for j in range(len(tab_list)):
    count = 0
    for i in range(len(samples)):
        if samples[i][1] == tab_sub_list[j]:
            ID_ser_dic.update({tab_list[j]:i})
            count += 1
    if count == 0:
        ID_ser_dic.update({tab_list[j]: -99})
MB_info_a = []
MB_info_b = []
for tar_in in range(len(tab_list)):       
    t_name = tab_list[tar_in]
    ser = ID_ser_dic[t_name]
#    print ser
    if ser!=-99 and float(samples[ser][10]) != 0:
        FWMH_a = float(samples[ser][8])
        logLHadr = float(samples[ser][6])
#        cal_logMa = 6.71+0.48*(logLHadr-42)+2.12*np.log10(FWMH_a/1000)  # as used in Andreas
        cal_logMa = 6.459+0.55*(logLHadr-42)+2.*np.log10(FWMH_a/1000)  # as used in H0liCOW 7 and McGill
        MB_info_a.append([FWMH_a, logLHadr, cal_logMa])
    else:
        MB_info_a.append([-99,-99,-99])
    if ser!=-99 and float(samples[ser][21]) != 0:
        FWMH_b = float(samples[ser][19])
        logL5100dr = float(samples[ser][16])
#        cal_logMb = 6.91+0.5*(logL5100dr-44)+2.*np.log10(FWMH_b/1000)  # as used in Andreas
        cal_logMb = 6.882+0.518*(logL5100dr-44)+2.*np.log10(FWMH_b/1000)        # calibrated in H0liCOW 7
        MB_info_b.append([FWMH_b, logL5100dr, cal_logMb])
    else:
        MB_info_b.append([-99,-99,-99])
        
for i in range(len(tab_list)):
    print(tab_list[i], MB_info_a[i][0], MB_info_a[i][1], round(MB_info_a[i][2],3), MB_info_b[i][0],MB_info_b[i][1],round(MB_info_b[i][2],3))
    