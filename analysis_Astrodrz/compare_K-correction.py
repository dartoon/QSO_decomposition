#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:29:58 2018

@author: Dartoon

To do the K-correction and check compared to the ACS
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import re
from math import log10
from adjustText import adjust_text   # avoid the overlapping while ploting
import sys
sys.path.insert(0,'../py_tools')
from filter_info import filt_info, redshift_info
import pickle

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

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
    chisq_list, inf_list, best_PSF_id = [],[], []
    flux_dict, FWHM_dict, locs_dict, filter_dict, id_stars_dict=pickle.load(open('PSFs_lib_dict','rb'))
    for j in range(len(ID)):
        if cam == 'WFC3':
            filt = filt_info[ID[j]]
            f = open("{0}/analysis/fit_result_each/each_PSF_fit_qso.txt".format(ID[j]),"r")
            if sort_PSF_by == 'fit_result_each':
                f_ps = open("{0}/analysis/fit_result_each/each_PSF_fit_qso.txt".format(ID[j]),"r")
            elif sort_PSF_by == 'fit_ps_each':
                f_ps = open("{0}/analysis/fit_ps_each/each_PSF_fit_ps.txt".format(ID[j]),"r")
        elif cam == 'ACS':
            filt = 'F814w'
            folder = "fit_result_each_fix"
            f = open("../analysis_ACS/{0}/first_analysis/{1}/each_PSF_fit_qso.txt".format(ID[j],folder),"r")
            f_ps = open("../analysis_ACS/{0}/first_analysis/{1}/each_PSF_fit_qso.txt".format(ID[j],folder),"r")
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
#    return ratio_results, Re_results, n_results, host_amp_results, total_flux_results, np.asarray(mag_result), np.asarray(mag_lh_result)
    return np.asarray(mag_result), np.asarray(mag_lh_result)

from scipy.integrate import quad
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)
h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]

from dmag import k_corr_R
#ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
#'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724'\
#]
ID = ['CID1174','CID216', 'CID50','CID70','XID2138','CID3242',\
'LID1273','XID2202','CID206','CID543','LID1538','XID2396','CID452',\
'LID360','CID237','CID454','CID607','CID3570', 'CID597', 'CID1281','CID255']
mags_obs_IR, mag_IR_err = load_result(ID, count_rank=8, cam = 'WFC3')
mags_obs_UV, mag_UV_err = load_result(ID, count_rank=8, cam = 'ACS')
zs = np.array([redshift_info[ID[i]] for i in range(len(ID))])
dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #in pc
dm_k_IR, dm_k_UV = [],[]
galay_temp = '1Gyrs'
for i in range(len(ID)):
    dm_k_IR.append(k_corr_R(redshift_info[ID[i]],filt_info[ID[i]], galaxy_age = galay_temp))
    dm_k_UV.append(k_corr_R(redshift_info[ID[i]],'F814w', galaxy_age = galay_temp))
dm_k_IR = np.asarray(dm_k_IR)   # Get the k-correction for each target as an array
dm_k_UV = np.asarray(dm_k_UV)
mag_k_corrected_IR=mags_obs_IR-5*(np.log10(dl)-1) + dm_k_IR 
mag_k_corrected_UV=mags_obs_UV-5*(np.log10(dl)-1) + dm_k_UV 
#for i in range(len(ID)):
#    print ID[i], mag_k_corrected_IR[i], mag_k_corrected_UV[i]
#if needed lumi_s = 0.4*(4.61-mag_k_corrected)

plt.figure(figsize=(10, 10))
plt.errorbar(mag_k_corrected_IR, mag_k_corrected_UV,xerr= mag_IR_err.T, yerr=mag_UV_err.T,color='red', fmt='o',ecolor='gray' )
#texts = []
#for i in range(len(ID)):
#    texts.append(plt.text(mag_k_corrected_IR[i], mag_k_corrected_UV[i], ID[i], fontsize=17))
#adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))
x=np.linspace(-26,-20,15)
y = x
plt.plot(x,y,'b',linewidth=4.0,alpha=0.5)
#plt.title('Galaxy age: '+galay_temp, fontsize=24)
plt.xlabel('Absolute rest-frame R mag by WFC3',fontsize=30)
plt.ylabel('Absolute rest-frame R mag by ACS', fontsize=30)
val_min, val_max = np.min([mag_k_corrected_UV, mag_k_corrected_IR]), np.max([mag_k_corrected_UV, mag_k_corrected_IR])
plt.xlim(val_min-0.5, val_max+0.5)
plt.ylim(val_min-0.5, val_max+0.5)
plt.tick_params(labelsize=20)     
plt.savefig('comp_gtemp_{0}.pdf'.format(galay_temp))
plt.show()

#%%
#Plot the color and compare with the Tommaso's color value:
#Plot TT's figure:
from dmag import k_corr_R
plt.figure(figsize=(10, 8))
z_d = np.linspace(1,2,40)
k_c_140_5gy = k_corr_R(z_d, filt = 'F140w', galaxy_age = '5Gyrs')
k_c_814_5gy = k_corr_R(z_d, filt = 'F814w', galaxy_age = '5Gyrs')
k_c_125_5gy = k_corr_R(z_d, filt = 'F125w', galaxy_age = '5Gyrs')
k_c_140_1gy = k_corr_R(z_d, filt = 'F140w', galaxy_age = '1Gyrs')
k_c_814_1gy = k_corr_R(z_d, filt = 'F814w', galaxy_age = '1Gyrs')
k_c_125_1gy = k_corr_R(z_d, filt = 'F125w', galaxy_age = '1Gyrs')
plt.plot(z_d[z_d>1.44],(k_c_140_5gy-k_c_814_5gy)[z_d>1.44], label = '5gy', c='r')
plt.plot(z_d[z_d<1.44],(k_c_125_5gy-k_c_814_5gy)[z_d<1.44], c='r')
plt.plot(z_d[z_d>1.44],(k_c_140_1gy-k_c_814_1gy)[z_d>1.44], label = '1gy', c='b')
plt.plot(z_d[z_d<1.44],(k_c_125_1gy-k_c_814_1gy)[z_d<1.44], c='b')
plt.legend(fontsize=25)
texts = []
for i in range(len(ID)):
    if zs[i] > 1.44 :
        plt.plot(zs[i], mags_obs_UV[i] - mags_obs_IR[i], 'go')
    else:
        plt.plot(zs[i], mags_obs_UV[i] - mags_obs_IR[i], 'ko')
    texts.append(plt.text(zs[i], mags_obs_UV[i] - mags_obs_IR[i], ID[i], fontsize=17))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='green'))  
plt.xlabel('redshift',fontsize=35)
plt.ylabel('Observed Color', fontsize=35)
plt.tick_params(labelsize=20)    
plt.show()


    