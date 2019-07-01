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

def load_flux(ID, flt = 'WFC3', folder='../',result_folder = 'analysis'):
    '''
    Load the name as a list based on the ID list 
    
    Parameter
    --------
        list: The list of the file name
        
    Return
    --------
        A list of the zs
    '''
    _, _, _, host_amp_results, _, _, _, _ = load_result(ID, flt = flt, folder=folder,result_folder=result_folder)
    return host_amp_results

from dmag import pass_dmag
from scipy.integrate import quad
h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)

def load_host_p(ID, folder='../', temp='1Gyrs', dm = 0):
    '''
    Return the host properties including
    1. host LR, 2. M_star, 3. Host Rest frame R absolute mags
    '''
    zs = np.asarray(load_zs(ID))
    mags = np.array(load_mag(ID, folder = folder)[0])
    from dmag import k_corr_R
    from filter_info import filt_info
    dm_k_R = []
    for i in range(len(zs)):
        dm_k_R.append(k_corr_R(zs[i],filt_info[ID[i]], galaxy_age = temp))
    dm_k_R = np.asarray(dm_k_R) # Get the k-correction for each target as an array
    dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #in pc
    host_Mags = mags -5*(np.log10(dl)-1) + dm_k_R + dm*pass_dmag(zs) # This is in AB system
    host_LR = 10 ** (0.4*(4.61-host_Mags)) #LR in AB
    host_Mags_Vega = host_Mags - 0.21  # Transfer to Vega system
    host_LR_Vega = 10 ** (0.4*(4.43-host_Mags_Vega)) #LR in Vega
    Mstar = np.log10(host_LR_Vega * 0.54 * 0.684 * 1.4191)  
    return np.log10(host_LR), Mstar, host_Mags

def load_MBH(ID, MB_ID, if_reportHb=0,folder = '..', return_line = 0):
    #ID change:
    #XID2202 to LID1622
    #XID2138 to LID1820
    #XID2396 to LID1878
    #CDFS321 to ECDFS321
    f = open(folder+"/M_BH_relation/fmos_MBH_table","r")
    with f as g:
        lines = g.readlines()
#    porp_list = lines[0].replace('#','').split(' ')
    samples = [lines[i].split(' ') for i in range(1,len(lines))]
    ID_ser_dic =  {}
    
    for j in range(len(MB_ID)):
        count = 0
        for i in range(len(samples)):
            if samples[i][1] == MB_ID[j]:
                ID_ser_dic.update({ID[j]:i})
                count += 1
        if count == 0:
            ID_ser_dic.update({ID[j]: -99})
    MBs = []
    LogLa_list, FWHMa_list =[], []
#    CDFS_FWHMa = {'CDFS-1': 5449.4022,'CDFS-229': 2254.0105, 'CDFS-724': 3351.852239}  # by Hyewon
#    CDFS_logLHadr = {'CDFS-1': 43.08,'CDFS-229': 43.30, 'CDFS-724': 42.561413} # by Hyewon
    CDFS_FWHMa = {'CDFS-1': 2000,'CDFS-229': 2190, 'CDFS-724': 2541}  # by Malte
    CDFS_logLHadr = {'CDFS-1': 43.02,'CDFS-229': 43.60, 'CDFS-724': 42.95} # by Malte    
    for tar_in in range(len(ID)):       
        t_name = ID[tar_in]
        ser = ID_ser_dic[t_name]
    #    print ser
        if ser!=-99 and float(samples[ser][10]) != 0:
            FWHM_a = float(samples[ser][8])
            logLHadr = float(samples[ser][6])
            cal_logMa = 6.71+0.48*(logLHadr-42)+2.12*np.log10(FWHM_a/1000)  # as used in Andreas
#            cal_logMa = 6.301+0.55*(logLHadr-42)+2.06*np.log10(FWHM_a/1000)  # as used in GREENE & HO 2005
#            cal_logMa = 6.459+0.55*(logLHadr-42)+2.*np.log10(FWHM_a/1000)  # as used in H0liCOW 7 and McGill
            mbh = cal_logMa
            MBs.append(mbh)
            LogLa_list.append(logLHadr)
            FWHMa_list.append(FWHM_a)
        if ser!=-99 and float(samples[ser][21]) != 0 and if_reportHb==1:
#            if if_reportHb ==1:
            print "use Hb for", ID[tar_in]
            FWHM_b = float(samples[ser][19])
            logL5100dr = float(samples[ser][16])
            cal_logMb = 6.91+0.5*(logL5100dr-44)+2.*np.log10(FWHM_b/1000)  # as used in Andreas
    #        cal_logMb = 6.882+0.518*(logL5100dr-44)+2.*np.log10(FWHM_b/1000)        # calibrated in H0liCOW 7
            mbh = (cal_logMa + cal_logMb)/2
    #        print mbh
            MBs[-1] = mbh
        if ser==-99 and t_name in ['CDFS-1','CDFS-229', 'CDFS-724']:
            FWHM_a = float(CDFS_FWHMa[t_name])
            logLHadr = float(CDFS_logLHadr[t_name])
            cal_logMa = 6.71+0.48*(logLHadr-42)+2.12*np.log10(FWHM_a/1000)  # as used in Andreas
            MBs.append(cal_logMa)
            LogLa_list.append(logLHadr)
            FWHMa_list.append(FWHM_a)            
        elif ser==-99:
            MBs.append(-99)
            LogLa_list.append(-99)
            FWHMa_list.append(-99)    
    if return_line ==1:
        return np.asarray(MBs), np.asarray(LogLa_list), np.asarray(FWHMa_list)
    else:
        return np.asarray(MBs)


LR_error = {'CID1174': [-0.11, +0.15],
'CID1281': [-0.11, +0.15],
'CID206': [-0.23, +0.52],
'CID216': [-0.03, +0.03],
'CID237': [-0.09, +0.11],
'CID255': [-0.11, +0.15],
'CID3242': [-0.11, +0.14],
'CID3570': [-0.02, +0.02],
'CID452': [-0.03, +0.03],
'CID454': [-0.04, +0.04],
'CID50': [-0.19, +0.35],
'CID543': [-0.12, +0.16],
'CID597': [-0.15, +0.22],
'CID607': [-0.15, +0.23],
'CID70': [-0.1, +0.12],
'LID1273': [-0.07, +0.09],
'LID1538': [-0.08, +0.09],
'LID360': [-0.05, +0.06],
'XID2138': [-0.06, +0.07],
'XID2202': [-0.1, +0.12],
'XID2396': [-0.16, +0.26],
'CDFS-1': [-0.12, +0.16],
'CDFS-229': [-0.05, +0.06],
'CDFS-321': [-0.17, +0.28],
'CDFS-724': [-0.15, +0.23],
'ECDFS-358': [-0.1, +0.12],
'SXDS-X1136': [-0.08, +0.09],
'SXDS-X50': [-0.09, +0.11],
'SXDS-X717': [-0.06, +0.07],
'SXDS-X735': [-0.1, +0.13],
'SXDS-X763': [-0.22, +0.47],
'SXDS-X969': [-0.14, +0.21]}

Mstar_error = {'CID1174': [-0.15, +0.18],
'CID1281': [-0.15, +0.18],
'CID206': [-0.25, +0.53],
'CID216': [-0.1, +0.1],
'CID237': [-0.13, +0.14],
'CID255': [-0.15, +0.18],
'CID3242': [-0.15, +0.17],
'CID3570': [-0.1, +0.1],
'CID452': [-0.1, +0.1],
'CID454': [-0.1, +0.1],
'CID50': [-0.21, +0.36],
'CID543': [-0.15, +0.19],
'CID597': [-0.18, +0.24],
'CID607': [-0.18, +0.25],
'CID70': [-0.14, +0.16],
'LID1273': [-0.12, +0.13],
'LID1538': [-0.12, +0.13],
'LID360': [-0.11, +0.11],
'XID2138': [-0.12, +0.12],
'XID2202': [-0.14, +0.16],
'XID2396': [-0.19, +0.28],
'CDFS-1': [-0.15, +0.19],
'CDFS-229': [-0.11, +0.12],
'CDFS-321': [-0.2, +0.3],
'CDFS-724': [-0.18, +0.25],
'ECDFS-358': [-0.14, +0.16],
'SXDS-X1136': [-0.13, +0.14],
'SXDS-X50': [-0.13, +0.15],
'SXDS-X717': [-0.12, +0.12],
'SXDS-X735': [-0.14, +0.17],
'SXDS-X763': [-0.24, +0.48],
'SXDS-X969': [-0.17, +0.23]}

def load_err(prop, ID):
    if prop=='Mstar':
        un_dict = Mstar_error
    elif prop == 'LR':
        un_dict = LR_error
    errs = [un_dict[ID[i]] for i in range(len(ID))]
    return np.array(errs)
    
def load_Lbol(ID, folder='../'):
    ext_ID = {'XID2202':'LID1622', 'XID2138':'LID1820', 'XID2396':'LID1878', 'CDFS-321':'ECDFS-321'}
    f_mbh = open(folder+"M_BH_relation/fmos_MBH_table","r")
    with f_mbh as g:
        lines = g.readlines()
    samples = [lines[i].split(' ') for i in range(1,len(lines))]
    #outliers = ['CDFS-1', 'SXDS-X1136', 'SXDS-X763', 'CDFS-724']
    ID_ser_dic =  {}
    import copy
    tab_sub_list = copy.deepcopy(ID)
    for i in range(len(tab_sub_list)):
        if tab_sub_list[i] in ext_ID.keys():
            tab_sub_list[i] = ext_ID[tab_sub_list[i]]
    for j in range(len(ID)):
        count = 0
        for i in range(len(samples)):
            if samples[i][1] == tab_sub_list[j]:
                ID_ser_dic.update({ID[j]:i})
                count += 1
        if count == 0:
            ID_ser_dic.update({ID[j]: -99})
    samples_Lbol = np.loadtxt(folder+"/M_BH_relation/fmos_BH_Lbol")
    MB_Lbol_info = []
    CDFS_Lbol = {'CDFS-1': 45.89,'CDFS-229': 45.68, 'CDFS-724': 44.95}
    for tar_in in range(len(ID)):       
        t_name = ID[tar_in]
        ser = ID_ser_dic[t_name]
        if ser!=-99 and float(samples_Lbol[ser, 1]) != 0:
            logLbol = float(samples_Lbol[ser, 1])
            MB_Lbol_info.append(logLbol)
        elif ser==-99 and t_name in CDFS_Lbol.keys():
            logLbol = float(CDFS_Lbol[t_name])
            MB_Lbol_info.append(logLbol)      
    return MB_Lbol_info
#ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
#'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
#MB_ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'ECDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','LID1820','LID1622',\
#'LID1878', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
#MBs = load_MBH(ID, MB_ID)

##test:
#ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
#'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
#zs_ = load_zs(ID)
#mag, mag_lh_result = load_mag(ID)
#errs = load_err(prop='Mstar', ID=ID)