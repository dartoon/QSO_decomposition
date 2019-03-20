#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 21:46:04 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import re
import matplotlib

import os
path = os.getcwd()
ID = path.split('/')[-2]

import sys
sys.path.insert(0,'../../../py_tools')
from filter_info import filt_info
filt = filt_info[ID]

del filt_info[ID]

f = open("fit_result_each/each_PSF_fit_qso.txt","r")
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

import pickle
PSFs_dict = {}
for key in filt_info.keys():
    if filt_info[key] == filt:
        PSFs, _=pickle.load(open('../../{0}/analysis/{0}_PSFs_QSO'.format(key),'rb'))
        PSFs_dict.update({'{0}'.format(key):PSFs})

flux_dict, FWHM_dict, locs_dict, filter_dict, id_stars_dict=pickle.load(open('../../PSFs_lib_dict','rb'))

#for i in range(len(PSF_id)):
#    print PSF_id[i], Chisq[i], host_flux_ratio[i], Re[i], S_n_list[i], filter_list[i]

total_flux = np.asarray(QSO_amp) + np.asarray(host_amp)

sort_Chisq = np.argsort(np.asarray(Chisq))
count_n = 8
Chisq_best = Chisq[sort_Chisq[0]]
Chisq_last= Chisq[sort_Chisq[count_n-1]]
inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
weight = np.zeros(len(Chisq))
for i in sort_Chisq[:count_n]:
    weight[i] = np.exp(-1/2. * (Chisq[i]-Chisq_best)/(Chisq_best* inf_alp))
for i in sort_Chisq:
    print PSF_id[i], Chisq[i], weight[i], host_flux_ratio[i], Re[i], S_n_list[i], round(total_flux[i],3),  round(host_amp[i],3),\
round(flux_dict[PSF_id[i]],3), round(FWHM_dict[PSF_id[i]],3), "[{0},{1}]".format(int(round(locs_dict[PSF_id[i]][0])) , int(round(locs_dict[PSF_id[i]][1]))), round(id_stars_dict[PSF_id[i]],3)

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


print "Weighted: HOST_Ratio, Reff, Sersic_n, total_flux, host_flux: {0}+-{1} {2}+-{3} {4}+-{5} {6}+-{7} {8}+-{9}".format(\
round(weighted_host_ratio,3), round(rms_host_ratio,3),\
round(weighted_Re,3), round(rms_Re,3), round(weighted_index,3),round(rms_index,3),\
round(weighted_total_flux,3), round(rms_total_flux,3),\
round(weighted_host_flux,3), round(rms_host_flux,3))