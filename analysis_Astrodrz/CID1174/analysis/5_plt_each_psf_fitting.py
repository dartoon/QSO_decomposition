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

#f0 = open("fit_result_each/each_PSF_fit_ps.txt","r")
#string0 = f0.read()
#print string0

f = open("fit_result_each/each_PSF_fit_qso.txt","r")
string = f.read()

#labels = re.findall(r"PSF\d+", string)
S_n_list = re.findall(r"n_sersic':(.*?),",string)
Re = re.findall(r"R_sersic':(.*?),",string)
host_flux_ratio = re.findall(r"host_flux_ratio_percent':(.*?)}",string)
Chisq = re.findall(r"redu_Chisq':(.*?),",string)

S_n_list = [float(value) for value in S_n_list]
Re = [float(i) for i in Re]
host_flux_ratio = [float(i) for i in host_flux_ratio]
Chisq = [float(i) for i in Chisq]

from filter_info import filt_info
import pickle
PSFs_dict = {}
for key in filt_info.keys():
    PSFs, _=pickle.load(open('../../{0}/analysis/{0}_PSFs_QSO'.format(key),'rb'))
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

for i in range(len(PSF_id)):
    print PSF_id[i], Chisq[i], host_flux_ratio[i], Re[i], S_n_list[i], filter_list[i]