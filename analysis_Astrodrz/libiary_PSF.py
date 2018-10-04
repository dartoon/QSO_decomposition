#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 16:45:01 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle

filt_info = {'CDFS-1': 'F140w', 'CDFS-229': 'F125w', 'CDFS-724': 'F125w',\
'CID1174': 'F140w', 'CID206': 'F140w', 'CID216': 'F140w', 'CID3242': 'F140w',\
'CID3570': 'F125w', 'CID452': 'F125w', 'CID454': 'F140w', 'CID50': 'F125w',\
'CID607': 'F125w', 'CID70': 'F140w', 'LID1273': 'F140w', 'LID360': 'F140w',\
'XID2138': 'F140w', 'XID2202': 'F140w', 'XID2396': 'F140w', 'ECDFS-358': 'F140w',\
'SXDS-X1136': 'F125w', 'SXDS-X50': 'F125w', 'SXDS-X735': 'F140w',\
'CID543': 'F125w', 'LID1538': 'F140w', 'CID237': 'F140w', 'SXDS-X717': 'F125w',\
'SXDS-X763': 'F125w', 'SXDS-X969': 'F140w', 'CDFS-321': 'F140w'}


PSFs_dict = {}
QSOs_dict = {}
for key in filt_info.keys():
    ID = key
    PSFs, QSOs=pickle.load(open('{0}/analysis/{0}_PSFs_QSO'.format(ID),'rb'))
    PSFs_dict.update({'{0}'.format(ID):PSFs})
    QSOs_dict.update({'{0}'.format(ID):QSOs})

PSF_list = []
PSF_msk_list = []
QSO_list = []
for ID in PSFs_dict.keys():
    psfs_dict = PSFs_dict[ID]
    psfs = [psfs_dict[i][0] for i in range(len(psfs_dict))]
    psfs_msk = [psfs_dict[i][3] for i in range(len(psfs_dict))]
    PSF_list += psfs
    PSF_msk_list += psfs_msk
    
    qso = QSOs_dict[ID]
    QSO_list += [qso[0]]
#    if_ind += 
    
PSFs_fluxs=[]
for i in range(len(PSF_list)):
    flux = np.sum(PSF_list[i] * PSF_msk_list[i])
#    if flux< 1000:
    PSFs_fluxs.append([flux])
    
plt.hist(np.asarray(PSFs_fluxs))
plt.show()

QSOs_flux = [np.sum(QSO_list[i]) for i in range(len(QSO_list))]
plt.hist(np.asarray(QSOs_flux))
plt.show()

import sys
sys.path.insert(0,'../py_tools')
from flux_profile import profiles_compare
prf_list = PSF_list
scal_list = np.ones(len(PSF_list))
prf_name_list = ["PSF{0}".format(i) for i in range(len(PSF_list))]
fig_pro_compare = profiles_compare(prf_list, scal_list, prf_name_list=prf_name_list,norm_pix = 5.0,
                                   gridspace = 'log',if_annuli=True)