#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 21:00:09 2019

@author: Dartoon

Compare the fitting before and after deep
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib as mat
import matplotlib.lines as mlines
from matplotlib import colors
mat.rcParams['font.family'] = 'STIXGeneral'
import matplotlib as mpl
mpl.rc('image', cmap='jet')

import sys
sys.path.insert(0,'../../py_tools')
from dmag import pass_dmag

#==============================================================================
# My new inference
#==============================================================================
from load_result import load_host_ratio

ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
mags_before = np.array(load_host_ratio(ID, folder = '../',result_folder='analysis'))[:,0]
mags_after = np.array(load_host_ratio(ID, folder = '../',result_folder='deep_analysis'))[:,0]

import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'

plt.figure(figsize=(10,6))
common_params = dict(label=('Previous analysis','After upgrading the fitting'))
plt.hist((mags_before, mags_after), **common_params)
plt.xlabel('Host flux ratio (%)',fontsize=25)
plt.ylabel('#', fontsize=25)
plt.tick_params(labelsize=25)
plt.legend(fontsize=15)
#plt.savefig('hist_compare.pdf')
plt.show()

