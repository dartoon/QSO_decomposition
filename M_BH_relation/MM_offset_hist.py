#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 15:50:27 2018

@author: Dartoon
"""
import numpy as np
np.set_printoptions(precision=4)
import matplotlib.pyplot as plt
import matplotlib as mat
import matplotlib.lines as mlines
from matplotlib import colors
mat.rcParams['font.family'] = 'STIXGeneral'
host=plt.figure(figsize=(14.5,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label

import matplotlib as mpl
mpl.rc('image', cmap='jet')
import sys
sys.path.insert(0,'../py_tools')

########## input local data ####
#==============================================================================
# The seleting for dm and host_total and dmag are in this local
#==============================================================================
from local_MM_vz import *
plt.close()
#%%
#==============================================================================
# My new inference
#==============================================================================
from load_result import load_host_p, load_MBH, load_err
from load_result import load_zs, load_n
ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']
MB_ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'ECDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','LID1820','LID1622',\
'LID1878', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']
zs = np.asarray(load_zs(ID))
Mstar = load_host_p(ID)[1]
MBs = load_MBH(ID,MB_ID,if_reportHb=0)
#Mstar_err = load_err(prop = 'Mstar', ID=ID)
#yerr_highz = [(Mstar_err[:,0]**2+0.4**2)**0.5, (Mstar_err[:,1]**2+0.4**2)**0.5]
#%%
#local_offset = y -(m_ml*x+b_ml)
local_offset_bloc = bloc[:,3]- (m_ml*bloc[:,1]+b_ml)
local_offset_hloc = hloc[:,3]- (m_ml*hloc[:,1]+b_ml)

highz_offset = MBs-(m_ml*Mstar+b_ml)

Mstar_err = load_err(prop = 'Mstar', ID=ID)
yerr_highz = [(Mstar_err[:,0]**2+0.4**2)**0.5, (Mstar_err[:,1]**2+0.4**2)**0.5]
yerr_hz = (yerr_highz[0]+ yerr_highz[1])/2

yerr_hz = yerr_hz#[highz_offset<0.9]
highz_offset = highz_offset#[highz_offset<0.9]

plt.figure(figsize=(8,6))
#high0, x0, _ = plt.hist(local_offset,normed=True, histtype=u'step', linestyle=('dashed'),
#         label=('local sample'), linewidth = 2, color='gray')
plt.hist(local_offset_bloc,normed=True, histtype=u'step', linestyle=('dashed'),
         label=('Local by Bennert+11'), linewidth = 2, color='gray')
plt.hist(local_offset_hloc,normed=True, histtype=u'step', linestyle=(':'),
         label=("Local by H&R"), linewidth = 2, color='gray')
high1, x1, _ = plt.hist(highz_offset,normed=True, histtype=u'step',
         label=('High z sample'), linewidth = 2, color='black')

#x0_m = np.median(local_offset)
#high_m0 = high0[np.where(abs(x0_m-x0) == abs(x0_m-x0).min())[0][0]]
#x1_m = np.median(highz_offset)
x1_m = np.mean(highz_offset)
high_m1 = high1[np.where(abs(x1_m-x1) == abs(x1_m-x1).min())[0][0]]

x1_mean = np.mean(highz_offset)
x1_weight =  np.sum(np.asarray(highz_offset)*yerr_hz) / np.sum(yerr_hz)
xi_weight_err = np.sqrt(np.sum((np.asarray(highz_offset)-x1_weight)**2*yerr_hz) / np.sum(yerr_hz))
#plt.plot(np.linspace(0,high_m0)*0 , np.linspace(0,high_m0),
#         linewidth = 2,color='gray',linestyle=('dashed'),alpha=0.3)


plt.plot(np.linspace(0,high_m1)*0+x1_m , np.linspace(0,high_m1),
         linewidth = 2, color='black',alpha=0.3)
#plt.plot(np.linspace(0,high_m1)*0+x1_m-0.21 , np.linspace(0,high_m1)*2,
#         linewidth = 2, color='red',alpha=0.5)
plt.plot(np.linspace(0,high_m1)*0 + 0.21 , np.linspace(0,high_m1)*2,
         linewidth = 2, color='red',alpha=0.5)


plt.text((x1_m)*0.85+0.05, 0.15, '{0}'.format(round(x1_m,2)), color='gray',fontsize=18)
plt.text((x1_m)*0+0.02, 0.15, '0.21', color='red',fontsize=18,alpha=0.5)
#plt.text(np.median(x1_m)-0.21, 0.25, '{0}'.format(round(np.median(x1_m)-0.21,2)), color='red',fontsize=18)
plt.xlabel("$\Delta$log($M_{BH}$)",fontsize=27)
plt.ylabel("Number Density",fontsize=27)
plt.ylim([0, 1.8])
plt.tick_params(labelsize=20)
plt.legend(prop={'size':14})
#plt.savefig('hist_offset.pdf')
plt.show()

