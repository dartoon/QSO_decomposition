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
plt.figure(figsize=(11.5,12))


import matplotlib as mpl
mpl.rc('image', cmap='jet')

import sys
sys.path.insert(0,'../py_tools')
from dmag import pass_dmag
########## input local data ####
#==============================================================================
# The seleting for dm and host_total and dmag are in this local
#==============================================================================
from local_MagMBH import *
#==============================================================================
#input Park's data 
#==============================================================================
################ bulge or total relaiton? #################
#host= input('with bulge or total relation??? 0 = bulge;   1= total:')
####### input Park data ####
#######in AB system, V band#######
if host == 0:
   f0 ='../M_BH_relation/data/Park_b'
if host == 1:
   f0 ='../M_BH_relation/data/Park_t'
pk = np.loadtxt(f0)  #0 redshift; 1 BH mass; 2 Lv,obs,solar=0.4*(4.83-M_v)
pk[:,2]=4.83-pk[:,2]/0.4  # tansfer from L to Magnitude 
pk[:,2]=pk[:,2]-0.46  # transfer from V to R; 0.46 is the mean value of the data from Taka
pk[:,2]=pk[:,2]+dm*pass_dmag(pk[:,0])  #evolution of stellar population makes the mag fainter.
#pk[:,2]=0.4*(4.61-pk[:,2])
plt.scatter(pk[:,2],pk[:,1],c='darkseagreen',marker="^",s=180,zorder=100, alpha = 0.7, edgecolors='white')
tx, ty = -24.5,7.0
plt.text(tx, ty, "intermediate\n  sample\nuncertainties",  fontsize=20)
plt.errorbar(tx+0.2,ty-0.05, xerr=0.2/0.4, yerr=0.4, color='darkseagreen',ecolor='black', fmt='^',zorder=-500,markersize=10)

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
host_n = np.array(load_n(ID, folder = '../'))[:,0]
lumi_s = load_host_p(ID, dm = dm)[0] #!!! This dm is important 
Mag_s = 4.61 - lumi_s/0.4
MBs = load_MBH(ID,MB_ID,if_reportHb=0   )
LR_err = load_err(prop = 'LR', ID=ID)
LR_err = LR_err/0.4
#plt.scatter(Mag_s,MBs,c='tomato',s=580,marker="*",zorder=100, edgecolors='k')
#plt.errorbar(Mag_s,MBs, xerr=[np.abs(LR_err)[:,0], np.abs(LR_err)[:,1]], yerr=0.4, color='blue',ecolor='orange', fmt='.',zorder=-500,markersize=1)
#==============================================================================
# The colorbar label setting up
#==============================================================================

if dm ==0:
    plt.title(r"M$_{\rm BH}-$Mag$_{\rm R}$ relation",fontsize=35)
#elif dm ==1:
#    plt.title("Evolution-corrected $M_{BH}-L_{host}$ relation",fontsize=35)

plt.xlabel(r"R band magnitude",fontsize=35)
plt.ylabel(r"log(M$_{\rm BH}$/M$_{\odot})$",fontsize=35)
#plt.xticks(np.arange(5,16.5,0.5))
plt.yticks(np.arange(6,12.0,0.5))
plt.axis([-26,-19.5,6.,10.25])
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)
plt.gca().invert_xaxis()

#lens = mlines.Line2D([], [], color='cyan', ls='', marker='o', markersize=9)
#peng = mlines.Line2D([], [], color='cyan', ls='', marker='s', markersize=9)
park = mlines.Line2D([], [], color='darkseagreen', ls='', marker='^', markersize=9)
#new_sample = mlines.Line2D([], [], color='tomato', ls='', marker='*', markersize=20,markeredgecolor='k')

plt.legend([Pkc, park],[
"Local AGNs by Bennert+10",
"AGNs of z > 0.3 by Park+2015",
],scatterpoints=1,numpoints=1,loc=1,prop={'size':24},ncol=1)
#else:
#    plt.legend([HE0435,RXJ1131,park,Pkc, new_sample],['HE0435',\
#    'RXJ1131',\
#    "AGNs from P15",
#    #'M$_{BH}$ calibrated from Mg$_{II}$',\
#    #"M$_{BH}$ calibrated from H${\\beta}$",'M$_{BH}$ calibrated from C$_{IV}$'
#    "local AGNs",
#    "our new samples"\
#    ],scatterpoints=1,numpoints=1,loc=2,prop={'size':24},ncol=2)
#plt.savefig('MBH-Mag_without_32QSO.pdf')
plt.show()
