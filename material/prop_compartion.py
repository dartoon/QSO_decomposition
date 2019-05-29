
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 16:11:41 2019

@author: Dartoon

Make the properties comparsion
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'../py_tools')
from load_result import load_MBH, load_host_ratio,load_host_p, load_zs,load_err

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

MBs = load_MBH(ID,MB_ID,if_reportHb=0)
host_ratio = np.asarray(load_host_ratio(ID))
zs = np.asarray(load_zs(ID))
Mstar = load_host_p(ID)[1]
Mstar_err = load_err(prop = 'Mstar', ID=ID)
#%%
tab_list = ID
ext_ID = {'XID2202':'LID1622', 'XID2138':'LID1820', 'XID2396':'LID1878', 'CDFS-321':'ECDFS-321'}
f_mbh = open("../M_BH_relation/fmos_MBH_table","r")
with f_mbh as g:
    lines = g.readlines()
porp_list = lines[0].replace('#','').split(' ')
samples = [lines[i].split(' ') for i in range(1,len(lines))]
#outliers = ['CDFS-1', 'SXDS-X1136', 'SXDS-X763', 'CDFS-724']
outliers = ['CDFS-1', 'SXDS-X763', 'CDFS-724']
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

samples_Lbol = np.loadtxt("../M_BH_relation/fmos_BH_Lbol")
MB_Lbol_info = []
CDFS_Lbol = {'CDFS-1': 45.89,'CDFS-229': 45.68, 'CDFS-724': 44.95}
for tar_in in range(len(tab_list)):       
    t_name = tab_list[tar_in]
    ser = ID_ser_dic[t_name]
    if ser!=-99 and float(samples_Lbol[ser, 1]) != 0:
        logLbol = float(samples_Lbol[ser, 1])
        MB_Lbol_info.append(logLbol)
    elif ser==-99 and t_name in CDFS_Lbol.keys():
        logLbol = float(CDFS_Lbol[t_name])
        MB_Lbol_info.append(logLbol)    
##        print t_name, MB_Lbol_info[-1]
#    if ser!=-99 and float(samples_Lbol[ser, 3]) != 0:
#        logLbol = float(samples_Lbol[ser, 3])
##        print "use Hb for", tab_list[tar_in], logLbol, MB_Lbol_info[-1]
#        MB_Lbol_info[-1] = (logLbol + MB_Lbol_info[-1])/2
        
MBH = load_MBH(tab_list,tab_sub_list, if_reportHb=0)   
print  "Calculate the Eddington ratio:"
logLedd = 38 + np.log10(1.2) + MBH
logEdd = MB_Lbol_info - logLedd
#%%Make the plot 
plt.figure(figsize=(8,6))
plt.errorbar(zs,host_ratio[:,0], yerr=host_ratio[:,1],fmt='.',color='gray',markersize=15)
plt.xlabel("z",fontsize=20)
plt.ylabel("host flux ratio",fontsize=20)
plt.grid(linestyle='--')
plt.tick_params(labelsize=15)
plt.show()

#%%Make the plot
plt.figure(figsize=(8,6))
plt.errorbar(MBs,host_ratio[:,0], xerr=0.4,fmt='.',color='gray',markersize=15)
plt.xlabel("MBH",fontsize=20)
plt.ylabel("host flux ratio",fontsize=20)
plt.grid(linestyle='--')
plt.tick_params(labelsize=15)
plt.show()

#%%Make the plot
plt.figure(figsize=(8,6))
plt.errorbar(logEdd,host_ratio[:,0], xerr=0,fmt='.',color='gray',markersize=15)
plt.xlabel("logEdd",fontsize=20)
plt.ylabel("host flux ratio",fontsize=20)
plt.grid(linestyle='--')
plt.tick_params(labelsize=15)
plt.show()

#%%Make the plot
plt.figure(figsize=(8,6))
plt.errorbar(zs,MBs, yerr=0.4,fmt='.',color='gray',markersize=15)
plt.xlabel("zs",fontsize=20)
plt.ylabel("MBs",fontsize=20)
plt.grid(linestyle='--')
plt.tick_params(labelsize=15)
plt.show()

#%%Make the plot
plt.figure(figsize=(8,6))
plt.errorbar(zs,Mstar, yerr=[np.abs(Mstar_err)[:,0], np.abs(Mstar_err)[:,1]],fmt='.',color='gray',markersize=15)
plt.xlabel("zs",fontsize=20)
plt.ylabel("Mstar",fontsize=20)
plt.grid(linestyle='--')
plt.tick_params(labelsize=15)
plt.show()