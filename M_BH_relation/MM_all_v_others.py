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
import matplotlib as mpl
mpl.rc('image', cmap='jet')
import sys
sys.path.insert(0,'../py_tools')

#The mid-z data
f0 ='data/SS13_MM.txt'
ss = np.loadtxt(f0)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
f1 ='data/B11_MM.txt'
b11 = np.loadtxt(f1)[:,1:]  #0 redshift; 1 M*; 2 BH mass;

#%%
#==============================================================================
# My new inference
#==============================================================================
from load_result import load_host_p, load_MBH, load_err, load_re
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
Mstar = load_host_p(ID)[1]
MBs, La, FWHMa = load_MBH(ID,MB_ID,if_reportHb=0, return_line=1)
Mstar_err = load_err(prop = 'Mstar', ID=ID)
yerr_highz = [(Mstar_err[:,0]**2+0.4**2)**0.5, (Mstar_err[:,1]**2+0.4**2)**0.5]
Re_results = np.asarray(load_re(ID))
n_results = np.asarray(load_n(ID))
#%%Plot sig_MBH as MBH
########## input local data ####
host=plt.figure(figsize=(13,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label
from local_MM_v_others import *

#select_eff = np.loadtxt('../comparsion_Simulation_CANDEL/Select_effect/Andreas/hst_study_bias_logm.txt')
#plt.scatter(select_eff[:,0],select_eff[:,1],c='blue',
#            s=100   ,marker=".",zorder=100)

plt.scatter(MBs,MBs-(m_ml*Mstar+b_ml),c='tomato',
            s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
plt.errorbar(MBs,MBs-(m_ml*Mstar+b_ml),
             yerr= yerr_highz,
             color='tomato',ecolor='orange', fmt='.',markersize=1)    

plt.xlabel("log$(M_{BH}/M_{\odot})$",fontsize=35)
new_sample = mlines.Line2D([], [], color='tomato', ls='', marker='*', markersize=20,markeredgecolor='k')

plt.yticks(np.arange(-5.5,6,0.5))
plt.axis([7,9.6,-2.0,3.5])
plt.ylabel("$\Delta$log$M_{BH}$ (vs $M_*$)",fontsize=35)
plt.grid()
plt.tick_params(labelsize=25)

plt.legend([Bkc, Hkc, new_sample],[
'Local by Bennert+11',\
"Local by H&R",
"This work"
],scatterpoints=1,numpoints=1,loc=2,prop={'size':23},ncol=3,handletextpad=0)
#plt.savefig("offset_vs_MBH.pdf")
plt.show()

"""
#%%Plot sig_MBH as FWHM(Ha)
host=plt.figure(figsize=(13,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label
xl = np.linspace(-100, 100, 100)
plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)
plt.scatter(La,MBs-(m_ml*Mstar+b_ml),c='tomato',
            s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
plt.errorbar(La,MBs-(m_ml*Mstar+b_ml),
             yerr= yerr_highz,
             color='tomato',ecolor='orange', fmt='.',markersize=1)    

plt.xlabel("L_Halpha",fontsize=35)
plt.axis([int(La.min()),int(La.max())+1,-2.0,3.5])
plt.ylabel("$\Delta$log$M_{BH}$ (vs $M_*$)",fontsize=35)
plt.grid()
plt.tick_params(labelsize=25)
plt.show()


#%%Plot sig_MBH as FWHM(Ha)
host=plt.figure(figsize=(13,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label
xl = np.linspace(1500, 6000, 100)
plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)
plt.scatter(FWHMa,MBs-(m_ml*Mstar+b_ml),c='tomato',
            s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
plt.errorbar(FWHMa,MBs-(m_ml*Mstar+b_ml),
             yerr= yerr_highz,
             color='tomato',ecolor='orange', fmt='.',markersize=1)    

plt.xlabel("FWHM_Halpha",fontsize=35)
plt.axis([1800,5500,-2.0,3.5])
plt.ylabel("$\Delta$log$M_{BH}$ (vs $M_*$)",fontsize=35)
plt.grid()
plt.tick_params(labelsize=25)
plt.show()

#%%Plot sig_MBH as FWHM(Ha)
host=plt.figure(figsize=(13,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label
xl = np.linspace(1500, 6000, 100)
plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)
plt.scatter(FWHMa,MBs-(m_ml*Mstar+b_ml),c='tomato',
            s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
plt.errorbar(FWHMa,MBs-(m_ml*Mstar+b_ml),
             yerr= yerr_highz,
             color='tomato',ecolor='orange', fmt='.',markersize=1)    

plt.xlabel("FWHM_Halpha",fontsize=35)
plt.axis([1800,5500,-2.0,3.5])
plt.ylabel("$\Delta$log$M_{BH}$ (vs $M_*$)",fontsize=35)
plt.grid()
plt.tick_params(labelsize=25)
plt.show()

#%%
tab_list = ID
ext_ID = {'XID2202':'LID1622', 'XID2138':'LID1820', 'XID2396':'LID1878', 'CDFS-321':'ECDFS-321'}
f_mbh = open("../M_BH_relation/fmos_MBH_table","r")
with f_mbh as g:
    lines = g.readlines()
porp_list = lines[0].replace('#','').split(' ')
samples = [lines[i].split(' ') for i in range(1,len(lines))]
outliers = ['CDFS-1', 'SXDS-X1136', 'SXDS-X763', 'CDFS-724']
#outliers = ['CDFS-1', 'SXDS-X763', 'CDFS-724']
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
#for i in range(len(tab_list)):
#    print tab_list[i], round(MBH[i],3), round(logEdd[i], 3)
host=plt.figure(figsize=(13,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label
xl = np.linspace(-100, 100, 100)
plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)
plt.scatter(logEdd,MBs-(m_ml*Mstar+b_ml),c='tomato',
            s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
plt.errorbar(logEdd,MBs-(m_ml*Mstar+b_ml),
             yerr= yerr_highz,
             color='tomato',ecolor='orange', fmt='.',markersize=1)    

plt.xlabel("log($L_{bol}/L_{Edd}$)",fontsize=35)
plt.yticks(np.arange(-5.5,6,0.5))
plt.xlim([-1.5,0])
plt.ylim([-2.0,3.5])
plt.ylabel("$\Delta$log$M_{BH}$ (vs $M_*$)",fontsize=35)
plt.grid()
plt.tick_params(labelsize=25)
#plt.savefig("offset_vs_Edd.pdf")
plt.show()   

#%%Plot sig_MBH color in 2-D
host=plt.figure(figsize=(13,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label
xl = np.linspace(1500, 6000, 100)
plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)

import copy, matplotlib
my_cmap = copy.copy(matplotlib.cm.get_cmap('Oranges',10)) # copy the default cmap
plt.scatter(MBs,logEdd,c=MBs-(m_ml*Mstar+b_ml),
            s=580,marker="*",zorder=300, vmin=-0.5, vmax=1., edgecolors='k', cmap=my_cmap)
cl=plt.colorbar()
cl.set_label("$\Delta$log$M_{BH}$ (vs $M_*$)", size=20)
x_cline = np.linspace(7.1, 10)
y_line_0 = x_cline*0 - 1.5
plt.plot(x_cline, y_line_0,'k--')
y_line_1 = x_cline*0 + 0.0
plt.plot(x_cline, y_line_1,'k--')
y_cline = np.linspace(-20, 1, 10)
x_line = y_cline*0 + 8.6
plt.plot(x_line, y_cline,'r--')
x_cline = np.linspace(7.1, 8.5)
y_line_3 = -1.1*(x_cline-7.5) -0.5
plt.plot(x_cline, y_line_3,'k--')

plt.xlabel("log$(M_{BH}/M_{\odot})$",fontsize=35)
plt.ylabel("log$L_{bol}/L_{Edd}$", fontsize=35) 
#plt.yticks(np.arange(-5.5,6,0.5))
plt.axis([7,9.6,-1.9, 0.3])
plt.grid()
plt.tick_params(labelsize=25)
plt.show()

#%%Plot sig_MBH as Re_result
host=plt.figure(figsize=(13,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label
xl = np.linspace(0, 6000, 100)
plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)
plt.scatter(Re_results[:,0],MBs-(m_ml*Mstar+b_ml),c='tomato',
            s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
plt.errorbar(Re_results[:,0],MBs-(m_ml*Mstar+b_ml),
             yerr= yerr_highz, xerr=Re_results[:,1],
             color='tomato',ecolor='orange', fmt='.',markersize=1)    

plt.xlabel("Reff (arcsec)",fontsize=35)
plt.axis([0.05,1.0,-2.0,3.5])
plt.ylabel("$\Delta$log$M_{BH}$ (vs $M_*$)",fontsize=35)
plt.grid()
plt.tick_params(labelsize=25)
plt.show()

#%%Plot sig_MBH as Re_result(Ha)
host=plt.figure(figsize=(13,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label
xl = np.linspace(0, 6000, 100)
plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)
plt.scatter(n_results[:,0],MBs-(m_ml*Mstar+b_ml),c='tomato',
            s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
plt.errorbar(n_results[:,0],MBs-(m_ml*Mstar+b_ml),
             yerr= yerr_highz, xerr=n_results[:,1],
             color='tomato',ecolor='orange', fmt='.',markersize=1)    

plt.xlabel("Sersic",fontsize=35)
plt.axis([0.0,5.0,-2.0,3.5])
plt.ylabel("$\Delta$log$M_{BH}$ (vs $M_*$)",fontsize=35)
plt.grid()
plt.tick_params(labelsize=25)
plt.show()
"""