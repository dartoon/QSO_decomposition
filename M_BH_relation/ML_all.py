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
from local_ML import *
#==============================================================================
#input Park's data 
#==============================================================================
################ bulge or total relaiton? #################
#host= input('with bulge or total relation??? 0 = bulge;   1= total:')
####### input Park data ####
#######in AB system, V band#######
if host == 0:
   f0 ='data/Park_b'
if host == 1:
   f0 ='data/Park_t'
pk = np.loadtxt(f0)  #0 redshift; 1 BH mass; 2 Lv,obs,solar=0.4*(4.83-M_v)
pk[:,2]=4.83-pk[:,2]/0.4  # tansfer from L to Magnitude 
pk[:,2]=pk[:,2]-0.46  # transfer from V to R; 0.46 is the mean value of the data from Taka
pk[:,2]=pk[:,2]+dm*pass_dmag(pk[:,0])  #evolution of stellar population makes the mag fainter.
pk[:,2]=0.4*(4.61-pk[:,2])
plt.scatter(pk[:,2],pk[:,1],c='darkseagreen',marker="^",s=80,zorder=100, alpha = 0.7, edgecolors='white')
tx, ty = 9.05,9.1
plt.text(tx, ty, "intermediate\n  sample\nuncertainties",  fontsize=20)
plt.errorbar(tx+0.2,ty-0.05, xerr=0.2, yerr=0.4, color='darkseagreen',ecolor='black', fmt='^',zorder=-500,markersize=10)

##==============================================================================
## input Peng's data
##==============================================================================
#inp_peng = 0
#from cal_peng import *
#
#Mg_cali=vec_Mgcal(L_Mg,FWHM_Mg)
#Mg[:,3]=Mg[:,3]+0.21+dm*pass_dmag(Mg[:,1])
#Hb_cali=vec_Hcal(L_hb,FWHM_hb)
#H[:,3]=H[:,3]+0.21+dm*pass_dmag(H[:,1])   #change Mg from Vega to AB
#C_cali=vec_Ccal(L_c,FWHM_c)
#C[:,3]=C[:,3]+0.21+dm*pass_dmag(C[:,1])   #change Mg from Vega to AB
###################plot#############################
################plot the data in Peng################
##For lensed samples
#C[:,3],Mg[:,3],H[:,3]=0.4*(4.61-C[:,3]),0.4*(4.61-Mg[:,3]),0.4*(4.61-H[:,3])  #transfer from Mag to Luminosity
#
#if inp_peng == 1:
#    plt.scatter(C[:21,3],C_cali[:21],c=C[:21,1],s=80,marker="s",zorder=100, vmin=0.3, vmax=5, edgecolors='k', alpha = 0.7)
#    plt.scatter(Mg[:18,3],Mg_cali[:18],c=Mg[:18,1],s=80,marker="s",zorder=100, vmin=0.3, vmax=5, edgecolors='k', alpha = 0.7)
#    plt.scatter(H[:,3],Hb_cali,c=H[:,1],s=80,marker="s",zorder=100, vmin=0.3, vmax=5, edgecolors='k', alpha = 0.7)  #zorder is the position at z foreground.
#    #For non-lensed samples
#    plt.scatter(C[21:,3],C_cali[21:],c=C[21:,1],s=80,marker="s",zorder=100, vmin=0.3, vmax=5, edgecolors='k', alpha = 0.7)
#    plt.scatter(Mg[18:,3],Mg_cali[18:],c=Mg[18:,1],s=80,marker="s",zorder=100, vmin=0.3, vmax=5, edgecolors='k', alpha = 0.7)
##==============================================================================
## Two sample in H0licow VII: i.e. HE0435 and RXJ1131
##==============================================================================
############HE0435 information using Mg#############
#z_0435,L_0435,FWHM_0435=1.69,46-np.log10(5.15),4930
#M_0435=Mgcal(L_0435,FWHM_0435)
#mag_0435=21.75 -1.2514 + dm*pass_dmag(z_0435) #to vega for K correction, -1.2514 is the filter at F160w
#Mag_0435=mag_0435-(20.09-(-23.40))  #23.40 and 20.09 are in Vege system, the K-correction value is the same as AB.
#Mag_0435=Mag_0435+0.21     #change Mg from Vega to AB in R band
#L_0435=0.4*(4.61-Mag_0435)
##print "L0435=",Mag_0435
#plt.scatter(L_0435,M_0435,c=1.69,s=880,marker="*",zorder=100, vmin=0.3, vmax=5, edgecolors='k')
#HE0435= mlines.Line2D([], [], color='m', ls='', marker='*', markersize=28, markeredgecolor='k')
#
############RXJ1131 information using Mg#############
#z_1131,L_1131,FWHM_1131=0.65,45-np.log10(5.15),5630
#M_Mg_1131=Mgcal(L_1131,FWHM_1131)
############RXJ1131 information using Hb#############
#z_1131,L_1131,FWHM_1131=0.65,45-np.log10(3.81),4545.
#M_Hb_1131=Hcal(L_1131,FWHM_1131)
#from RXJ_Mr import Mr_bulge, Mr_disk ###1=bulge; 2=disk, K_correction inferred by Takahiro. For the details see RXJ_Mr.py
#Mag_1131_1= Mr_bulge   
#Mag_1131_2= Mr_disk
##print "L_bulge_1131=",0.4*(4.61-Mag_1131_1),"L_disk_1131=",0.4*(4.61-Mag_1131_2) 
#Mag_1131_tot=-2.5*np.log10(10**((Mag_1131_1+dm*pass_dmag(0.65))/-2.5)+10**((Mag_1131_2+dm*pass_dmag(0.65))/-2.5))  # with brightness correction
#Mag_1131_b=Mag_1131_1+dm*pass_dmag(0.65)  #-21.85 for the bulge, -23.183 for the disk
#M_1131=np.log10((10**M_Hb_1131+10**M_Mg_1131)/2)
#if host==0:
#    Mag_1131=Mag_1131_b
#if host==1:
#    Mag_1131=Mag_1131_tot
#L_1131 = 0.4*(4.61-Mag_1131)
#plt.scatter(L_1131,M_1131,c=0.65,s=880,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
#RXJ1131= mlines.Line2D([], [], color='lightpink', ls='', marker='*', markersize=28, markeredgecolor='k')

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
MBs = load_MBH(ID,MB_ID,if_reportHb=0   )
plt.scatter(lumi_s,MBs,c='tomato',s=580,marker="*",zorder=100, edgecolors='k')
LR_err = load_err(prop = 'LR', ID=ID)
#for i in range(len(lumi_s)):
plt.errorbar(lumi_s,MBs, xerr=[np.abs(LR_err)[:,0], np.abs(LR_err)[:,1]], yerr=0.4, color='blue',ecolor='orange', fmt='.',zorder=-500,markersize=1)
#==============================================================================
# The colorbar label setting up
#==============================================================================
#cl=plt.colorbar()          #cl take the inforamtion from the lastest plt
#cl.set_label('Source redshift',rotation=270,size=20)
#cl.ax.get_yaxis().labelpad=35     #the distance of the colorbar titel from bar
#cl.ax.tick_params(labelsize=30)   #the labe size

if dm ==0:
    plt.title(r"M$_{\rm BH}-$L$_{\rm host}$ relation",fontsize=35)
#elif dm ==1:
#    plt.title("Evolution-corrected $M_{BH}-L_{host}$ relation",fontsize=35)

plt.xlabel(r"log(L$_{\rm R,host}$/L$_{\odot})$",fontsize=35)
plt.ylabel(r"log(M$_{\rm BH}$/M$_{\odot})$",fontsize=35)
plt.xticks(np.arange(5,16.5,0.5))
plt.yticks(np.arange(6,16.5,0.5))
plt.axis([9.0,12.,6.25,10.25])
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)

#lens = mlines.Line2D([], [], color='cyan', ls='', marker='o', markersize=9)
#peng = mlines.Line2D([], [], color='cyan', ls='', marker='s', markersize=9)
park = mlines.Line2D([], [], color='darkseagreen', ls='', marker='^', markersize=9)
new_sample = mlines.Line2D([], [], color='tomato', ls='', marker='*', markersize=20,markeredgecolor='k')

plt.legend([Pkc, park, new_sample],[
"Local AGNs by Bennert+10",
"Intermediate redshift AGNs",
"This work"\
],scatterpoints=1,numpoints=1,loc=2,prop={'size':24},ncol=1)
#else:
#    plt.legend([HE0435,RXJ1131,park,Pkc, new_sample],['HE0435',\
#    'RXJ1131',\
#    "AGNs from P15",
#    #'M$_{BH}$ calibrated from Mg$_{II}$',\
#    #"M$_{BH}$ calibrated from H${\\beta}$",'M$_{BH}$ calibrated from C$_{IV}$'
#    "local AGNs",
#    "our new samples"\
#    ],scatterpoints=1,numpoints=1,loc=2,prop={'size':24},ncol=2)
if dm ==0:
    plt.savefig("MBH-L_obs.pdf")
#elif dm ==1:
#     plt.savefig("MBH-L_ev.pdf")
plt.show()
