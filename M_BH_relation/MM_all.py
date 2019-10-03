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
########## input local data ####
#==============================================================================
# The seleting for dm and host_total and dmag are in this local
#==============================================================================
from local_MMstar import *
#==============================================================================
#input SS13 and B11 data 
#==============================================================================
################ bulge or total relaiton? #################
####### input SS13 data ####
#######in AB system, V band#######
inp_SS13 = 1
f0 ='data/SS13_MM.txt'
ss = np.loadtxt(f0)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
if inp_SS13 ==1:
    plt.scatter(ss[:,1],ss[:,2],c='darkseagreen',marker="^",s=180,zorder=100, edgecolors='white')
  
inp_b11= 1
f1 ='data/B11_MM.txt'
b11 = np.loadtxt(f1)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
if inp_b11 ==1:
    plt.scatter(b11[:,1],b11[:,2],c='darkseagreen',marker="^",s=180,zorder=100, edgecolors='white')

inp_Cis= 1
f2 = 'data/Cisternas_data.txt'
cis11 = np.loadtxt(f2)  #0 redshift;
if inp_Cis ==1:
    plt.scatter(cis11[:,2],cis11[:,1],c='darkseagreen',marker="^",s=180,zorder=100, edgecolors='white')

inp_Knud= 0
f3 = 'data/high_edd_agn.txt'
Knud = np.loadtxt(f3)[:,2:]  # 0 redshift; 1 L_bol; 2 M_BH; 3 M_acc; 4 M_*
if inp_Knud ==1:
    plt.scatter(Knud[:,4], Knud[:,2],c='blue',marker="o",s=180,zorder=101, edgecolors='white')


tx, ty = 11.85, 7.
plt.text(tx, ty, "intermediate\n  sample\nuncertainties",  fontsize=20)
plt.errorbar(tx+0.2,ty-0.08, xerr=0.2, yerr=0.4, color='darkseagreen',ecolor='black', fmt='^',zorder=-500,markersize=10)
#%%
#==============================================================================
# My new inference
#==============================================================================
from scipy.integrate import quad
from load_result import load_host_p, load_MBH, load_err
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)

h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]

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
MBs = load_MBH(ID,MB_ID)
plt.scatter(Mstar,MBs,c='tomato',s=580,marker="*",zorder=100, edgecolors='k')
Mstar_err = load_err(prop = 'Mstar', ID=ID)
#for i in range(len(lumi_s)):
plt.errorbar(Mstar,MBs, xerr=[np.abs(Mstar_err)[:,0], np.abs(Mstar_err)[:,1]], yerr=0.4, color='blue',ecolor='orange', fmt='.',zorder=-500,markersize=1)

#==============================================================================
# The colorbar label setting up
#==============================================================================
#cl=plt.colorbar()          #cl take the inforamtion from the lastest plt
#cl.set_label('Source redshift',rotation=270,size=20)
#cl.ax.get_yaxis().labelpad=35     #the distance of the colorbar titel from bar
#cl.ax.tick_params(labelsize=30)   #the labe size

plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
#plt.xticks(np.arange(5,16.5,0.5))
#plt.yticks(np.arange(6,16.5,0.5))
plt.xlim(9,12.5)
plt.ylim(6.0,10)
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)

Bkc=mlines.Line2D([], [], color='gray', ls='', marker='.', markersize=15)
Hkc=mlines.Line2D([], [], color='black', ls='', marker='.', markersize=15)
SS13 = mlines.Line2D([], [], color='darkseagreen', ls='', marker='^', markersize=13)
new_sample = mlines.Line2D([], [], color='tomato', ls='', marker='*', markersize=16,markeredgecolor='k')
#unlens = mlines.Line2D([], [], color='cyan', ls='', marker='s', markersize=9)

plt.legend([Bkc,Hkc,SS13,new_sample],[
'Local by Bennert+11',\
"Local by H&R",
"intermediate redshift AGNs",
"This work"\
],scatterpoints=1,numpoints=1,loc=2,prop={'size':20},ncol=2)
#plt.savefig("MBH-Mstar.pdf")
plt.show()
