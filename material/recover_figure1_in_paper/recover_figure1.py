#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 10:30:19 2019

@author: Dartoon

Recover the plot 

The four outliers are ['CDFS-1', 'SXDS-X1136', 'SXDS-X763', 'CDFS-724']
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
#from creat_table_for_manscript import load_result_here

import sys
sys.path.insert(0,'../../py_tools')
from load_result import load_MBH, load_Lbol
import copy
tab_list = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 'CID3570',
            'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 'CID70', 'LID1273',
            'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396', 'CDFS-1', 'CDFS-229',
            'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']
ext_ID = {'XID2202':'LID1622', 'XID2138':'LID1820', 'XID2396':'LID1878', 'CDFS-321':'ECDFS-321'}
MB_Lbol_info = load_Lbol(tab_list, folder = '../../') 
outliers = ['CDFS-1', 'SXDS-X1136', 'SXDS-X763', 'CDFS-724']
#outliers = ['CDFS-1', 'SXDS-X763', 'CDFS-724']
tab_sub_list = copy.deepcopy(tab_list)
for i in range(len(tab_sub_list)):
    if tab_sub_list[i] in ext_ID.keys():
        tab_sub_list[i] = ext_ID[tab_sub_list[i]]
MBH = load_MBH(tab_list,tab_sub_list, if_reportHb=0, folder = '../../')   
print  "Calculate the Eddington ratio:"
logLedd = 38 + np.log10(1.2) + MBH
logEdd = MB_Lbol_info - logLedd
#for i in range(len(tab_list)):
#    print tab_list[i], round(MBH[i],3), round(logEdd[i], 3)
    
#%%Define the rho functions:
#calculate BHMF
def rho_bh(mbh, alpha=-1.50, beta=0.96, mbh_star= 9.09):
    c = 1
    rho = c * (10**mbh/10**mbh_star)**(alpha+1)*np.exp(-(10**mbh/10**mbh_star)**beta)
    return np.log10(rho)

#calculate the rho_lam in 2D
def rho_lam_2d(lam, mbh, alpha=-0.29 , lam_star=-1.19, k_lam=0.099, logMc=8.):
    '''
    all the lam and mbh in log 
    '''
    lamstar = lam_star + k_lam*(mbh-logMc)   #lam_star in table should be in log sample already.
    rho_lam_vs_Mbh = 1/np.log10(np.e) * (10**lam/10**lamstar)**(alpha+1) * np.exp(-(10**lam/10**lamstar))
    return np.log10(rho_lam_vs_Mbh)
vec_rho_lam_2d=np.vectorize(rho_lam_2d)
mbh = np.linspace(7.4,10, 100)
lam = np.linspace(-2, 0.25, 101)
rho_lam_vs_mbh_list = []
for i in range(len(mbh)):
    rho_lam_vs_mbh_list.append(vec_rho_lam_2d(lam, mbh[i]))
rho_lam_vs_mbh = np.asarray(rho_lam_vs_mbh_list)
import matplotlib
my_cmap = copy.copy(matplotlib.cm.get_cmap('Oranges',10)) # copy the default cmap
rho_2d = np.asarray([(rho_lam_vs_mbh[i]+rho_bh(mbh[i])) for i in range(len(mbh))])
for i in range(len(mbh)):
#    print mbh[i], rho_bh(mbh[i])
    plt.scatter(lam*0+mbh[i], lam, c= rho_2d[i], s = 140, vmin = np.min(rho_2d), vmax = np.max(rho_2d),
                   marker='s', alpha=0.9, edgecolors='none', cmap = my_cmap)
cl=plt.colorbar()
cl.set_label('Value in Col. 4', size=20)
plt.show()

#%% To plot the main figure
fig = plt.figure(1, figsize=(8,8))
from mpl_toolkits.axes_grid1 import make_axes_locatable
marker =  ["o" , "v" , "s"]
colors =  ['r','g','b']
lables = ['ID', 'SXDS', 'CDFS']
# the scatter plot:
ax = plt.subplot(111)
for i in range(len(tab_list)):
    k = np.where([lables[j] in tab_list[i] for j in range(3)])[0][0]
    ax.scatter(MBH[i], logEdd[i],marker=marker[k], color=colors[k]) 
#    if tab_list[i] in outliers:
#        ax.scatter(MBH[i], logEdd[i],marker=marker[k], color='k', zorder=1) 
x_cline = np.linspace(7.1, 10)
y_line_0 = x_cline*0 - 1.5
ax.plot(x_cline, y_line_0,'k--')
y_line_1 = x_cline*0 + 0.0
ax.plot(x_cline, y_line_1,'k--')
y_cline = np.linspace(-20, 1, 10)
x_line = y_cline*0 + 8.6
ax.plot(x_line, y_cline,'r--')
x_cline = np.linspace(7.1, 8.5)
y_line_3 = -1.1*(x_cline-7.5) -0.5
#ax.plot(x_cline, y_line_3,'k--')
plt.xlabel("log$(M_{BH}/M_{\odot})$",fontsize=25)
plt.ylabel("log$L_{bol}/L_{Edd}$", fontsize=25) 
plt.tick_params(labelsize=15)
plt.xticks(np.arange(7.5,9.6,0.5))
other_data = np.loadtxt('Peng_Decarli_data.txt')
Peng_data = other_data[other_data[:,0]==0][:,1:]
Peng_data[:,0] = np.log10(Peng_data[:,0]*10**9.)
Peng_data[:,1] = np.log10(Peng_data[:,1])
Peng_data = Peng_data[(Peng_data[:,2]>1.)*(Peng_data[:,2]<20)]
Decarli_data = other_data[other_data[:,0]==1][:,1:]
Decarli_data = Decarli_data[(Decarli_data[:,2]>1.)*(Decarli_data[:,2]<20)]
ax.scatter(Peng_data[:,0], Peng_data[:,1], marker='s', edgecolors='gray', facecolors= 'none')
ax.scatter(Decarli_data[:,0], Decarli_data[:,1],marker='o', edgecolors='gray', facecolors = 'none')
import matplotlib.lines as mlines
cosmos=mlines.Line2D([], [], color=colors[0], ls='', marker=marker[0], markersize=7)
sxds=mlines.Line2D([], [], color=colors[1], ls='', marker=marker[1], markersize=7)
cdfs = mlines.Line2D([], [], color=colors[2], ls='', marker=marker[2], markersize=7)
peng = mlines.Line2D([], [],  color='none',markeredgecolor='gray', ls='', marker='s', markersize=7)
decarli = mlines.Line2D([], [], color='none', markeredgecolor='gray', ls='', marker='o', markersize=7)
plt.legend([cosmos,peng,sxds, decarli, cdfs],["COSMOS", "Peng 2006a", "SXDS", "Decarli 2010", "CDFS"],
           scatterpoints=1,numpoints=1,loc=3,prop={'size':13},ncol=3)

#%% create new axes on the right and on the top of the current axes
# The first argument of the new_vertical(new_horizontal) method is
# the height (width) of the axes to be created in inches.
divider = make_axes_locatable(ax)
#For plot the BHMF:
mbh_star = 9.09
alpha = -1.50
beta = 0.96
mbh_x = mbh
#mbh_rho = rho_bh(mbh_x)
mbh_rho = np.log10(np.sum(10**(rho_2d), axis=1)) # !!! Note the way to combining the likelihood.
#normize mbh_rho:
mbh_rho = 10**mbh_rho
mbh_rho /= np.sum(mbh_rho)* (mbh_x[-1]-mbh_x[1])/len(mbh_x)
mbh_rho = np.log10(mbh_rho)
ax_x = divider.append_axes("top", 1.5, pad=0.1, sharex=ax)
ax_x.plot(mbh_x,mbh_rho, 'k')

#lam_reg = [-1.5,0]
#ind0, ind1= np.where(lam>=-1.5)[0][0], np.where(lam<=0.)[0][-1]
#mbh_rho_select = np.log10(np.sum(10**(rho_2d[:,ind0:ind1]), axis=1)) # !!! Note the way to combining the likelihood.
##normize mbh_rho_select:
#mbh_rho_select = 10**mbh_rho_select
#mbh_rho_select /= np.sum(mbh_rho_select)* (mbh_x[-1]-mbh_x[1])/len(mbh_x)
#mbh_rho_select = np.log10(mbh_rho_select)
#ax_x.plot(mbh_x,mbh_rho_select, 'b')

ax_x.plot(x_line, y_cline,'r--')
plt.yticks(np.arange(-10,10,1))
plt.ylim([int(mbh_rho.min()),int(mbh_rho.max())+1])
plt.tick_params(labelsize=15)
plt.ylabel("log$\phi(M_{BH})$", fontsize=25) 


#%%    
ax_y = divider.append_axes("right", 1.2, pad=0.1, sharey=ax)
#ax_y.hist(logEdd, orientation='horizontal')
rho_lam = np.log10(np.sum(10**(rho_2d), axis=0)) # !!! Note the way to combining the likelihood.
#normize rho_lam:
rho_lam = 10**rho_lam
rho_lam /= np.sum(rho_lam)* (lam[-1]-lam[1])/len(lam)
rho_lam = np.log10(rho_lam)
ax_y.plot(rho_lam, lam, 'k')#, orientation='horizontal')
x_cline = np.linspace(-10, 10)
y_line_0 = x_cline*0 - 1.5
ax_y.plot(x_cline, y_line_0, 'k--')
y_line_1 = x_cline*0 + 0.0
ax_y.plot(x_cline, y_line_1, 'k--')
plt.xticks(np.arange(-10,10,2))
plt.xlim([int(rho_lam.min()),int(rho_lam.max())+1])
plt.tick_params(labelsize=15)
plt.xlabel("log$\phi(\lambda))$", fontsize=25) 

#%%
# the xaxis of ax_x and yaxis of ax_y are shared with ax,
# thus there is no need to manually adjust the xlim and ylim of these
# axis.

#ax_x.axis["bottom"].major_ticklabels.set_visible(False)
for tl in ax_x.get_xticklabels():
    tl.set_visible(False)
for tl in ax_y.get_yticklabels():
    tl.set_visible(False) 
ax.set_xlim([7.4, 9.7])
ax.set_ylim([-1.9, 0.3])
#plt.savefig("AGN_selection.pdf")
plt.show()

#%%Print them for Andreas
#for i in range(len(tab_list)):
#    print round(MBH[i],3), round(MB_Lbol_info[i],3), round(logEdd[i], 3)


