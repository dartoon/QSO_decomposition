#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 23:06:35 2019

@author: Dartoon

Read the inference from gamma
"""
import numpy as np
np.set_printoptions(precision=4)
import matplotlib.pyplot as plt
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
import matplotlib.lines as mlines
from matplotlib import colors
import pickle

gamma,Mstar_samp = pickle.load(open('gamma_Mbulge_result_(BT_min_0.1)_201910.pkl','rb'))
Mstar_samp = np.asarray(Mstar_samp)
hist_value = gamma #Mstar_samp[:,1]
plt.figure(figsize=(10,8))
l =np.percentile(hist_value,16,axis=0)
m =np.percentile(hist_value,50,axis=0)
h =np.percentile(hist_value,84,axis=0)
yl = np.linspace(0,4,100)
x_l = yl*0+ l
x_m = yl*0 + m
x_h = yl*0 + h
plt.plot(x_l, yl, 'k--',linewidth=4.0)
plt.plot(x_m, yl, 'blue',linewidth=4.0)
plt.plot(x_h, yl, 'k--',linewidth=4.0)
plt.axis([m-4*(m-l),m+(h-m)*4,0,1.5])
plt.hist(hist_value,bins=30, normed=True)
plt.xlabel('$\gamma$',fontsize=35)
plt.ylabel('PDF',fontsize=35)
plt.tick_params(labelsize=25)
plt.yticks([])
plt.savefig("gamma_hist.pdf")
plt.show()

#%%
ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']

tab_list = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 'CID3570',
            'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 'CID70', 'LID1273',
            'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396', 'CDFS-1', 'CDFS-229',
            'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']

Mstar_err = []
Mstar = []
#plot the estimated bulge relation:
for i in range(len(ID)):
    hist_value = Mstar_samp[:,i]
    l =np.percentile(hist_value,16,axis=0)
    m =np.percentile(hist_value,50,axis=0)
    h =np.percentile(hist_value,84,axis=0)
    Mstar.append(m)
    Mstar_err.append([m-l, h-m])
Mstar = np.asarray(Mstar)
Mstar_err = np.asarray(Mstar_err)

for target in tab_list:
    i = [i for i in range(len(ID)) if target == ID[i]][0]
    hist_value = Mstar_samp[:,i]
    l =np.percentile(hist_value,16,axis=0)
    m =np.percentile(hist_value,50,axis=0)
    h =np.percentile(hist_value,84,axis=0)
#    print target, "{0}\substack+{1}-{2}".format(round(m,2), round(h-m,2), round(m-l,2))
 
#%% Plot the bulge offset for BH relation
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
#==============================================================================
#input Park's data 
#==============================================================================
################ bulge or total relaiton? #################
#host= input('with bulge or total relation??? 0 = bulge;   1= total:')
####### input Park data ####
########in AB system, V band#######
f0 ='data/SS13_MM.txt'
ss = np.loadtxt(f0)[:,1:]  #0 redshift; 1 M*; 2 BH mass; 3 B/T


f1 ='data/B11_MM.txt'
b11 = np.loadtxt(f1)[:,1:]  #0 redshift; 1 M*; 2 BH mass; 3 M*_bulge
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
MBs = load_MBH(ID,MB_ID,if_reportHb=0)

if style ==0:
    plt.scatter(np.log10(1+ss[:,0]), 10** (ss[:,2]-np.log10(10**(ss[:,1])*ss[:,3])), c='darkseagreen',
                s=180,marker="^", zorder=100,vmin=0.3, vmax=2, edgecolors='white')
    
    plt.scatter(np.log10(1+b11[:,0]), 10** (b11[:,2]-b11[:,3]), c='darkseagreen',
                s=180,marker="^", zorder=100,vmin=0.3, vmax=2, edgecolors='white')    
    
    plt.scatter(np.log10(1+zs), 10**(MBs-Mstar),c='tomato',
                s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')

    plt.scatter(np.log10(1+zs[9]), np.median(np.concatenate([10**(hloc[:,3]-hloc[:,1]),10**(bloc[:,3]-bloc[:,1])])),
                facecolors='none', s=280,marker="o",zorder=900, vmin=0.3, vmax=5, edgecolors='blue', linewidth='6', alpha=0.5)       
    plt.scatter(np.log10(1+zs[9]), np.median(np.concatenate([10**(hloc[:,3]-hloc[:,1]+0.21),10**(bloc[:,3]-bloc[:,1]+0.21)])),
                facecolors='none', s=280,marker="o",zorder=900, vmin=0.3, vmax=5, edgecolors='blue', linewidth='6', alpha=0.5)      
    plt.arrow(np.log10(1+zs[9]),np.median(np.concatenate([10**(hloc[:,3]-hloc[:,1]),10**(bloc[:,3]-bloc[:,1])])), 0, +0.0009,
              zorder=800, head_length= 0.001374/8,head_width= 0.005,fc='k',ec='k')
    
#    
if style ==1:  
    plt.errorbar(np.log10(1+ss[:,0]),ss[:,2]-(m_ml* np.log10(10**(ss[:,1])*ss[:,3]) +b_ml),yerr=(0.4**2+0.2**2)**0.5,fmt='^',color='darkseagreen',markersize=9)
    plt.errorbar(np.log10(1+b11[:,0]),b11[:,2]-(m_ml*b11[:,3]+b_ml),yerr=(0.4**2+0.2**2)**0.5,fmt='^',color='darkseagreen',markersize=9)  
    samples = gamma

    yerr_highz = [(Mstar_err[:,0]**2+0.4**2)**0.5, (Mstar_err[:,1]**2+0.4**2)**0.5]
    for i in range(100):
        posi=np.random.uniform(16,84)
        b=np.percentile(samples,posi,axis=0)#[0]
        plt.plot(xl, xl*0+xl*b, color="lightgray", alpha=0.2,linewidth=7.0,zorder=-1)
    plt.plot(xl, xl*0+xl*np.percentile(gamma,50,axis=0), color="red", linewidth=4.0,zorder=0)        
    plt.plot(xl, xl*0, color="black", linewidth=2.0,zorder=0)
    plt.scatter(np.log10(1+zs),MBs-(m_ml*Mstar+b_ml),c='tomato',
                s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
    plt.errorbar(np.log10(1+zs),MBs-(m_ml*Mstar+b_ml),
                 yerr= yerr_highz,
                 color='tomato',ecolor='orange', fmt='.',markersize=1) 
    
    
    plt.scatter(np.log10(1+zs[9]),0,facecolors='none',
                s=280,marker="o",zorder=900, edgecolors='blue', linewidth='6', alpha=0.5)    
    plt.scatter(np.log10(1+zs[9]),0+0.21,facecolors='none',
                s=280,marker="o",zorder=900, edgecolors='blue', linewidth='6', alpha=0.5)    
    plt.arrow(np.log10(1+zs[9]),0, 0, +0.19, zorder=800, head_length=0.05,head_width=0.005,fc='k',ec='k')
       
#plt.text(0.15, -1.75, "$\Delta$log$M_{BH}$=$(%s\pm%s)$log$(1+z)$"%(1.86,0.30),color='blue',fontsize=25)
    value,sig=np.percentile(gamma,50,axis=0),round((np.percentile(gamma,84,axis=0)-np.percentile(gamma,16,axis=0))/2,2)
#    print value, sig

plt.xlabel("log(1+z)",fontsize=45)
new_sample = mlines.Line2D([], [], color='tomato', ls='', marker='*', markersize=20,markeredgecolor='k')
plt.xticks(np.arange(-0.1,1,0.1))
xl=-0.01
xh=np.log10(1+2.5)
if style ==0:
    ax.set_yscale('log')
    plt.axis([xl,xh,0,0.7])
    plt.ylabel(r"M$_{\rm BH}$/M$_{*, {\rm bulge}}$",fontsize=45)
if style ==1:
    plt.yticks(np.arange(-5.5,6,0.5))
    plt.axis([xl,xh,-2.0,3.5])
    plt.ylim([-2.0,3.5])
    plt.ylabel(r"$\Delta$logM$_{\rm BH}$ (vs  M$_{*, {\rm bulge}}$)",fontsize=45)
plt.grid()
plt.tick_params(labelsize=35)
ax2=ax.twiny()
tticks=np.array([10**xl-1,0.5,1,1.5,2,10**xh-1])
ax2.set_xticks([np.log(t+1) for t in tticks])  # for the entire scale
ax2.set_xticklabels([0,0.5,1,1.5,2,2.5])  # 0 actuall is corresponds to 10**-0.01-1
ax2.set_xlabel('z',fontsize=45)
plt.tick_params(labelsize=35)
SS13 = mlines.Line2D([], [], color='darkseagreen', ls='', marker='^', markersize=8)
plt.legend([Bkc, Hkc, SS13, new_sample],[
'Local by Bennert+11',\
"Local by H&R",
"Intermediate redshift AGNs",
"High-z bulge"
],scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.savefig("MBH-Mbulge-style{0}.pdf".format(style))
plt.show()

 #%%
##calcualte the mean offset, need to run something?#!!!
#weighted_offset = np.sum(np.asarray(MBs-(m_ml*Mstar+b_ml))*yerr_highz) / np.sum(yerr_highz)                              
#rms_offset = np.sqrt(np.sum((np.asarray(MBs-(m_ml*Mstar+b_ml))-weighted_offset)**2*yerr_highz) / np.sum(yerr_highz))
##rms_offset = np.sqrt(np.mean((np.asarray(y_cosmos)-weighted_offset)**2)) #No error bar
#print weighted_offset, rms_offset/np.sqrt(32.)


#%% Plot Figure12-b:
#The histogram of the M_BH/M_star.
high_ratio = []
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111)
for i in range(len(Mstar_samp)):
    high_ratio.append(10**(MBs-Mstar_samp[i,:])) 
dis0 = np.asarray(high_ratio)
dis0 = dis0.flatten()
high0, x0, _ = plt.hist(np.log10(dis0),normed=True , histtype=u'step',
         label=('High-z sample'), linewidth = 4, color='red')

dis1 = np.concatenate([10**(hloc[:,3]-hloc[:,1]), 10**(bloc[:,3]-bloc[:,1])])

high1, x1, _ = plt.hist(np.log10(dis1),normed=True, histtype=u'step',
         label=('Local sample'), linewidth = 2, color='gray',linestyle=('dashed'))

x0_m = np.median(np.log10(dis0))
high_m0 = high0[np.where(abs(x0_m-x0) == abs(x0_m-x0).min())[0][0]-1]
#plt.plot(np.linspace(0,high_m0)*0+np.median(x0_m) , np.linspace(0,high_m0), linewidth = 4,color='orange',linestyle=('dashed'))

x1_m = np.median(np.log10(dis1))
high_m1 = high1[np.where(abs(x1_m-x1) == abs(x1_m-x1).min())[0][0]-1]
#plt.plot(np.linspace(0,high_m1)*0+np.median(x1_m) , np.linspace(0,high_m1), linewidth = 4, color='green')
#
#plt.text(np.median(x0_m)-0.2, high_m0*1.05, '{0}'.format(round(np.median(x0_m),1)), color='orange',fontsize=25)
#plt.text(np.median(x1_m)-0.2, high_m1*1.05, '{0}'.format(round(np.median(x1_m),1)), color='green',fontsize=25)

plt.xlabel(r"M$_{\rm BH}$/M$_{*, {\rm bulge}}$",fontsize=40)
plt.ylabel("PDF",fontsize=40)
plt.tick_params(labelsize=30)
plt.legend(prop={'size':30})

fig.canvas.draw()
labels = [item.get_text().encode('ascii', 'replace').replace('?','-') for item in ax.get_xticklabels()]
for i in range(len(labels)-2):
    labels[i+1] = '10$^{'+ labels[i+1] + '}$'
plt.yticks([])
ax.set_xticklabels(labels)

plt.savefig('realizations_M_BH_M_star.pdf')
plt.show()
    

