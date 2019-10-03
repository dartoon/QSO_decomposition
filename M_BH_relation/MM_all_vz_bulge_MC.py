#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 15:50:27 2018

@author: Dartoon

Note: Using MC to infer the gamma value.
"""
import numpy as np
np.set_printoptions(precision=4)
import matplotlib.pyplot as plt
import matplotlib as mat
import matplotlib.lines as mlines
from matplotlib import colors
import pickle


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
#==============================================================================
#input Park's data 
#==============================================================================
################ bulge or total relaiton? #################
#host= input('with bulge or total relation??? 0 = bulge;   1= total:')
####### input Park data ####
########in AB system, V band#######
#f0 ='data/SS13_MM.txt'
#ss = np.loadtxt(f0)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
#
#f1 ='data/B11_MM.txt'
#b11 = np.loadtxt(f1)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
#
#f2 = 'data/Cisternas_data.txt'
#cis11 = np.loadtxt(f2)  #0 redshift; 1 BH mass; 2 M*;

def lnlike(theta, x, y, yerr):
    b, sint= theta
    model = b*x
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf
def lnprior(theta):
    b, sint	 = theta
    if -10 < b < 10.0 and 0 < sint < 10:
        return 0.0
    return -np.inf
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)
import scipy.optimize as op
import emcee
def find_n(array,value):           #get the corresponding b for a given m 
    idx= (np.abs(array-value)).argmin()
    return array[idx]

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

##exclude 6 outliers: X763, CID543, LID360, XID2396, CDFS-229 and CDFS-724.
#ID = ['CDFS-1', 'CID70',  'SXDS-X735', 'CDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X969','XID2138','XID2202',\
# 'CID206', 'ECDFS-358', 'CID597','CID1281','CID255']
#MB_ID = ['CDFS-1', 'CID70',  'SXDS-X735', 'ECDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X969','LID1820','LID1622',\
#'CID206', 'ECDFS-358', 'CID597','CID1281','CID255']

zs = np.asarray(load_zs(ID))
host_n_inf = np.array(load_n(ID, folder = '../'))
#Mstar = load_host_p(ID)[1]
n_BT_data = pickle.load(open('../Comparsion/CANDELS_catalog/bulge_disk_fit/Sersic_BT_data.pkl','rb'))
MBs = load_MBH(ID,MB_ID,if_reportHb=0)
#    #Plot the median circle.
#    plt.scatter(np.log10(1+zs[?]),(MBs[?]-(m_ml*Mstar[?]+b_ml)),facecolors='none',
#                s=180,marker="o",zorder=900, vmin=0.3, vmax=5, edgecolors='green', linewidth='4',)    
#    plt.scatter(np.log10(1+zs[?]),(MBs[?]-(m_ml*Mstar[?]+b_ml))-0.21,facecolors='none',
#                s=180,marker="o",zorder=900, vmin=0.3, vmax=5, edgecolors='green', linewidth='4',)    
#####fit the evolution##########
################################
import time
gamma = [] 
Mstar_samp = []
Mstar_err = load_err(prop = 'Mstar', ID=ID)
t1 = time.time()
for seed in range(5000):
    host_n = host_n_inf[:,0] + np.random.normal(0, host_n_inf[:,1]) # Randomly sampling the host_n
    host_n[host_n<0.3]=0.3
    host_n[host_n>7]=7    
    host_BT = np.zeros(len(host_n))
    for i in range(len(host_n)):
        n = host_n[i]
        if n <4:
            step = 0.1
        else:
            step = 0.5
        idx = [[n_BT_data[0]>n-step][0] * [n_BT_data[0]<n+step][0]][0]
        BT_list = n_BT_data[1][idx]
        host_BT[i] = np.random.choice(BT_list)
    host_BT[host_BT<0.1] = 0.1  #Set the lower limit?
    Mstar = np.log10((10.**load_host_p(ID)[1])*host_BT)  #Bulge Mstar
    yerr_highz = [(Mstar_err[:,0]**2+0.4**2)**0.5, (Mstar_err[:,1]**2+0.4**2)**0.5]
    z_cosmos, y_cosmos = zs, MBs-(m_ml*Mstar+b_ml)
    yerr_hz = (yerr_highz[0]+ yerr_highz[1])/2
    #    #if consider 32 AGN only:
    z=z_cosmos
    y=y_cosmos
    yerr = yerr_hz    
    yerr = np.sqrt(yerr**2 + sint_ml**2)
    #### fit with emcee ###############
    x=np.log10(1+z)
    y=y
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [1.8, 0.3], args=(x, y, yerr))
    b_ml_offset,_= result["x"]
    xp = np.array([5, 13])
    ndim, nwalkers = 2, 100
    pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
    sampler.run_mcmc(pos, 500)
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    value=round(b_ml_offset,2)
    #####################
    value,sig=b_ml_offset,(np.percentile(samples,84,axis=0)[0]-np.percentile(samples,16,axis=0)[0])/2.
    gamma.append(value + np.random.normal(0, sig))
    Mstar_samp.append(Mstar+np.random.normal(0, ((Mstar_err[:,0]**2+Mstar_err[:,1]**2)/2)**0.5 ))
    if seed/10 > (seed-1)/10:
        t2 = time.time()
        print "seed:", seed, "time spend:", round((t2-t1)/60), "mins", "total time should be:", (t2-t1)/60/(seed+1)*5000, 'mins'
#print t2-t1
#    value,sig = 1.71, 0.27
#    plt.text(0.15, -1.75, "$\Delta$log$M_{BH}$=$(%s\pm%s)$log$(1+z)$"%(value,sig),color='blue',fontsize=25)
b_ml_offset, _ =np.percentile(samples, 50,axis=0)
#print "lnlike=",lnlike(theta=[b_ml_offset, sint_mid],x=x, y=y, yerr=yerr)
xl = np.linspace(0, 5, 100)
b=np.percentile(samples,50,axis=0)[0]
#print samples[:,1][samples[:,0]==find_n(samples[:,0],m)]
for i in range(100):
    posi=np.random.uniform(16,84)
    b=np.percentile(samples,posi,axis=0)[0]    
    plt.plot(xl, xl*0+xl*b, color="lightgray", alpha=0.2,linewidth=7.0,zorder=-1)
plt.plot(xl, xl*0+xl*b_ml_offset, color="red", linewidth=4.0,zorder=0)
plt.plot(xl, xl*0, color="black", linewidth=2.0,zorder=0)
plt.scatter(np.log10(1+zs),MBs-(m_ml*Mstar+b_ml),c='tomato',
            s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
plt.errorbar(np.log10(1+zs),MBs-(m_ml*Mstar+b_ml),
             yerr= yerr_highz,
             color='tomato',ecolor='orange', fmt='.',markersize=1)    
plt.xlabel("log$(1+z)$",fontsize=35)
new_sample = mlines.Line2D([], [], color='tomato', ls='', marker='*', markersize=20,markeredgecolor='k')
plt.xticks(np.arange(-0.1,1,0.1))
xl=-0.01
xh=np.log10(1+2.5)
plt.yticks(np.arange(-5.5,6,0.5))
plt.axis([xl,xh,-2.0,3.5])
plt.ylim([-2.0,3.5])
plt.ylabel("$\Delta$log$M_{BH}$ (vs $ M_{*,bulge}$)",fontsize=35)
plt.grid()
plt.tick_params(labelsize=25)
ax2=ax.twiny()
tticks=np.array([10**xl-1,0.5,1,1.5,2,10**xh-1])
ax2.set_xticks([np.log(t+1) for t in tticks])  # for the entire scale
ax2.set_xticklabels([0,0.5,1,1.5,2,2.5])  # 0 actuall is corresponds to 10**-0.01-1
ax2.set_xlabel('$z$',fontsize=35)
plt.tick_params(labelsize=25)
SS13 = mlines.Line2D([], [], color='darkseagreen', ls='', marker='^', markersize=8)
plt.legend([Bkc, Hkc, new_sample],[
'Local by Bennert+11',\
"Local by H&R",
"High-z bulge"
],scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
#plt.savefig("MBH-Mbulge-style{0}.pdf".format(style))
plt.show()

filename = "gamma_Mbulge_result_(BT_min_0.1)_201910.pkl"
pickle.dump([gamma,Mstar_samp], open(filename, 'wb'))

#%%
gamma,Mstar_samp = pickle.load(open(filename,'rb'))
Mstar_samp = np.asarray(Mstar_samp)
hist_value = gamma#Mstar_samp[:,1]
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
#plt.savefig("gamma_hist.pdf")
plt.show()
#m, l, h = 1.8614898041641283, 1.544294753221521, 2.1670754259051903

##%%Print the Mstar:
#tab_list = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 'CID3570',
#            'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 'CID70', 'LID1273',
#            'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396', 'CDFS-1', 'CDFS-229',
#            'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']
#for i in range(len(tab_list)):
#    j = [jj for jj in range(len(ID)) if tab_list[i] == ID[jj]][0]
##    print tab_list[i], ID[j]
#    print tab_list[i], round(Mstar[j],2)
#
