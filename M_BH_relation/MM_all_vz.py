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
#==============================================================================
#input Park's data 
#==============================================================================
################ bulge or total relaiton? #################
#host= input('with bulge or total relation??? 0 = bulge;   1= total:')
####### input Park data ####
#######in AB system, V band#######
f0 ='data/SS13_MM.txt'
ss = np.loadtxt(f0)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
#if inp_SS13 ==1:
#    plt.scatter(ss[:,1],ss[:,2],c=ss[:,0],marker="^",s=180,zorder=100,vmin=0.5, vmax=2, edgecolors='white')

f1 ='data/B11_MM.txt'
b11 = np.loadtxt(f1)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
#if inp_b11 ==1:
#    plt.scatter(b11[:,1],b11[:,2],c=b11[:,0],marker="^",s=180,zorder=100,vmin=0.3, vmax=2, edgecolors='white')


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
host_n = np.array(load_n(ID, folder = '../'))[:,0]
Mstar = load_host_p(ID)[1]
MBs = load_MBH(ID,MB_ID,if_reportHb=0)
Mstar_err = load_err(prop = 'Mstar', ID=ID)
yerr_highz = [(Mstar_err[:,0]**2+0.4**2)**0.5, (Mstar_err[:,1]**2+0.4**2)**0.5]

#plt.scatter(Mstar,MBs,c=zs,s=880,marker="*",zorder=100, vmin=0.3, vmax=2, edgecolors='k')

#plt.errorbar(np.log10(1+zs[MBs!=-99]),MBs[MBs!=-99]-(m_ml*lumi_s[MBs!=-99]+b_ml),yerr=(0.4**2+0.2**2)**0.5,fmt='x',color='royalblue',markersize=28,zorder=250)#, mec='k')
if style ==0:
    plt.scatter(np.log10(1+ss[:,0]), 10** (ss[:,2]-ss[:,1]), c='darkseagreen',
                s=180,marker="^", zorder=100,vmin=0.3, vmax=2, edgecolors='white')
    
    plt.scatter(np.log10(1+b11[:,0]), 10** (b11[:,2]-b11[:,1]), c='darkseagreen',
                s=180,marker="^", zorder=100,vmin=0.3, vmax=2, edgecolors='white')    
    
    plt.scatter(np.log10(1+zs[MBs!=-99]), 10**( MBs[MBs!=-99]- Mstar[MBs!=-99]),c='tomato',
                s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
if style ==1:
#    plt.scatter(np.log10(1+ss[:,0]),ss[:,2]-(m_ml*ss[:,1]+b_ml), c='darkseagreen',
#                s=180,marker="^", zorder=100,vmin=0.3, vmax=2, edgecolors='white')
    plt.errorbar(np.log10(1+ss[:,0]),ss[:,2]-(m_ml*ss[:,1]+b_ml),yerr=(0.4**2+0.2**2)**0.5,fmt='^',color='darkseagreen',markersize=9)

#    plt.scatter(np.log10(1+b11[:,0]),b11[:,2]-(m_ml*b11[:,1]+b_ml), c='darkseagreen',
#                s=180,marker="^", zorder=100,vmin=0.3, vmax=2, edgecolors='white')    
    plt.errorbar(np.log10(1+b11[:,0]),b11[:,2]-(m_ml*b11[:,1]+b_ml),yerr=(0.4**2+0.2**2)**0.5,fmt='^',color='darkseagreen',markersize=9)  #used to be tomato
    plt.scatter(np.log10(1+zs[MBs!=-99]),MBs[MBs!=-99]-(m_ml*Mstar[MBs!=-99]+b_ml),c='tomato',
                s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
    plt.errorbar(np.log10(1+zs[MBs!=-99]),MBs[MBs!=-99]-(m_ml*Mstar[MBs!=-99]+b_ml),
                 yerr= yerr_highz,
                 color='tomato',ecolor='orange', fmt='.',markersize=1)    
    
#    #Plot the median circle.
#    plt.scatter(np.log10(1+zs[9]),(MBs[9]-(m_ml*Mstar[9]+b_ml)),facecolors='none',
#                s=180,marker="o",zorder=900, vmin=0.3, vmax=5, edgecolors='green', linewidth='4',)    
#    plt.scatter(np.log10(1+zs[9]),(MBs[9]-(m_ml*Mstar[9]+b_ml))-0.21,facecolors='none',
#                s=180,marker="o",zorder=900, vmin=0.3, vmax=5, edgecolors='green', linewidth='4',)    
#    plt.arrow(np.log10(1+zs[9]),(MBs[9]-(m_ml*Mstar[9]+b_ml)), 0, -0.2, zorder=900, head_length=0.05,head_width=0.01,fc='k',ec='k')
    
    #####fit the evolution##########
    ################################
    z_ss, y_ss = ss[:,0], ss[:,2]-(m_ml*ss[:,1]+b_ml)
    
    z_b11, y_b11 = b11[:,0], b11[:,2]-(m_ml*b11[:,1]+b_ml)
    
    z_cosmos, y_cosmos = zs[MBs!=-99], MBs[MBs!=-99]-(m_ml*Mstar[MBs!=-99]+b_ml)
                                  
    
    z=np.concatenate((z_ss, z_b11, z_cosmos),axis=0)
    y=np.concatenate((y_ss, y_b11, y_cosmos),axis=0)

    yerr_imd= np.zeros(len(z_ss)+len(z_b11))+(0.4**2+0.2**2)**0.5   # the error for the fitting
    yerr_hz = (yerr_highz[0]+ yerr_highz[1])/2
    yerr = np.concatenate((yerr_imd, yerr_hz),axis=0)
    
#    #if consider 32 AGN only:
#    z=z_cosmos
#    y=y_cosmos
#    yerr = yerr_hz    
    
    yerr = np.sqrt(yerr**2 + sint_ml**2)
    
    #### fit with emcee ###############
    x=np.log10(1+z)
    y=y
    
    def lnlike(theta, x, y, yerr):
        b, sint= theta
        model = b*x
        sigma2 = (yerr**2 + sint**2)
        if sint>=0 :
          return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
        else:
          return -np.inf
    
    import scipy.optimize as op
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [1.8, 0.3], args=(x, y, yerr))
    b_ml,_= result["x"]
    #print b_ml, sint_ml, "ka=",lnlike(theta=[b_ml, sint_ml],x=loc[:,0], y=loc[:,1], yerr=loc[:,2])
    
    xp = np.array([5, 13])
    #plt.plot(xp, m_ml*xp+b_ml, 'r-')
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
    ndim, nwalkers = 2, 100
    pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    import emcee
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
    sampler.run_mcmc(pos, 500)
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    
    b_ml, _ =np.percentile(samples, 50,axis=0)
    #print "lnlike=",lnlike(theta=[b_ml, sint_mid],x=x, y=y, yerr=yerr)
    xl = np.linspace(0, 5, 100)
    plt.plot(xl, xl*0+xl*b_ml, color="red", linewidth=4.0,zorder=0)
    
    plt.plot(xl, xl*0, color="black", linewidth=2.0,zorder=0)
    def find_n(array,value):           #get the corresponding b for a given m 
        idx= (np.abs(array-value)).argmin()
        return array[idx]
    b=np.percentile(samples,50,axis=0)[0]
    #print samples[:,1][samples[:,0]==find_n(samples[:,0],m)]
    for i in range(100):
        posi=np.random.uniform(16,84)
        b=np.percentile(samples,posi,axis=0)[0]    
        #print b
        plt.plot(xl, xl*0+xl*b, color="lightgray", alpha=0.2,linewidth=7.0,zorder=-1)
    value=round(b_ml,2)
    #####################
    value,sig=round(b_ml,2),round((np.percentile(samples,84,axis=0)[0]-np.percentile(samples,16,axis=0)[0])/2,2)
    print value,sig
    plt.text(0.15, -1.75, "$\Delta$log$M_{BH}$=$(%s\pm%s)$log$(1+z)$"%(value,sig),color='blue',fontsize=25)
    #plt.text(0.15, -1.75, "$M_{BH} $VS$ L_{host}\propto (1+z)^{%s\pm%s}$"%(value,sig),color='blue',fontsize=25)


plt.xlabel("log$(1+z)$",fontsize=35)
new_sample = mlines.Line2D([], [], color='tomato', ls='', marker='*', markersize=20,markeredgecolor='k')

plt.xticks(np.arange(-0.1,1,0.1))
xl=-0.01
xh=np.log10(1+2.5)
if style ==0:
    ax.set_yscale('log')
    plt.axis([xl,xh,0,0.5])
    plt.ylabel("$M_{BH}/M_*$",fontsize=35)
if style ==1:
    plt.yticks(np.arange(-5.5,6,0.5))
    plt.axis([xl,xh,-2.0,3.5])
    plt.ylim([-2.0,3.5])
    plt.ylabel("$\Delta$log$M_{BH}$ (vs $M_*$)",fontsize=35)
plt.grid()
plt.tick_params(labelsize=25)

ax2=ax.twiny()
tticks=np.array([10**xl-1,0.5,1,1.5,2,10**xh-1])
ax2.set_xticks([np.log(t+1) for t in tticks])  # for the entire scale
ax2.set_xticklabels([0,0.5,1,1.5,2,2.5])  # 0 actuall is corresponds to 10**-0.01-1
ax2.set_xlabel('$z$',fontsize=35)
plt.tick_params(labelsize=25)

SS13 = mlines.Line2D([], [], color='darkseagreen', ls='', marker='^', markersize=8)

plt.legend([Bkc, Hkc, SS13, new_sample],[
'Local by Bennert+11',\
"Local by H&R",
"Intermediate redshift AGNs",
"This work"
],scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
#plt.savefig("MBH-Mstar-vz_style{0}.pdf".format(style))
plt.show()
