#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 20:15:08 2019

@author: Dartoon

Plot the MBHII MBH-Mr relation and fit their scatter.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy
import sys

import matplotlib.lines as mlines
from matplotlib import colors
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

sys.path.insert(0,'../../py_tools')
from load_result import load_MBH, load_host_p, load_err

import matplotlib as mpl
#data=np.loadtxt('../AGN_host_inference_for_npy.txt')
#print data
h=0.7

tab_list = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 'CID3570',
            'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 'CID70', 'LID1273',
            'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396', 'CDFS-1', 'CDFS-229',
            'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']
ext_ID = {'XID2202':'LID1622', 'XID2138':'LID1820', 'XID2396':'LID1878', 'CDFS-321':'ECDFS-321'}
tab_sub_list = copy.deepcopy(tab_list)
for i in range(len(tab_sub_list)):
    if tab_sub_list[i] in ext_ID.keys():
        tab_sub_list[i] = ext_ID[tab_sub_list[i]]
Lr, M_star, M_r = load_host_p(tab_list,temp='3Gyrs', folder = '../../')  # M_star by Chabrier 

bh_mass_obs = load_MBH(tab_list,tab_sub_list, if_reportHb=0, folder = '../../') #+np.log10(h)
stellar_mass_obs=  np.log10(10**M_star / 0.54 * 0.91) #+np.log10(h)
stellar_mass_obs_err= load_err(prop='Mstar', ID=tab_list)

#M_r_obs= M_r
#M_r_obs_err = load_err(prop='LR', ID=tab_list) / 0.4
#r_band_magnitudes_overall=np.loadtxt('../Aklant/new_sample/log10_host_r_mag_full_population.txt')
#r_band_magnitudes_selected=np.loadtxt('../Aklant/new_sample/log10_host_r_mag_selected_population.txt')

bhmass_overall=np.loadtxt('../Aklant/new_sample/log10_bh_mass_full_population.txt') - np.log10(h) 
bhmass_selected=np.loadtxt('../Aklant/new_sample/log10_bh_mass_selected_population.txt') - np.log10(h) 

mstar_overall=np.loadtxt('../Aklant/new_sample/log10_stellar_mass_full_population.txt') - np.log10(h) 
mstar_selected=np.loadtxt('../Aklant/new_sample/log10_stellar_mass_selected_population.txt') - np.log10(h) 

#%% Plot it out:
def lfit(x,m,c):
    return m*x+c

import scipy.optimize
f,ax=plt.subplots(1,1,figsize=(8,8))

obj=ax
redshift=1.5
#panel2=obj.hist2d(mstar_overall,bhmass_overall,
#                  norm=mpl.colors.LogNorm(),cmap='Blues_r',bins=50,zorder=0,alpha=0.5)
panel2=obj.scatter(mstar_overall,bhmass_overall,c='gray',alpha=0.5,label='Simulated population')

####Fit the overall sample (x as function of y):
#fit=scipy.optimize.curve_fit(lfit,-r_band_magnitudes_overall,log10_bhmass_overall)
#fit_err=np.sqrt(np.diag(fit[1]))
#r_band_space=np.linspace(20,26,100)
#lmbh_space=lfit(r_band_space,fit[0][0],fit[0][1])
#plt.plot(r_band_space,lmbh_space,color='blue',linewidth=3)
###Not corrected yet:
##r_band_space_ub=lfit(lmbh_space,fit[0][0]+fit_err[0]/2,fit[0][1]+fit_err[1]/2)
##r_band_space_lb=lfit(lmbh_space,fit[0][0]-fit_err[0]/2,fit[0][1]-fit_err[1]/2)
##plt.fill_betweenx(lmbh_space,r_band_space_lb,r_band_space_ub,color='blue',alpha=0.3)

#cbar=f.colorbar(panel2[3],ax=obj)
#cbar.ax.tick_params(labelsize=30) 

obj.errorbar(mstar_selected,bhmass_selected,zorder=1,
             color='red',label='Selected population',linestyle=' ',marker='o',ms=10,mec='k')
obj.errorbar(stellar_mass_obs,bh_mass_obs, 
#             xerr = [abs(M_r_obs_err[:,0]),abs(M_r_obs_err[:,1])], yerr=np.ones(len(bh_mass_obs))*0.4, 
             zorder=1,color='green',label='Observed population',
             linestyle=' ',marker='o',ms=10,mec='k')

fit=scipy.optimize.curve_fit(lfit,mstar_selected, bhmass_selected)
fit_err=np.sqrt(np.diag(fit[1]))
mstar_space=np.linspace(9,13,100)
lmbh_space=lfit(mstar_space,fit[0][0],fit[0][1])
plt.plot(mstar_space,lmbh_space,color='red',linewidth=3)
lmbh_space_ub=lfit(mstar_space,fit[0][0]-fit_err[0]/2,fit[0][1]-fit_err[1]/2)
lmbh_space_lb=lfit(mstar_space,fit[0][0]+fit_err[0]/2,fit[0][1]+fit_err[1]/2)
plt.fill_between(mstar_space,lmbh_space_lb,lmbh_space_ub,color='red',alpha=0.15)

def lfit_fixm(x,c):
    return fit[0][0]*x+c
fit_fixm=scipy.optimize.curve_fit(lfit_fixm,stellar_mass_obs, bh_mass_obs)
#fit_err=np.sqrt(np.diag(fit[1]))
lmbh_space=lfit_fixm(mstar_space,fit_fixm[0])
plt.plot(mstar_space,lmbh_space,color='green',linewidth=3)

print "mismatch:", fit_fixm[0]- fit[0][1]

#import linmix
#x = mstar_selected
#xsig = np.zeros(len(bhmass_selected))
#y = bhmass_selected
#ysig = np.zeros(len(bhmass_selected))
#lm = linmix.LinMix(x, y, xsig=xsig, ysig=ysig, K=3)
#lm.run_mcmc(silent=True)
#alpha = lm.chain['alpha'].mean()
#beta = lm.chain['beta'].mean()
#xs = np.arange(-26,-19)
#ys = alpha + xs * beta
##plt.plot(xs, ys, color='red',linewidth=3)
#print "intrinsic scatter:", np.sqrt(lm.chain['sigsqr'].mean())
###Don't know how to fix the slope value...
#x = stellar_mass_obs
#xsig = (abs(stellar_mass_obs_err[:,0]) + abs(stellar_mass_obs_err[:,1]))/2
#y = bh_mass_obs
#ysig = np.ones(len(bh_mass_obs))*0.4
#lm_obs = linmix.LinMix(x, y, xsig=xsig, ysig=ysig, K=3)
#lm_obs.run_mcmc(silent=True)
#alpha_obs = lm_obs.chain['alpha'].mean()
#beta_obs = lm_obs.chain['beta'].mean()
#xs = np.arange(-26,-19)
#ys = alpha_obs + xs * beta_obs
##plt.plot(xs, ys, color='green',linewidth=3)
#print "intrinsic scatter:", np.sqrt(lm_obs.chain['sigsqr'].mean())


obj.set_yticks([7.5,8.0,8.5,9.0])
obj.set_xticks([10,10.5,11,11.5,12])
#obj.set_xticklabels(['-18','-20','-22','-24','-26'])
ax.set_xlim(9.7,12)  # decreasing time
ax.set_ylim(7.2, 9.4)  # decreasing time

obj.tick_params(labelsize=30)
#ax.set_rasterized(True)
obj.set_ylabel('log(M$_{BH}$/M$_{\odot}$)',fontsize=30)
obj.set_xlabel('log(M$_{*}$/M$_{\odot}$)',fontsize=30)
obj.legend(loc='upper left',fontsize=17,numpoints=1)
plt.savefig("MBII_MM_.pdf")
plt.show()