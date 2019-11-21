#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 08:03:39 2019

@author: Dartoon

my_selection_for_SAM's sample
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../../py_tools')
import copy
from load_result import load_MBH, load_host_p, load_err, load_Lbol
import matplotlib as mpl
from scipy import stats
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
np.random.seed(5)  #12345

data_overall = np.loadtxt('../Nicola/realization/lagn_lhost1.dat')
#data_overall = np.loadtxt('../Nicola/realization/lagn_lhost_2D.dat')

bhmass_overall=data_overall[:,-1]
Lbol_overall=np.log10(data_overall[:,-2]*10**45.)
mstar_overall=data_overall[:,2]
magr_overall=data_overall[:,0]

logLedd_overall = 38. + np.log10(1.2) + bhmass_overall
#logLedd_selected = 38. + np.log10(1.2) + bhmass_selected

Eddr_overall = Lbol_overall-logLedd_overall

###Add noise to the data: 
#Noise level: MBH 0.4dex, mag_R 0.3mag, M* 0.17dex, Lbol 0.03dex
dMBH, dmag, dMstar, dLbol= 0.4, 0.3, 0.17, 0.03
#dMBH, dmag, dMstar, dLbol= 0.00004, 0.00003, 0.000017, 0.000003

bhmass_overall_noi = bhmass_overall + np.random.normal(0, dMBH, size=bhmass_overall.shape)
mstar_overall_noi = mstar_overall + np.random.normal(0, dMstar, size=mstar_overall.shape)
magr_overall_noi = magr_overall + np.random.normal(0, dmag, size=magr_overall.shape)
Lbol_overall_noi = Lbol_overall + np.random.normal(0, dLbol, size= Lbol_overall.shape)

logLedd_overall_noi = 38. + np.log10(1.2) + bhmass_overall_noi
Eddr_overall_noi = Lbol_overall_noi - logLedd_overall_noi

###Select sample:
select_window = (bhmass_overall_noi>7.7)*(bhmass_overall_noi<8.6)*(Eddr_overall_noi<0.0)*\
(Eddr_overall_noi > -1.1*(bhmass_overall_noi-7.5)-0.5 )* (Eddr_overall_noi>-1.5)
bhmass_select_noi = bhmass_overall_noi[select_window]
mstar_select_noi = mstar_overall_noi[select_window]
magr_select_noi = magr_overall_noi[select_window]
Lbol_select_noi = Lbol_overall_noi[select_window]
Eddr_select_noi = Eddr_overall_noi[select_window]

#Now can just plot the figure1 from figures_producer.py
plt.figure(figsize=(11,9))
plt.hist2d(bhmass_overall,Eddr_overall,norm=mpl.colors.LogNorm(),cmap='copper',bins=50,zorder=0,alpha=0.5)
plt.errorbar(bhmass_overall, Eddr_overall, c='gray',linestyle=' ',marker='.',ms=10,mec='k')

cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=30) 
plt.errorbar(bhmass_select_noi, Eddr_select_noi, c='green',linestyle=' ',marker='o',ms=10,mec='k')
xspace = np.linspace(6,10)
plt.plot(xspace, 0*xspace,'k--',linewidth=3)
plt.plot(xspace, 0*xspace-1.5,'k--',linewidth=3)
y_line3 = -1.1*(xspace-7.5) -0.5
plt.plot(xspace, y_line3,'k--',linewidth=3)
yspace = np.linspace(-5,2)
plt.plot(yspace*0+7.7, yspace,'k--',linewidth=3)
plt.plot(xspace*0+8.6, yspace,'k--',linewidth=3)

#Sign to the data to fit
bhmass_selected = bhmass_select_noi
mstar_selected = mstar_select_noi
r_band_magnitudes_selected = magr_select_noi
Lbol_selected = Lbol_select_noi
#logEddR_selected = Eddr_select_noi


# Import data and set up function:
tab_list = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 'CID3570',
            'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 'CID70', 'LID1273',
            'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396', 'CDFS-1', 'CDFS-229',
            'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']
ext_ID = {'XID2202':'LID1622', 'XID2138':'LID1820', 'XID2396':'LID1878', 'CDFS-321':'ECDFS-321'}
tab_sub_list = copy.deepcopy(tab_list)
for i in range(len(tab_sub_list)):
    if tab_sub_list[i] in ext_ID.keys():
        tab_sub_list[i] = ext_ID[tab_sub_list[i]]
Lr, M_star, M_r = load_host_p(tab_list, folder = '../../')  # M_star by Chabrier 

M_r_obs= M_r
M_r_obs_err = load_err(prop='LR', ID=tab_list) / 0.4

bh_mass_obs = load_MBH(tab_list,tab_sub_list, if_reportHb=0, folder = '../../') #+np.log10(h)
stellar_mass_obs=  np.log10(10**M_star / 0.54 * 0.91) #+np.log10(h)
stellar_mass_obs_err= load_err(prop='Mstar', ID=tab_list)

logLbol_obs = load_Lbol(ID=tab_list, folder = '../../')
logLedd_obs = 38. + np.log10(1.2) + bh_mass_obs
logEddR_obs = logLbol_obs - logLedd_obs

plt.errorbar(bh_mass_obs, logEddR_obs, c='orange',linestyle=' ',marker='o',ms=10,mec='k')
plt.xlim([7.15,9.15])
plt.ylim([-3,1])
xfill = np.linspace(7.7, 8.6)
yfill_sline = -1.1*(xfill-7.5) -0.5
y_sline1 = xfill*0
y_sline2 = xfill*0-1.5
y4 = np.maximum(yfill_sline, y_sline2)
plt.fill_between(xfill, y4, y2=0, color='green', alpha='0.5', zorder=-1)
plt.tick_params(labelsize=30)
plt.ylabel(r"log(L$_{\rm bol}$/L$_{\rm Edd}$)",fontsize=30)
plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
#plt.savefig('SAM_selectfunc.pdf')
plt.show()

#%%Plot MM data
def lfit(x,m,c):
    return m*x+c

import scipy.optimize


f,ax=plt.subplots(1,1,figsize=(11,10))

obj=ax
redshift=1.5
panel2=obj.hist2d(mstar_overall,bhmass_overall,
                  norm=mpl.colors.LogNorm(),cmap='copper',bins=50,zorder=0,alpha=0.5)
#panel2=obj.scatter(mstar_overall,bhmass_overall,c='gray',alpha=0.5,label='Simulated population')

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

cbar=f.colorbar(panel2[3],ax=obj)
cbar.ax.tick_params(labelsize=30) 

obj.errorbar(mstar_selected,bhmass_selected,zorder=1,
             color='green',label='SAM population',linestyle=' ',marker='o',ms=10,mec='k')
obj.errorbar(stellar_mass_obs,bh_mass_obs, 
#             xerr = [abs(M_r_obs_err[:,0]),abs(M_r_obs_err[:,1])], yerr=np.ones(len(bh_mass_obs))*0.4, 
             zorder=100,color='orange',label='Observed population',
             linestyle=' ',marker='o',ms=10,mec='k')

##Fit y as function of x
#x, y = mstar_selected, bhmass_selected
#fit=scipy.optimize.curve_fit(lfit, x, y)
#fit_err=np.sqrt(np.diag(fit[1]))
#x_space=np.linspace(-50,50,100)
#y_space=lfit(x_space,fit[0][0],fit[0][1])
#plt.plot(x_space,y_space,color='green',linewidth=3)
#y_space_ub=lfit(x_space,fit[0][0]-fit_err[0]/2,fit[0][1]-fit_err[1]/2)
#y_space_lb=lfit(x_space,fit[0][0]+fit_err[0]/2,fit[0][1]+fit_err[1]/2)
#plt.fill_between(x_space,y_space_lb,y_space_ub,color='green',alpha=0.15)
#def lfit_fixm(x,c):
#    m_0 = fit[0][0]
#    return m_0*x+c
#x_obs, y_obs = stellar_mass_obs, bh_mass_obs
#fit_fixm=scipy.optimize.curve_fit(lfit_fixm, x_obs, y_obs)
##fit_err=np.sqrt(np.diag(fit[1]))
#y_obs_space=lfit_fixm(x_space,fit_fixm[0])
#plt.plot(x_space,y_obs_space,color='green',linewidth=3)
#print "mismatch:", fit_fixm[0]- fit[0][1]

#Fit x as function of y
x, y = mstar_selected, bhmass_selected
fit=scipy.optimize.curve_fit(lfit, y, x)
fit_err=np.sqrt(np.diag(fit[1]))
y_space=np.linspace(-50,50,100)
x_space=lfit(y_space,fit[0][0],fit[0][1])   #y_space become to x_space
plt.plot(x_space, y_space, color='green',linewidth=3)
x_space_ub=lfit(y_space,fit[0][0]+fit_err[0]/2,fit[0][1]+fit_err[1]/2)
x_space_lb=lfit(y_space,fit[0][0]-fit_err[0]/2,fit[0][1]-fit_err[1]/2)
#plt.fill_betweenx(y_space,x_space_lb,x_space_ub,color='green',alpha=0.15)
def lfit_fixm(y,c):
    m_0 = fit[0][0]
    return m_0*y+c
x_obs, y_obs = stellar_mass_obs, bh_mass_obs
fit_fixm=scipy.optimize.curve_fit(lfit_fixm, y_obs, x_obs)
#fit_err=np.sqrt(np.diag(fit[1]))
x_obs_space=lfit_fixm(y_space,fit_fixm[0])
plt.plot(x_obs_space,y_space,color='orange',linewidth=3)
print "mismatch:", fit_fixm[0]- fit[0][1]  #In BH mass offset space

# Plot the fitting scatter 
x_space_ub=lfit(y_space,fit[0][0],fit[0][1]+0.72)     #Values are measured in the next box
x_space_lb=lfit(y_space,fit[0][0],fit[0][1]-0.72)
plt.fill_betweenx(y_space,x_space_lb,x_space_ub,color='green',alpha=0.35)
#print "mismatch:", fit_fixm[0]- fit[0][1]  #In BH mass offset space
x_space_ub=lfit_fixm(y_space,fit_fixm[0]+0.3247)
x_space_lb=lfit_fixm(y_space,fit_fixm[0]-0.3247)
plt.fill_betweenx(y_space,x_space_lb,x_space_ub,color='orange',alpha=0.35)

obj.set_yticks([7.5,8.0,8.5,9.0])
obj.set_xticks([10,10.5,11,11.5,12])
#obj.set_xticklabels(['-18','-20','-22','-24','-26'])
ax.set_xlim(9.7,11.9)  #
ax.set_ylim(7.2, 9.4)  #
obj.tick_params(labelsize=30)
#ax.set_rasterized(True)
obj.set_ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
obj.set_xlabel('log(M$_{*}$/M$_{\odot}$)',fontsize=30)
obj.legend(loc='upper left',fontsize=21,numpoints=1)
#plt.savefig("SAM_MM_consider_nois.pdf")
plt.show()

#%%
'''
import linmix
x = mstar_selected
xsig = np.zeros(len(bhmass_selected))+dMstar
y = bhmass_selected
ysig = np.zeros(len(bhmass_selected))+dMBH
lm = linmix.LinMix(x, y, xsig=xsig, ysig=ysig, K=3)
lm.run_mcmc(silent=True)
alpha = lm.chain['alpha'].mean()
beta = lm.chain['beta'].mean()
xs = np.arange(-26,-19)
ys = alpha + xs * beta
#plt.plot(xs, ys, color='green',linewidth=3)
print "intrinsic scatter:", np.sqrt(lm.chain['sigsqr'].mean()), np.sqrt(lm.chain['sigsqr'].std())
##Don't know how to fix the slope value...
x = stellar_mass_obs
xsig = (abs(stellar_mass_obs_err[:,0]) + abs(stellar_mass_obs_err[:,1]))/2
y = bh_mass_obs
ysig = np.ones(len(bh_mass_obs))*0.4
lm_obs = linmix.LinMix(x, y, xsig=xsig, ysig=ysig, K=3)
lm_obs.run_mcmc(silent=True)
alpha_obs = lm_obs.chain['alpha'].mean()
beta_obs = lm_obs.chain['beta'].mean()
xs = np.arange(-26,-19)
ys = alpha_obs + xs * beta_obs
#plt.plot(xs, ys, color='green',linewidth=3)
print "intrinsic scatter:", np.sqrt(lm_obs.chain['sigsqr'].mean()), np.sqrt(lm_obs.chain['sigsqr'].std())
'''
#%%Plot the 1-D scatter for MM.
plt.figure(figsize=(8,7))
plt.hist(stellar_mass_obs - lfit_fixm(bh_mass_obs,fit_fixm[0]), histtype=u'step',normed=True,
         label=('HST sample'), linewidth = 2, color='orange')
plt.hist(mstar_selected - lfit(bhmass_selected,fit[0][0],fit[0][1]),histtype=u'step',normed=True,
         label=('SAM sample'), linewidth = 2, color='green')
plt.title(r"The offset comparison for the M$_{\rm BH}$-M$_{*}$ relation", fontsize = 20)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.yticks([])
plt.xlabel('$\Delta$log(M$_{*}$/M$_{\odot}$)',fontsize=30)
#plt.savefig('comp_scatter_MM.pdf')
plt.show()

print "sim scatter:", np.std(mstar_selected - lfit(bhmass_selected,fit[0][0],fit[0][1]))
print "obs scatter:", np.std(stellar_mass_obs - lfit_fixm(bh_mass_obs,fit_fixm[0]))
print "KS scatter:", stats.ks_2samp((mstar_selected - lfit(bhmass_selected,fit[0][0],fit[0][1])),
                                    (stellar_mass_obs - lfit_fixm(bh_mass_obs,fit_fixm[0]))).pvalue

#%%Plot ML data
f,ax=plt.subplots(1,1,figsize=(11,10))

obj=ax
redshift=1.5
panel2=obj.hist2d(magr_overall,bhmass_overall,
                  norm=mpl.colors.LogNorm(),cmap='copper',bins=50,zorder=0,alpha=0.5)
#panel2=obj.scatter(mstar_overall,bhmass_overall,c='gray',alpha=0.5,label='Simulated population')

####Fit the overall sample (x as function of y):
#fit_1=scipy.optimize.curve_fit(lfit,-r_band_magnitudes_overall,log10_bhmass_overall)
#fit_err=np.sqrt(np.diag(fit_1[1]))
#r_band_space=np.linspace(20,26,100)
#lmbh_space=lfit(r_band_space,fit_1[0][0],fit_1[0][1])
#plt.plot(r_band_space,lmbh_space,color='blue',linewidth=3)
###Not corrected yet:
##r_band_space_ub=lfit(lmbh_space,fit_1[0][0]+fit_err[0]/2,fit_1[0][1]+fit_err[1]/2)
##r_band_space_lb=lfit(lmbh_space,fit_1[0][0]-fit_err[0]/2,fit_1[0][1]-fit_err[1]/2)
##plt.fill_betweenx(lmbh_space,r_band_space_lb,r_band_space_ub,color='blue',alpha=0.3)

cbar=f.colorbar(panel2[3],ax=obj)
cbar.ax.tick_params(labelsize=30) 

obj.errorbar(r_band_magnitudes_selected,bhmass_selected,zorder=1,
             color='green',label='SAM population',linestyle=' ',marker='o',ms=10,mec='k')
obj.errorbar(M_r_obs,bh_mass_obs, 
#             xerr = [abs(M_r_obs_err[:,0]),abs(M_r_obs_err[:,1])], yerr=np.ones(len(bh_mass_obs))*0.4, 
             zorder=100,color='orange',label='Observed population',
             linestyle=' ',marker='o',ms=10,mec='k')

##Fit y as function of x
#x, y = r_band_magnitudes_selected, bhmass_selected
#fit_1=scipy.optimize.curve_fit(lfit,x, y)
#fit_err_1=np.sqrt(np.diag(fit_1[1]))
#x_space=np.linspace(-26,-18,100)
#y_space=lfit(x_space,fit_1[0][0],fit_1[0][1])
#plt.plot(x_space,y_space,color='green',linewidth=3)
#lmbh_space_ub=lfit(x_space,fit_1[0][0]-fit_err_1[0]/2,fit_1[0][1]+fit_err_1[1]/2)
#lmbh_space_lb=lfit(x_space,fit_1[0][0]+fit_err_1[0]/2,fit_1[0][1]-fit_err_1[1]/2)
#plt.fill_between(x_space,lmbh_space_lb,lmbh_space_ub,color='green',alpha=0.15)
#def lfit_fixm_1(x,c):
#    return fit_1[0][0]*x+c
#fit_fixm_1=scipy.optimize.curve_fit(lfit_fixm_1,M_r_obs, bh_mass_obs)
##fit_err_1=np.sqrt(np.diag(fit_1[1]))
#lmbh_space=lfit_fixm_1(x_space,fit_fixm_1[0])
#plt.plot(x_space,lmbh_space,color='green',linewidth=3)
#print "mismatch:", fit_fixm_1[0]- fit_1[0][1]
#Fit x as function of y
x, y = r_band_magnitudes_selected, bhmass_selected
fit_1=scipy.optimize.curve_fit(lfit,y, x)
fit_err_1=np.sqrt(np.diag(fit_1[1]))
y_space=np.linspace(5,13,100)
x_space=lfit(y_space,fit_1[0][0],fit_1[0][1])
plt.plot(x_space,y_space,color='green',linewidth=3)
x_space_ub=lfit(y_space,fit_1[0][0]+fit_err_1[0]/2,fit_1[0][1]+fit_err_1[1]/2)
x_space_lb=lfit(y_space,fit_1[0][0]-fit_err_1[0]/2,fit_1[0][1]-fit_err_1[1]/2)
#plt.fill_betweenx(y_space,x_space_lb,x_space_ub,color='green',alpha=0.15)
def lfit_fixm_1(y,c):
    return fit_1[0][0]*y+c
x_obs, y_obs = M_r_obs, bh_mass_obs
fit_fixm_1=scipy.optimize.curve_fit(lfit_fixm_1, y_obs, x_obs)
#fit_err_1=np.sqrt(np.diag(fit_1[1]))
x_space_obs=lfit_fixm_1(y_space,fit_fixm_1[0])
plt.plot(x_space_obs,y_space,color='orange',linewidth=3)

# Plot the fitting scatter 
#plt.plot(x_space, y_space, color='steelblue',linewidth=3)
x_space_ub=lfit(y_space,fit_1[0][0],fit_1[0][1]+1.47775)
x_space_lb=lfit(y_space,fit_1[0][0],fit_1[0][1]-1.47775)
plt.fill_betweenx(y_space,x_space_lb,x_space_ub,color='green',alpha=0.35)
#print "mismatch:", fit_fixm[0]- fit[0][1]  #In BH mass offset space
x_space_ub=lfit_fixm_1(y_space,fit_fixm_1[0]+0.7222)
x_space_lb=lfit_fixm_1(y_space,fit_fixm_1[0]-0.7222)
plt.fill_betweenx(y_space,x_space_lb,x_space_ub,color='orange',alpha=0.35)


print "\n\nPlot M-Mag relation:"
print "mismatch:", fit_fixm_1[0]- fit_1[0][1]

obj.set_yticks([7.5,8.0,8.5,9.0])
obj.set_xticks([-20,-21, -22, -23, -24, -25])
#obj.set_xticklabels(['-18','-20','-22','-24','-26'])
ax.set_xlim(-19.8, -25.5)  # 
ax.set_ylim(7.2, 9.4)  # 

obj.tick_params(labelsize=30)
#ax.set_rasterized(True)
obj.set_ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
obj.set_xlabel('R band magnitude',fontsize=30)
obj.legend(loc='upper left',fontsize=21,numpoints=1)
#plt.savefig("SAM_ML_consider_nois.pdf")
plt.show()

#%%Plot the 1-D scatter for ML.
plt.figure(figsize=(8,7))
plt.hist(M_r_obs - lfit_fixm_1(bh_mass_obs,fit_fixm_1[0]), histtype=u'step',normed=True,
         label=('HST sample'), linewidth = 2, color='orange')
plt.hist(r_band_magnitudes_selected - lfit(bhmass_selected,fit_1[0][0],fit_1[0][1]),histtype=u'step',normed=True,
         label=('SAM sample'), linewidth = 2, color='green')
plt.title(r"The offset comparison for the M$_{\rm BH}$-mag relation", fontsize = 20)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.yticks([])
plt.xlabel(r"$\Delta$magnitude",fontsize=30)
#plt.savefig('comp_scatter_ML.pdf')
plt.show()

print "sim scatter:", np.std(r_band_magnitudes_selected - lfit(bhmass_selected,fit_1[0][0],fit_1[0][1]))
print "obs scatter:", np.std(M_r_obs - lfit_fixm_1(bh_mass_obs,fit_fixm_1[0]))
print "KS:", stats.ks_2samp((r_band_magnitudes_selected - lfit(bhmass_selected,fit_1[0][0],fit_1[0][1])),
                                    (M_r_obs - lfit_fixm_1(bh_mass_obs,fit_fixm_1[0]))).pvalue



##%%Plot the 1-D hist for Mstar, R_Mag and MBH and do the K-S test in 1D.
#
#plt.figure(figsize=(8,6))
#plt.hist(bhmass_selected ,histtype=u'step',normed=True,
#         label=('SAM BH sample'), linewidth = 2, color='orange')
#plt.hist(bh_mass_obs , histtype=u'step',normed=True,
#         label=('HST BH sample'), linewidth = 2, color='green')
#plt.tick_params(labelsize=20)
#plt.legend(prop={'size':20})
#plt.yticks([])
#plt.show()
#print stats.ks_2samp(bhmass_selected, bh_mass_obs).pvalue
#
#plt.figure(figsize=(8,6))
#plt.hist(mstar_selected ,histtype=u'step',normed=True,
#         label=('SAM M* sample'), linewidth = 2, color='orange')
#plt.hist(stellar_mass_obs , histtype=u'step',normed=True,
#         label=('HST M* sample'), linewidth = 2, color='green')
#plt.tick_params(labelsize=20)
#plt.legend(prop={'size':20})
#plt.yticks([])
#plt.show()
#print stats.ks_2samp(mstar_selected, stellar_mass_obs).pvalue
#
#plt.figure(figsize=(8,6))
#plt.hist(r_band_magnitudes_selected ,histtype=u'step',normed=True,
#         label=('SAM MagR sample'), linewidth = 2, color='orange')
#plt.hist(M_r_obs , histtype=u'step',normed=True,
#         label=('HST MagR sample'), linewidth = 2, color='green')
#plt.tick_params(labelsize=20)
#plt.legend(prop={'size':20})
#plt.yticks([])
#plt.show()
#print stats.ks_2samp(r_band_magnitudes_selected, M_r_obs).pvalue


##%% Estiamte the in
#import linmix
#x = mstar_overall
#xsig = np.zeros(len(mstar_overall))
#y = bhmass_overall
#ysig = np.zeros(len(mstar_overall))
#lm = linmix.LinMix(x, y, xsig=xsig, ysig=ysig, K=3)
#lm.run_mcmc(silent=True)
#alpha = lm.chain['alpha'].mean()
#beta = lm.chain['beta'].mean()
#xs = np.arange(5,15)
#ys = alpha + xs * beta
#plt.scatter(x, y)
#plt.plot(xs, ys, color='green',linewidth=3)
#print "intrinsic scatter:", np.sqrt(lm.chain['sigsqr'].mean()), np.sqrt(lm.chain['sigsqr'].std())
#
##%%To plot data and plot together with MBII
#
#plt.figure(figsize=(8,7))
#plt.hist(M_r_obs - lfit_fixm_1(bh_mass_obs,fit_fixm_1[0]), histtype=u'step',normed=True,
#         label=('Observed sample'), linewidth = 2, color='orange')
#plt.hist(MBII_scatter_ML,histtype=u'step',normed=True,  #Run the "my_selection_plots_MBII.py" first
#         label=('MBII sample'), linewidth = 2, color='steelblue')
#plt.hist(r_band_magnitudes_selected - lfit(bhmass_selected,fit_1[0][0],fit_1[0][1]),histtype=u'step',normed=True,
#         label=('SAM sample'), linewidth = 2, color='green')
#plt.title(r"Scatter comparison for the M$_{\rm BH}$-mag relation", fontsize = 20)
#plt.tick_params(labelsize=20)
#plt.legend(prop={'size':20})
#plt.yticks([])
#plt.xlabel(r"$\Delta$magnitude",fontsize=30)
#plt.savefig('comp_scatter_ML.pdf')
#plt.show()
#
#
#
#plt.figure(figsize=(8,7))
#plt.hist(stellar_mass_obs - lfit_fixm(bh_mass_obs,fit_fixm[0]), histtype=u'step',normed=True,
#         label=('Observed sample'), linewidth = 2, color='orange')
#plt.hist(MBII_scatter_MM,histtype=u'step',normed=True,  #Run the "my_selection_plots_MBII.py" first
#         label=('MBII sample'), linewidth = 2, color='steelblue')
#plt.hist(mstar_selected - lfit(bhmass_selected,fit[0][0],fit[0][1]),histtype=u'step',normed=True,
#         label=('SAM sample'), linewidth = 2, color='green')
#plt.title(r"Scatter comparison for the M$_{\rm BH}$-M$_{*}$ relation", fontsize = 20)
#plt.tick_params(labelsize=20)
#plt.legend(prop={'size':20},loc=2)
#plt.yticks([])
#plt.xlabel('$\Delta$log(M$_{*}$/M$_{\odot}$)',fontsize=30)
#plt.savefig('comp_scatter_MM.pdf')
#plt.show()