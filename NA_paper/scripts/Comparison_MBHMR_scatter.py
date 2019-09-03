#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 08:45:53 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

#%%Define the function to fit theas linear
def lfit(x,m,c):
    return m*x+c

import scipy.optimize as op
def lnlike(theta, x, y, err):
    a, b, sint= theta
    model = a * x + b
    sigma2 = (err**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(sigma2))) 
    else:
      return -np.inf
nll = lambda *args: -lnlike(*args)
#result = op.minimize(nll, [1.036, -1.947, 0.3], args=(x, y, yerr))
#m_ml, b_ml,sint_ml= result["x"]

#from scipy.integrate import quad
#h0=70.             #km/s/Mpc
#om=0.3
#c=299790.        #speed of light [km/s]
#def EZ(z,om):
#      return 1/np.sqrt(om*(1+z)**3+(1-om))
#def EE(z):
#      return quad(EZ, 0, z , args=(om))[0]
#vec_EE=np.vectorize(EE)

#%%Calculate the scatter for the Aklant MBH-Mr
plt.figure(figsize=(8,6))
import scipy
h=0.7
#r_band_magnitudes_overall=np.load('Aklant/Aklant_scripts/r_band_magnitude_overall_population.npy')
r_band_magnitudes_selected=np.load('../Aklant/Aklant_scripts/r_band_magnitude_selected_population.npy')
L_r = (0.4*(4.61-r_band_magnitudes_selected)) #LR in AB

##log10_bhmass_overall=np.load('Aklant/Aklant_scripts/log10_bhmass_overall_population.npy')
##log10_bhmass_overall -= np.log10(h)
log10_bhmass_selected=np.load('../Aklant/Aklant_scripts/log10_bhmass_selected_population.npy')
log10_bhmass_selected -= np.log10(h)   #remove the [/h]
plt.errorbar(L_r,log10_bhmass_selected,zorder=1,
             color='red',label='Simulated population',linestyle=' ',marker='o',ms=10,mec='k')

##1. If fit using my way:
#result = op.minimize(nll, [0.58, 1.68, 0.3], 
#                     args=(L_r, log10_bhmass_selected, 0.5))
#a_alk, b_alk, sint_alk= result["x"]
##2. If fit y to x as function of f(y)
#fit=scipy.optimize.curve_fit(lfit, log10_bhmass_selected, L_r)
#y_line = np.linspace(7, 9, 20)
#x_line = lfit(y_line,fit[0][0],fit[0][1])
#a_alk, b_alk = fit[0][1],fit[0][0]
#print fit[0][1],fit[0][0]
##3. If fit x to y as function of f(x)
fit=scipy.optimize.curve_fit(lfit, L_r, log10_bhmass_selected)
x_line = np.linspace(9, 12,20)
a_alk, b_alk =  fit[0][0],fit[0][1]
#fit_err=np.sqrt(np.diag(fit[1]))  # To estimate the uncertainty of the fitting parameter (i.e. m and c in lfit)

y_line = x_line*a_alk + b_alk
plt.plot(x_line, y_line)
#Calculate the scatter:
print a_alk, b_alk
print "The scatter of Aklant sample:", np.sqrt(np.mean((log10_bhmass_selected - (L_r*a_alk + b_alk))**2))
#plt.title("Aklants sample")
plt.xlabel('L_R')
plt.ylabel('MBH')
plt.xlim(9.5,12)
plt.ylim(7.25, 8.75)
plt.show()
#%%
import linmix
x = L_r
xsig = L_r*0
y = log10_bhmass_selected
ysig = log10_bhmass_selected * 0

lm = linmix.LinMix(x, y, xsig=xsig, ysig=ysig, K=2)
lm.run_mcmc(silent=True)

#print("{}, {}".format(lm.chain['alpha'].mean(), lm.chain['alpha'].std()))
#print("{}, {}".format(lm.chain['beta'].mean(), lm.chain['beta'].std()))
#print("{}, {}".format(lm.chain['sigsqr'].mean(), lm.chain['sigsqr'].std()))

# The value of the intrinsic scatter:
sig_int = np.sqrt(lm.chain['sigsqr'].mean())
print sig_int

plt.errorbar(x,y, xerr=0, yerr=0, color='blue',ecolor='orange', fmt='.',zorder=-500,markersize=1)
for i in range(0, len(lm.chain), 25):
    xs = np.arange(8,13)
    ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
    plt.plot(xs, ys, color='r', alpha=0.02)
plt.xlim(9,12.5)
plt.ylim(6.0,10)
plt.show()




#%%
##%%
#import sys
#sys.path.insert(0,'../py_tools')
#from dmag import pass_dmag
#from load_result import load_host_p, load_MBH, load_err
#from load_result import load_zs, load_n
#ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
#'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']
#MB_ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'ECDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','LID1820','LID1622',\
#'LID1878', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']
#def lfit_fixm(x,c):
#    return a_alk*x+c
#
#zs = np.asarray(load_zs(ID))
#lumi_s = load_host_p(ID, dm = 0)[0] #!!! This dm is important 
#MBs = load_MBH(ID,MB_ID,if_reportHb=0   )
#LR_err = load_err(prop = 'LR', ID=ID)
#plt.errorbar(lumi_s,MBs, xerr=[np.abs(LR_err)[:,0], np.abs(LR_err)[:,1]], yerr=0.4, color='blue',ecolor='orange', fmt='.',zorder=-500,markersize=1)
#
#LR_err=(np.abs(LR_err)[:,0] + np.abs(LR_err)[:,1])/2.
#plt.scatter(lumi_s,MBs,c='green',s=580,marker="*",zorder=100, edgecolors='k')
#
#fit=scipy.optimize.curve_fit(lfit_fixm, lumi_s, MBs)
#x_line = np.linspace(9, 12,20)
#b_QSO =  fit[0][0]
#scatter =  np.sqrt(np.mean((MBs - (lumi_s*a_alk + b_QSO))**2))
#int_scatter = np.sqrt(scatter**2 - np.mean(LR_err**2) - 0.5**2)
#print "The scatter of QSO sample:", int_scatter
#
#plt.plot(x_line, y_line-b_alk + b_QSO)
#
#plt.xlabel('L_R')
#plt.ylabel('MBH')
#plt.show()
