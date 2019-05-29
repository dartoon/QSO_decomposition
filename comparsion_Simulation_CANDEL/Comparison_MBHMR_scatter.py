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
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf
  
nll = lambda *args: -lnlike(*args)
#result = op.minimize(nll, [1.036, -1.947, 0.3], args=(x, y, yerr))
#m_ml, b_ml,sint_ml= result["x"]

from scipy.integrate import quad
h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)

#%%Calculate the scatter for the Aklant MBH-Mr
import scipy
h=0.7
#r_band_magnitudes_overall=np.load('Aklant/Aklant_scripts/r_band_magnitude_overall_population.npy')
r_band_magnitudes_selected=np.load('Aklant/Aklant_scripts/r_band_magnitude_selected_population.npy')
L_r = (0.4*(4.61-r_band_magnitudes_selected)) #LR in AB

##log10_bhmass_overall=np.load('Aklant/Aklant_scripts/log10_bhmass_overall_population.npy')
##log10_bhmass_overall -= np.log10(h)
log10_bhmass_selected=np.load('Aklant/Aklant_scripts/log10_bhmass_selected_population.npy')
log10_bhmass_selected -= np.log10(h)
plt.errorbar(L_r,log10_bhmass_selected,zorder=1,
             color='red',label='Simulated population',linestyle=' ',marker='o',ms=10,mec='k')
#
#result = op.minimize(nll, [0.58, 1.68, 0.3], 
#                     args=(L_r, log10_bhmass_selected, 0.5))
#a_alk, b_alk, sint_alk= result["x"]

#fit=scipy.optimize.curve_fit(lfit, log10_bhmass_selected, L_r)
#y_line = np.linspace(7, 9, 20)
#x_line = lfit(y_line,fit[0][0],fit[0][1])
fit=scipy.optimize.curve_fit(lfit, L_r, log10_bhmass_selected)
x_line = np.linspace(9, 12,20)
y_line = x_line*fit[0][0] + fit[0][1]

fit_err=np.sqrt(np.diag(fit[1]))
#a_alk, b_alk = fit[0][0],fit[0][1]

plt.plot(x_line, y_line)
#ty1= y_line + fit_err
#ty2= y_line - sint_alk
#plt.fill_between(x_line,ty1,ty2,color='linen',zorder=-50)

#%%

