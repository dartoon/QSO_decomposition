#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 20:15:08 2019

@author: Dartoon

Plot the MBHII MBH-Mr relation and fit their scatter.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

h=0.7

bhmass_overall=np.loadtxt('../Aklant/new_sample/log10_bh_mass_full_population.txt') - np.log10(h) 
#bhmass_selected=np.loadtxt('../Aklant/new_sample/log10_bh_mass_selected_population.txt') - np.log10(h) 

mstar_overall=np.loadtxt('../Aklant/new_sample/log10_stellar_mass_full_population.txt') - np.log10(h) 
#mstar_selected=np.loadtxt('../Aklant/new_sample/log10_stellar_mass_selected_population.txt') - np.log10(h) 

magr_overall=np.loadtxt('../Aklant/new_sample/log10_host_r_mag_full_population.txt')
#r_band_magnitudes_selected=np.loadtxt('../Aklant/new_sample/log10_host_r_mag_selected_population.txt')

Lbol_overall=np.loadtxt('../Aklant/new_sample/log10_L_bol_full_population.txt') - np.log10(h) 
#Lbol_selected=np.loadtxt('../Aklant/new_sample/log10_L_bol_selected_population.txt') - np.log10(h) 

logLedd_overall = 38. + np.log10(1.2) + bhmass_overall
#logLedd_selected = 38. + np.log10(1.2) + bhmass_selected

Eddr_overall = Lbol_overall-logLedd_overall

#%% Estiamte the in
import linmix
x = mstar_overall
xsig = np.zeros(len(mstar_overall))
y = bhmass_overall
ysig = np.zeros(len(mstar_overall))
lm = linmix.LinMix(x, y, xsig=xsig, ysig=ysig, K=3)
lm.run_mcmc(silent=True)
alpha = lm.chain['alpha'].mean()
beta = lm.chain['beta'].mean()
xs = np.arange(5,15)
ys = alpha + xs * beta
plt.scatter(x, y)
plt.plot(xs, ys, color='red',linewidth=3)
print "intrinsic scatter:", np.sqrt(lm.chain['sigsqr'].mean()), np.sqrt(lm.chain['sigsqr'].std())
