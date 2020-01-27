#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 21:12:28 2019

@author: Dartoon

Produce figures
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib as mat
import matplotlib as mpl
mat.rcParams['font.family'] = 'STIXGeneral'

import sys
sys.path.insert(0,'../../py_tools')

h=0.7

bhmass_overall=np.loadtxt('../Aklant/new_sample/log10_bh_mass_full_population.txt') - np.log10(h) 
bhmass_selected=np.loadtxt('../Aklant/new_sample/log10_bh_mass_selected_population.txt') - np.log10(h) 

mstar_overall=np.loadtxt('../Aklant/new_sample/log10_stellar_mass_full_population.txt') - np.log10(h) 
mstar_selected=np.loadtxt('../Aklant/new_sample/log10_stellar_mass_selected_population.txt') - np.log10(h) 

Lbol_overall=np.loadtxt('../Aklant/new_sample/log10_L_bol_full_population.txt') - np.log10(h) 
Lbol_selected=np.loadtxt('../Aklant/new_sample/log10_L_bol_selected_population.txt') - np.log10(h) 

logLedd_overall = 38. + np.log10(1.2) + bhmass_overall
logLedd_selected = 38. + np.log10(1.2) + bhmass_selected

plt.figure(figsize=(11,9))
#plt.scatter(bhmass_selected, Lbol_selected-logLedd_selected)

plt.hist2d(bhmass_overall,Lbol_overall-logLedd_overall,
                  norm=mpl.colors.LogNorm(),cmap='copper',bins=50,zorder=0,alpha=0.5)
cbar = plt.colorbar()
#cbar=f.colorbar(panel2[3],ax=obj)
cbar.ax.tick_params(labelsize=30) 
#
xspace = np.linspace(6,10)
plt.plot(xspace, 0*xspace,'k--',linewidth=3)
plt.plot(xspace, 0*xspace-1.5,'k--',linewidth=3)
y_line3 = -1.1*(xspace-7.5) -0.5
plt.plot(xspace, y_line3,'k--',linewidth=3)

yspace = np.linspace(-5,2)
plt.plot(yspace*0+7.5, yspace,'k--',linewidth=3)
plt.plot(xspace*0+8.5, yspace,'k--',linewidth=3)

plt.xlim([7.15,9.15])
plt.ylim([-3,1])

xfill = np.linspace(7.5, 8.5)
yfill_sline = -1.1*(xfill-7.5) -0.5
y_sline1 = xfill*0
y_sline2 = xfill*0-1.5
y4 = np.maximum(yfill_sline, y_sline2)
plt.fill_between(xfill, y4, y2=0, color='red', alpha='0.5', zorder=-1)

plt.tick_params(labelsize=30)
#ax.set_rasterized(True)
plt.ylabel(r"log(L$_{\rm bol}$/L$_{\rm Edd}$)",fontsize=30)
plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.savefig('MBII_selectfunc.pdf')
plt.show()