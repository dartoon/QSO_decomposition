#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 11:21:07 2019

@author: Dartoon

The selection effect study by Xuheng
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy

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
mbh_grid = np.linspace(7.4,10, 100)
lam_grid = np.linspace(-2, 0.25, 101)
rho_lam_vs_mbh_list = []
for i in range(len(mbh_grid)):
    rho_lam_vs_mbh_list.append(vec_rho_lam_2d(lam_grid, mbh_grid[i]))
rho_lam_vs_mbh = np.asarray(rho_lam_vs_mbh_list)
import matplotlib
my_cmap = copy.copy(matplotlib.cm.get_cmap('Oranges',10)) # copy the default cmap
rho_2d = np.asarray([(rho_lam_vs_mbh[i]+rho_bh(mbh_grid[i])) for i in range(len(mbh_grid))])
for i in range(len(mbh_grid)):
#    print mbh[i], rho_bh(mbh[i])
    plt.scatter(lam_grid*0+mbh_grid[i], lam_grid, c= rho_2d[i], s = 140, vmin = np.min(rho_2d), vmax = np.max(rho_2d),
                   marker='s', alpha=0.9, edgecolors='none', cmap = my_cmap)
cl=plt.colorbar()
cl.set_label('Possibility density', size=20)
plt.show()

#%%Define the mock local sample and assuming the high redshift sample follows the sample.
#Demonstrate the way to generate the local sample:

def model(theta, x):
    """
    The corresponding MBH for a given M*(x)
    """
    m, b= theta
    model = m * x + b
    return model
def de_model(theta, y):
    """
    solve the model
    """
    m, b= theta
    x = (y-b) / m
    return x   

m_ml, b_ml, sint = [0.93, -1.96, 0.3]
x = np.linspace(8.5, 13, 200)
y = model([m_ml, b_ml], x)
plt.figure(1, figsize=(7,8))
plt.plot(x,y, 'k--')
#y_noised = model([m_ml, b_ml, sint], x, scatter = True)
#plt.plot(x,y_noised,'r')

#sampling the local sample:
local_x = np.random.uniform(9.5, 12.0, 30)
local_y = model([m_ml, b_ml], local_x)
local_x += np.random.normal(0, (0.15**2 + sint**2)**0.5, size=local_x.shape)
local_y += np.random.normal(0, 0.2, size=local_y.shape)
plt.errorbar(local_x, local_y, xerr=0.3, yerr=0.4, fmt='.',color='black')

#%%
#Simualte the high redshfit sample:
mbh_grid = np.linspace(7.5,9.5,100)
mbh_rho = rho_bh(mbh_grid)
#normize mbh_rho and calculate the CDF distribution for the mock data generate:
mbh_rho = 10**mbh_rho
mbh_grid_s = (mbh_grid[-1]-mbh_grid[1])/len(mbh_grid)
mbh_rho /= np.sum(mbh_rho) * mbh_grid_s
mbh_rho_int = np.array([np.sum((mbh_rho*mbh_grid_s)[:i]) for i in range(len(mbh_rho))])
# Random simulate the MBH given the BHMF
vol = 10000
mbh, lam = [], []
for i in range(vol):
    idx = np.sum(np.random.random()>mbh_rho_int)-1
    mbh.append(mbh_grid[idx]) #np.random.uniform(R[idx, 0],R[idx+1, 0])
for i in range(len(mbh)):
    lam_rho = 10**rho_lam_2d(lam_grid, mbh[i])
    lam_grid_s = (lam_grid[-1]-lam_grid[1])/len(lam_grid)  
    lam_rho /= np.sum(lam_rho) * lam_grid_s
    lam_rho_int = np.array([np.sum((lam_rho*lam_grid_s)[:i]) for i in range(len(lam_rho))])
    idx = np.sum(np.random.random()>mbh_rho_int)-1
    lam.append(lam_grid[idx])

mbh = np.array([mbh]) 
mbh += np.random.normal(0, 0.4, size=mbh.shape)
mstar = de_model([m_ml, b_ml], mbh)
mstar += np.random.normal(0, (0.3**2+ sint**2)**0.5, size=mstar.shape)
lam = np.array([lam])
bools = [ (mbh > 7.5) * (mbh < 8.5) * (lam < 0) * (lam >  -1.1*(mbh-7.5)-0.5)]
mbh_select = mbh[bools[0]]     
lam_select = lam[bools[0]]  
mstar_select = mstar[bools[0]]
plt.errorbar(mstar, mbh,xerr=0.0, yerr=0.0, fmt='.',color='gray',zorder=0)
plt.errorbar(mstar_select[:], mbh_select[:],xerr=0.0, yerr=0.0, fmt='.',color='orange',zorder=1)

plt.text(10., 6.5, "log$(M_{BH}/10^{7}M_{\odot})$=%s+%slog$(M_*/10^{10}M_{\odot})$"%(round(b_ml+m_ml*10-7,2),round(m_ml,2)),
         color='blue',fontsize=15)
plt.xlabel("log$(M_{*}/M_{\odot})$",fontsize=15)
plt.ylabel("log$(M_{BH}/M_{\odot})$", fontsize=15) 
plt.xlim([9.0,12.5])
plt.tick_params(labelsize=15)
plt.show()

#%%Selectoin window
plt.errorbar(mbh, lam,xerr=0.0, yerr=0.0, fmt='.',color='gray',zorder=0)
plt.errorbar(mbh_select, lam_select,xerr=0.0, yerr=0.0, fmt='.',color='orange',zorder=1)
plt.show()

