#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:13:37 2019

@author: Dartoon

Plot the n-B/T relation
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

data_all = np.loadtxt('bulge_disk_result.txt')
# 0ID, 1RA, 2DEC, 3FIELD(GDS:0, gdn:1, COSMOS:2, UDS:3, egs: 4), 4B_T_m, 5MAG_F125w, 6RE_F125w, 7single_Sersic_n_F125w, 8LogMass 
data = data_all[data_all[:,3] == 2]

#The COSMOS files in http://www.mpia.de/homes/vdwel/candels.html by van der Wel et al. (2012).
#   NUMBER         RA        DEC          f        mag       dmag         re        dre          n         dn          q         dq         pa        dpa          sn
galfit_loc = np.loadtxt('../cos_2epoch_wfc3_f160w_060mas_v1.0_galfit.cat')[:,[1,2]]  #RA, DEC
galfit = np.loadtxt('../cos_2epoch_wfc3_f160w_060mas_v1.0_galfit.cat')[:,[4,6,8,3]] # mag, re, n, flag

##The data from 3D HST
stellar_loc = np.loadtxt('../cosmos_3dhst.v4.1.cat')[:,[3,4]]  # RA, DEC
#stellar_flux_ap = np.loadtxt('../cosmos_3dhst.v4.1.cat')[:,[69,51,39]]  # flux, F140w, F125w, F814w.
stellar = np.loadtxt('../cosmos_3dhst.v4.1.fout')[:,[1,6,7]]  # redshift, stellar mass, star formation rate

def find_ind(inp_list, inp_ind, find_list):
    diff = np.sqrt(np.sum((find_list-inp_list[inp_ind])**2,axis=1)) * 3600 # Their difference in arcsec
    out_ind = np.where(diff==diff.min())[0][0]
    if diff.min() < 0.15:
        return out_ind, diff.min()
    else:
        return None, diff.min()
#print find_ind(galfit_loc, 0, stellar_loc)
lists = []
gal_list = range(len(galfit_loc))
for i in gal_list:
    find_dex_0, min_diff = find_ind(galfit_loc, i, data[:,1:3])
    find_dex_1, min_diff = find_ind(galfit_loc, i, stellar_loc)
    if find_dex_0 is not None and find_dex_1 is not None:
        lists.append([i, find_dex_0, find_dex_1])        
    
#%%  Coadd the new catalog:  
results= []  #Re, Stellar
for i in range(len(lists)):
    result_i=[]
    result_i.append(data[lists[i][1]][4])      #0: B_T_m
    result_i.append(data[lists[i][1]][5])      #1: MAG_F125w
    result_i.append(data[lists[i][1]][6])      #2: RE_F125w(pixel 0.06)
    result_i.append(data[lists[i][1]][7])      #3: single_Sersic_n_F125w
    result_i.append(data[lists[i][1]][8])      #4: LogMass
    result_i.append(galfit[lists[i][0]][0])      #5: galfit mag
    result_i.append(galfit[lists[i][0]][1])      #6: galfit re
    result_i.append(galfit[lists[i][0]][2])      #7: galfit n
    result_i.append(stellar[lists[i][2]][0])      #8: redshift
    result_i.append(stellar[lists[i][2]][1])      #9: stellar mass
    result_i.append(stellar[lists[i][2]][2])      #9: star formation rate
    results.append(result_i)  #Redshift, Reff(arcsec), Sersic_n, Stellar_mass, flag, SFR
results = np.asarray(results)    
results_all = copy.deepcopy(results)
        
#%%
results = results[results[:, 7]>0]   #Serisc n >0
results = results[results[:, 10]>-17]   #Serisc n >0

results = results[(abs(results[:, 7]-results[:, 3])<1)]  # #Stellar mass in right range
results = results[[results[:, 4]>9][0] * [results[:, 4]<11.5][0]]  #Stellar mass in right range
results = results[(abs(results[:, 9]-results[:, 4])<1)]  # #Stellar mass in right range

#%%
#Compare the fitted n difference
plt.figure(figsize=(11, 11))
plt.scatter(results[:, 3],results[:, 7],
            c='darkred',s=280,marker=".",zorder=100, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
plt.title('Sersic n comparison')
plt.xlabel("Dimauro et. al.")
plt.ylabel("van der Wel et. al.")
plt.close()

#Compare the fitted Reff difference
#plt.figure(figsize=(11, 11))
#plt.scatter(results[:, 2]*0.06,results[:, 6],
#            c='darkred',s=280,marker=".",zorder=100, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
#plt.title('Reff comparsion')
#plt.xlabel("Dimauro et. al.")
#plt.ylabel("van der Wel et. al.")
#plt.xlim(0,10)
#plt.ylim(0,10)
#plt.show()

##Compare the fitted Stellar mass difference
#plt.figure(figsize=(11, 11))
#plt.scatter(results[:, 4],results[:, 9],
#            c='darkred',s=280,marker=".",zorder=100, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
#plt.title('Stellar mass comparison')
#plt.xlabel("Dimauro et. al.")
#plt.ylabel("3D HST")
#plt.xlim(8.5,11.7)
#plt.ylim(8.5,11.7)
#plt.show()

#%%
#Plot Sersic_n, B/T relation
#results = results[[results[:, 8]>1.2][0] * [results[:, 8]<1.8][0]]

plt.figure(figsize=(11, 11))
plt.scatter((results[:, 3]+results[:, 7])/2.,results[:, 0],
            c='darkred',s=280,marker=".",zorder=100, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
#plt.title('', fontsize=27)
plt.xlabel("Sersic_n",fontsize=27)
plt.ylabel("B/T mass ratio",fontsize=27)
plt.tick_params(labelsize=20)
plt.close()

#%%
from scipy.interpolate import spline
n_list = np.logspace(np.log10(0.3),np.log10(7),41)
sersic_n_arr = (results[:, 3]+results[:, 7])/2  #Average value of Sersic_n
BT_mean_list, BT_median_list = [], []
for n in n_list:
    if n <4:
        step = 0.1
    else:
        step = 0.5
    idx = [[sersic_n_arr>n-step][0] * [sersic_n_arr<n+step][0]][0]
    BT_mean_list.append(np.mean(results[:, 0][idx]))
    BT_median_list.append(np.median(results[:, 0][idx]))
    
plt.figure(figsize=(11, 11))
plt.scatter((results[:, 3]+results[:, 7])/2.,results[:, 0],
            c='gray',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
n_line = np.linspace(n_list.min(),n_list.max(),300) #300 represents number of points to make between T.min and T.max
BTR_smooth_mean = spline(n_list,BT_mean_list,n_line)
BTR_smooth_median = spline(n_list,BT_median_list,n_line)

#plt.plot(n_list, BT_mean_list,zorder=1,linewidth=2.0)
plt.plot(n_line,BTR_smooth_mean,'r',zorder=1,linewidth=3.0,label = 'The averaged B/T relation')
#plt.plot(n_line,BTR_smooth_median,'b',zorder=1,linewidth=3.0,label = 'By median value')
#plt.title('', fontsize=27)
plt.xlabel("Sersic index",fontsize=27)
plt.ylabel("B/T stellar mass ratio",fontsize=27)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':28})
plt.savefig('BT_relation.pdf')
plt.show()

#import pickle
#pickle.dump([n_line,BTR_smooth_mean,BTR_smooth_median], open("n_BT_relation.pkl", 'wb'))

#import pickle
#pickle.dump([sersic_n_arr,results[:, 0]], open("Sersic_BT_data.pkl", 'wb'))

'''
#%%Plot B/T , n relation with Reff as color
from scipy.integrate import quad
h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)
da_result = 1/(1+results[:,8])*c*vec_EE(results[:,8])/h0  #in Mpc
Reff_kpc = da_result * 10 **3 * (results[:,6]/3600./180.*np.pi)

plt.figure(figsize=(13, 11))
plt.scatter((results[:, 3]+results[:, 7])/2.,results[:, 0],
            c=np.log10(Reff_kpc),s=280,marker=".",zorder=0, vmin=-0.5, vmax=1.5, edgecolors='white',alpha=0.7)

#plt.title('', fontsize=27)

cl=plt.colorbar()  
cl.set_label('log10(Reff) (kpc)',rotation=270,size=30)
plt.xlabel("Sersic_n",fontsize=27)
plt.ylabel("B/T mass ratio",fontsize=27)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':28})
plt.show()

#%%
plt.figure(figsize=(11, 11))
plt.scatter((results[:, 9]+results[:, 4])/2.,results[:, 0],
            c='green',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
#plt.title('', fontsize=27)
plt.xlabel("Stellar mass",fontsize=27)
plt.ylabel("B/T mass ratio",fontsize=27)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':28})
plt.show()
#%%
plt.figure(figsize=(11, 11))
plt.scatter(results[:, 10],results[:, 0],
            c='green',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
#plt.title('', fontsize=27)
plt.xlabel("sSFR (log[ssfr*yr])",fontsize=27)
plt.ylabel("B/T mass ratio",fontsize=27)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':28})
plt.show()
##%%
#
#plt.figure(figsize=(11, 11))
#plt.scatter((results[:, 9]+results[:, 4])/2.,results[:, 0],
#            c='green',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
##plt.title('', fontsize=27)
#plt.xlabel("sSFR",fontsize=27)
#plt.ylabel("B/T mass ratio",fontsize=27)
#plt.tick_params(labelsize=20)
#plt.legend(prop={'size':28})
#plt.show()
'''