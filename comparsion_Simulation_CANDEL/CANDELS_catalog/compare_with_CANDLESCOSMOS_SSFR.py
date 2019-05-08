#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 15:46:19 2019

@author: Dartoon

Read the CANDLES-COSMOS catalog and see if the stellar mass and Reff has any relations:

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
cmap = matplotlib.cm.get_cmap('viridis')
import matplotlib as mpl
mpl.rc('image', cmap='jet')

from scipy.integrate import quad
h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)

#The COSMOS files in http://www.mpia.de/homes/vdwel/candels.html by van der Wel et al. (2012).
#   NUMBER         RA        DEC          f        mag       dmag         re        dre          n         dn          q         dq         pa        dpa          sn
galfit_loc = np.loadtxt('cos_2epoch_wfc3_f160w_060mas_v1.0_galfit.cat')[:,[1,2]]  #RA, DEC
galfit = np.loadtxt('cos_2epoch_wfc3_f160w_060mas_v1.0_galfit.cat')[:,[4,6,8,3]] # mag, re, n, flag

##The data from 3D HST
stellar_loc = np.loadtxt('cosmos_3dhst.v4.1.cat')[:,[3,4]]  # RA, DEC
stellar_flux_ap = np.loadtxt('cosmos_3dhst.v4.1.cat')[:,[69,51,39]]  # flux, F140w, F125w, F814w.
stellar = np.loadtxt('cosmos_3dhst.v4.1.fout')[:,[1,6,8]]  # redshift, stellar mass, star formation rate

## Check if they are in a same field.
#plt.plot(galfit_loc[:,0], galfit_loc[:,1],'.')
#plt.plot(stellar_loc[:,0], stellar_loc[:,1],'r.')
#plt.show()

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
    find_dex, min_diff = find_ind(galfit_loc, i, stellar_loc)
    if find_dex is not None:
        lists.append([i, find_dex])
        
results= []  #Re, Stellar
for i in range(len(lists)):
    result0 = stellar[lists[i][1]][0]     #Redshift
    result1 = galfit[lists[i][0]][1]      #Reff(arcsec)
    result2 = galfit[lists[i][0]][2]      #Sersic n
    result3 = stellar[lists[i][1]][1]     #Stellar_mass
    result4 = galfit[lists[i][0]][3]      #flag
    result5 = stellar[lists[i][1]][2]     #Star formation rate
    result6 = stellar_flux_ap[lists[i][1]][0]
    result7 = stellar_flux_ap[lists[i][1]][1]
    result8 = stellar_flux_ap[lists[i][1]][2]
    results.append([result0, result1, result2, result3, result4, result5, result6, result7, result8])  #Redshift, Reff(arcsec), Sersic_n, Stellar_mass, flag, SFR
results = np.asarray(results)

#Clean up the sample
results = results[results[:,0] != -1] # Flag as good
results = results[results[:,1] != -999.0]
results = results[results[:,3] != 0] # Flag as good
results = results[np.nan_to_num(results[:,5]) != -99]
results = results[np.nan_to_num(results[:,5]) != 0]
results = results[np.nan_to_num(results[:,3]) != -1]
results = results[np.nan_to_num(results[:,3]) != 0]

results = results[(results[:,8]) != -99]
results = results[(results[:,7]) != -99]
results = results[(results[:,8]) != 0]
results = results[(results[:,7]) != 0]

da_result = 1/(1+results[:,0])*c*vec_EE(results[:,0])/h0  #in Mpc

#%% Plot and input my sample

relation = 0  # 0 M*- Reff ; 1 M*-Sersic_n
z_range = [1.2,1.7]

z_cut = ((results[:,0]>z_range[0]) * (results[:,0]<z_range[1]))

ssfr_break =-10.5
blue_galaxy = ([results[:,5]>ssfr_break])[0]
red_galaxy = ([results[:,5]<ssfr_break])[0]
#z_cut = ([results[:,0]>0]) 

Reff_kpc = da_result * 10 **3 * (results[:,1]/3600./180.*np.pi)
plt.figure(figsize=(15, 11))
if relation == 0:
    plt.scatter(results[:,3][z_cut*red_galaxy],np.log10(Reff_kpc[z_cut*red_galaxy]),
                c='darkred',s=280,marker=".",zorder=100, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
    plt.scatter(results[:,3][z_cut*blue_galaxy],np.log10(Reff_kpc[z_cut*blue_galaxy]),
                c='navy',s=280,marker=".",zorder=99, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
if relation ==1:
    plt.scatter(results[:,3][z_cut*red_galaxy],results[:,2][z_cut*red_galaxy],
                c='darkred',s=280,marker=".",zorder=100, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
    plt.scatter(results[:,3][z_cut*blue_galaxy],results[:,2][z_cut*blue_galaxy],
                c='navy',s=280,marker=".",zorder=99, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)

cmap_r = matplotlib.cm.get_cmap('jet_r')
#    
#if relation == 0:
#    plt.scatter(results[:,3][z_cut],np.log10(Reff_kpc[z_cut]),
#                c=results[:,5][z_cut] ,s=280,marker=".",zorder=100, vmin=-16, vmax=-4, alpha=0.6, edgecolors='white', cmap=cmap_r)
#
#elif relation == 1:
#    plt.scatter(results[:,3][z_cut ],results[:,2][z_cut ],
#                c=results[:,5][z_cut],s=280,marker=".",zorder=100, vmin=-16, vmax=-4, alpha=0.6, edgecolors='white', cmap=cmap_r)
#    

##plt.clim(0,6)
#cl=plt.colorbar()          #cl take the inforamtion from the lastest plt
#
##plt.ylim([0, 5])
#cl.set_label('SSFR',rotation=270,size=20)

import sys
sys.path.insert(0,'../../py_tools')
from load_result import load_zs, load_mag, load_re, load_n

ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
zs = np.asarray(load_zs(ID))
mags = np.array(load_mag(ID, folder = '../../')[0])
Reffs = np.array(load_re(ID, folder = '../../'))[:,0]
indexs = np.array(load_n(ID, folder = '../../'))[:,0]

dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #in pc
da=1/(1+zs)*c*vec_EE(zs)/h0   #in Mpc
ID_Reff_kpc = da * 10 **3 * (Reffs/3600./180.*np.pi)


from dmag import k_corr_R
from filter_info import filt_info
dm_k_R = []
for i in range(len(zs)):
    dm_k_R.append(k_corr_R(zs[i],filt_info[ID[i]], galaxy_age = '1Gyrs'))
dm_k_R = np.asarray(dm_k_R) # Get the k-correction for each target as an array

host_Mags = mags -5*(np.log10(dl)-1) + dm_k_R # This is in AB system
host_Mags = host_Mags - 0.21  # Transfer to Vega system
host_LR = 10 ** (0.4*(4.61-host_Mags))
Mstar = np.log10(host_LR * 0.54 * 0.684 * 1.4191)  

if relation == 0:
    plt.scatter(Mstar,np.log10(ID_Reff_kpc),s=180, c =zs,marker="s",zorder=100, vmin=1.2, vmax=1.8, edgecolors='k')
elif relation ==1:
    plt.scatter(Mstar,indexs,s=180, c =zs,marker="s",zorder=100, vmin=1.2, vmax=1.8, edgecolors='k')
cl=plt.colorbar()          #cl take the inforamtion from the lastest plt
plt.clim(1.2,1.8)
plt.xlim([8, 12.5])
cl.set_label('Source redshift',rotation=270,size=30)
cl.ax.get_yaxis().labelpad=35     #the distance of the colorbar titel from bar
cl.ax.tick_params(labelsize=30)   
plt.xlabel("log$(M_*)/M_{\odot}$",fontsize=35)
plt.tick_params(labelsize=25)
if relation ==0:
    plt.ylabel("log$(Reff)$ (kpc)",fontsize=35)
    plt.title('$M_* - Reff$ relation, sample redshift range {0}'.format(z_range), fontsize = 25)
    plt.ylim([-1.5, 2.6])
    plt.savefig('Mstar-Reff_z{0}-{1}.pdf'.format(z_range[0],z_range[1]))
elif relation ==1:
    plt.ylabel("Sersic n",fontsize=35)
    plt.title('$M_* -$ Sersic_n relation, sample redshift range {0}'.format(z_range), fontsize = 25)
    plt.ylim([0, 8.2])
    plt.savefig('Mstar-Sn_z{0}-{1}.pdf'.format(z_range[0],z_range[1]))  
plt.show()

#%% Compare hist
#Compare their hist
from scipy import stats

QSO_stellar, QSO_reff, QSO_n = Mstar, np.log10(ID_Reff_kpc), indexs

Candles_stellar, Candles_reff, Candles_n =  results[:,3][z_cut], np.log10(Reff_kpc[z_cut]), results[:,2][z_cut]
candles_cut = ((Candles_stellar>QSO_stellar.min()) * (Candles_stellar<QSO_stellar.max()))
Candles_stellar, Candles_reff, Candles_n  = Candles_stellar[candles_cut], Candles_reff[candles_cut], Candles_n[candles_cut]

plt.figure(figsize=(10,6))
#common_params = dict(#bins=20, 
#                     normed=True,
#                     label=('QSO sample Reff','CANDLES-COSMOS sample Reff'))
#plt.hist([QSO_reff, Candles_reff], **common_params)
high0, x0, _ =plt.hist(QSO_reff, normed=True, histtype=u'step', label=('QSO sample'), linewidth = 2, color='orange')
high1, x1, _ =plt.hist(Candles_reff, normed=True, histtype=u'step', label=('CANDLES sample'), linewidth = 2, color='green')
x0_m = np.median(QSO_reff)
high_m0 = high0[np.where(abs(x0_m-x0) == abs(x0_m-x0).min())[0][0]]
x1_m = np.median(Candles_reff)
high_m1 = high1[np.where(abs(x1_m-x1) == abs(x1_m-x1).min())[0][0]]
plt.plot(np.linspace(0,high_m0)*0+np.median(QSO_reff) , np.linspace(0,high_m0), linewidth = 4,color='orange')
plt.plot(np.linspace(0,high_m1)*0+np.median(Candles_reff) , np.linspace(0,high_m1), linewidth = 4, color='green')
plt.text(np.median(QSO_reff)-0.2, high_m0*1.05, '{0}'.format(round(np.median(QSO_reff),3)), color='orange',fontsize=25)
plt.text(np.median(Candles_reff)-0.2, high_m1*1.05, '{0}'.format(round(np.median(Candles_reff),3)), color='green',fontsize=25)
plt.xlabel("Reff (kpc)",fontsize=27)
plt.ylabel("Density",fontsize=27)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.savefig('Hist_Reff.pdf'.format(z_range[0],z_range[1]))  
plt.show()
#Plot step
#fig = plt.figure(figsize=(8,6))
#ax = fig.add_subplot(1,1,1)
#value0 = np.sort(QSO_reff)
#ax.step(value0,[i/float(len(value0)-1) for i in range(len(value0))], label = 'QSO sample Reff')
#value1 = np.sort(Candles_reff)
#ax.step(value1,[i/float(len(value1)-1) for i in range(len(value1))], label = 'CANDLES-COSMOS sample Reff')
#ax.legend(prop={'size': 16})
#plt.tick_params(labelsize=15)
#plt.show() 
#pvalue = stats.ks_2samp(QSO_reff, Candles_reff).pvalue
#print "p-value between QSO_reff and Candles_reff:", round(pvalue,3)
#
#
plt.figure(figsize=(10,6))
#common_params = dict(#bins=20, 
#                     normed=True,
#                     label=('QSO sample Sersic_n','CANDLES-COSMOS sample Sersic_n'))
#plt.hist([QSO_n, Candles_n], **common_params)
#plt.legend()
#plt.show()
high0, x0, _ =plt.hist(QSO_n, normed=True, histtype=u'step', label=('QSO sample'), linewidth = 2, color='orange')
high1, x1, _ =plt.hist(Candles_n, normed=True, histtype=u'step', label=('CANDLES sample'), linewidth = 2, color='green')
x0_m = np.median(QSO_n)
high_m0 = high0[np.where(abs(x0_m-x0) == abs(x0_m-x0).min())[0][0]-1]
x1_m = np.median(Candles_n)
high_m1 = high1[np.where(abs(x1_m-x1) == abs(x1_m-x1).min())[0][0]-1]
plt.plot(np.linspace(0,high_m0)*0+np.median(QSO_n) , np.linspace(0,high_m0), linewidth = 4,color='orange')
plt.plot(np.linspace(0,high_m1)*0+np.median(Candles_n) , np.linspace(0,high_m1), linewidth = 4, color='green')
plt.text(np.median(QSO_n)-0.2, high_m0*1.05, '{0}'.format(round(np.median(QSO_n),3)), color='orange',fontsize=25)
plt.text(np.median(Candles_n)-0.2, high_m1*1.05, '{0}'.format(round(np.median(Candles_n),3)), color='green',fontsize=25)
plt.xlabel("Sersic index",fontsize=27)
plt.ylabel("Density",fontsize=27)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.savefig('Hist_Sn.pdf'.format(z_range[0],z_range[1]))  
plt.show()

##Plot step
#fig = plt.figure(figsize=(8,6))
#ax = fig.add_subplot(1,1,1)
#value0 = np.sort(QSO_n)
#ax.step(value0,[i/float(len(value0)-1) for i in range(len(value0))], label = 'QSO sample Sersic_n')
#value1 = np.sort(Candles_n)
#ax.step(value1,[i/float(len(value1)-1) for i in range(len(value1))], label = 'CANDLES-COSMOS sample Sersic_n')
#ax.legend(prop={'size': 16})
#plt.tick_params(labelsize=15)
#plt.show() 
#pvalue = stats.ks_2samp(QSO_n, Candles_n).pvalue
#print "p-value between QSO_n, Candles_n:", round(pvalue,3)


