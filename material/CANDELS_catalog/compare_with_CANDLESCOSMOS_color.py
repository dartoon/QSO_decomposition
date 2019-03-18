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
import glob

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

relation = 1  # 0 M*- Reff ; 1 M*-Sersic_n
z_range = [1,2]

z_cut = ((results[:,0]>z_range[0]) * (results[:,0]<z_range[1]))

#ssfr_break =-10.5
ssfr_break =100

#blue_galaxy = ([results[:,5]>ssfr_break])[0]
red_galaxy = ([results[:,5]<ssfr_break])[0]
#z_cut = ([results[:,0]>0]) 

cmap_r = matplotlib.cm.get_cmap('jet')

plt.figure(figsize=(15, 11))
if relation == 0:
    Reff_kpc = da_result * 10 **3 * (results[:,1]/3600./180.*np.pi)
    plt.scatter(results[:,3][z_cut* red_galaxy],np.log10(Reff_kpc[z_cut * red_galaxy]),
                c=(results[:,7][z_cut* red_galaxy]/results[:,8][z_cut * red_galaxy]),s=280,marker=".",zorder=100, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)

elif relation == 1:
    plt.scatter(results[:,3][z_cut* red_galaxy],(results[:,7][z_cut* red_galaxy]/results[:,8][z_cut * red_galaxy]),
                c=(results[:,7][z_cut* red_galaxy]/results[:,8][z_cut * red_galaxy]),s=280,marker=".",zorder=100, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)

import sys
sys.path.insert(0,'../../py_tools')
from load_result import load_zs, load_mag, load_re, load_n, load_flux

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

host_flux_WFC3 = np.array(load_flux(ID, folder = '../../', flt = 'WFC3'))[:,0]
host_flux_ACS = []
for i in range(len(ID)):
    ifexit = glob.glob('../../analysis_ACS/{0}'.format(ID[i]))
    if ifexit!= []:
        host_flux_ACS.append(load_flux([ID[i]], folder = '../../', flt = 'ACS')[0][0])
    else:
        host_flux_ACS.append(-99)
host_flux_ACS = np.asarray(host_flux_ACS)

if relation == 0:
    plt.scatter(Mstar[host_flux_ACS>0],np.log10(ID_Reff_kpc)[host_flux_ACS>0],s=180, c =(host_flux_WFC3/host_flux_ACS)[host_flux_ACS>0],
                marker="s",zorder=100, vmin=0, vmax=7, edgecolors='white')
elif relation ==1:
    plt.scatter(Mstar[host_flux_ACS>0],(host_flux_WFC3/host_flux_ACS)[host_flux_ACS>0],s=180, c =(host_flux_WFC3/host_flux_ACS)[host_flux_ACS>0],
                marker="s",zorder=100, vmin=0, vmax=7, edgecolors='white')
cl=plt.colorbar()          #cl take the inforamtion from the lastest plt
#plt.clim(0,7)

plt.xlim([8, 12.5])
plt.xlabel("log$(M_*)/M_{\odot}$",fontsize=35)
plt.tick_params(labelsize=25)
if relation ==0:
    plt.ylabel("log$(Reff)$ (kpc)",fontsize=35)
    plt.title('$M_* - Reff$ relation, CANDLES sample redshift range {0}'.format(z_range), fontsize = 25)
    plt.ylim([-1.5, 2.6])
    cl.set_label('Flux, filter WFC3 / ACS',rotation=270,size=20)
    plt.savefig('Mstar-Reff_z{0}-{1}_color.pdf'.format(z_range[0],z_range[1]))
elif relation ==1:
    plt.ylabel("Flux, filter WFC3 / ACS",fontsize=35)
    plt.title('$M_* -$ color relation, CANDLES sample redshift range {0}'.format(z_range), fontsize = 25)
    plt.ylim([0, 20])
    cl.set_label('Flux, filter WFC3 / ACS',rotation=270,size=20)
    plt.savefig('Mstar-color{0}-{1}_Reff.pdf'.format(z_range[0],z_range[1]))
    
cl.ax.tick_params(labelsize=15)   #the labe size
    
    
plt.xlim([8, 12.5])
plt.show()
