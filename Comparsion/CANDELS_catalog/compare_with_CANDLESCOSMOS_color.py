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
    result6 = stellar_flux_ap[lists[i][1]][0] #F140w
    result7 = stellar_flux_ap[lists[i][1]][1] #F125w.
    result8 = stellar_flux_ap[lists[i][1]][2] #F814w
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

relation = 0  # 0 M*- Reff ; 1 M*-color; 2 M*-n; 3 Reff-n relation
z_range = [1.2,1.7]

z_cut = ((results[:,0]>z_range[0]) * (results[:,0]<z_range[1]))

#ssfr_break =-10.5
ssfr_break =100


def logR_mstar(mstar, logA, alpha):
    """
    The R_Mstar relation in
    http://adsabs.harvard.edu/abs/2014ApJ...788...28V
    mstar in unit of log(M*/Msun)
    """
#    mstar = np.log10(10**mstar/7.*10**10.)
#    logr = logA + alpha*(mstar) + 10*alpha + alpha*np.log10(1/7.)
    r = 10**logA*(10**mstar/(5.*10**10))**alpha
    logr = np.log10(r)
    return logr

#blue_galaxy = ([results[:,5]>ssfr_break])[0]
all_galaxy = ([results[:,5]<ssfr_break])[0]  #As all galaxy
#z_cut = ([results[:,0]>0]) 

cmap_r = matplotlib.cm.get_cmap('jet')

plt.figure(figsize=(15, 11))
Reff_kpc = da_result * 10 **3 * (results[:,1]/3600./180.*np.pi)
mstar_cut = [(results[:,3]>9.5) * (results[:,3]<11.5)][0]
if relation == 0:
    plt.scatter(results[:,3][z_cut* all_galaxy],np.log10(Reff_kpc[z_cut * all_galaxy]),
                c=(results[:,7][z_cut* all_galaxy]/results[:,8][z_cut * all_galaxy]),s=280,marker=".",zorder=90, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)
    mstar_line = np.linspace(10.5,11.5,20)
    plt.plot(mstar_line, logR_mstar(mstar_line,logA=0.155 , alpha=0.76), 'r')
    mstar_line = np.linspace(9,11.5,20)
    plt.plot(mstar_line, logR_mstar(mstar_line,logA=0.675 , alpha=0.23), 'b')
elif relation == 1:
    plt.scatter(results[:,3][z_cut* all_galaxy],(results[:,7][z_cut* all_galaxy]/results[:,8][z_cut * all_galaxy]),
                c=(results[:,7][z_cut* all_galaxy]/results[:,8][z_cut * all_galaxy]),s=280,marker=".",zorder=90, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)
elif relation == 2:
    plt.scatter(results[:,3][z_cut* all_galaxy],results[:,2][z_cut* all_galaxy],
                c=(results[:,7][z_cut* all_galaxy]/results[:,8][z_cut * all_galaxy]),s=280,marker=".",zorder=90, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)
elif relation == 3:
    plt.scatter(np.log10(Reff_kpc[z_cut * all_galaxy * mstar_cut]),results[:,2][z_cut* all_galaxy* mstar_cut],
                c=(results[:,7][z_cut* all_galaxy* mstar_cut]/results[:,8][z_cut * all_galaxy* mstar_cut]),s=280,marker=".",zorder=90, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)
elif relation == 4:
    plt.scatter(results[:,0][z_cut* all_galaxy* mstar_cut], np.log10(Reff_kpc[z_cut * all_galaxy * mstar_cut]),
                c=(results[:,7][z_cut* all_galaxy* mstar_cut]/results[:,8][z_cut * all_galaxy* mstar_cut]),s=280,marker=".",zorder=90, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)
   
    z_line = np.linspace(1,2,21)
    plt.plot(z_line, np.log10(8.9*(1+z_line)**(-0.75)), 'b')
    plt.plot(z_line, np.log10(5.6*(1+z_line)**(-1.48)), 'r')
    

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
Reffs_e = np.array(load_re(ID, folder = '../../'))[:,1]
indexs = np.array(load_n(ID, folder = '../../'))[:,0]

dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #in pc
da=1/(1+zs)*c*vec_EE(zs)/h0   #in Mpc
ID_Reff_kpc = da * 10 **3 * (Reffs/3600./180.*np.pi)
ID_Reff_kpc_e = da * 10 **3 * ((Reffs_e)/3600./180.*np.pi)

from load_result import load_host_p

#from dmag import k_corr_R
#from filter_info import filt_info
#dm_k_R = []
#for i in range(len(zs)):
#    dm_k_R.append(k_corr_R(zs[i],filt_info[ID[i]], galaxy_age = '1Gyrs'))
#dm_k_R = np.asarray(dm_k_R) # Get the k-correction for each target as an array
#
#host_Mags = mags -5*(np.log10(dl)-1) + dm_k_R # This is in AB system
#host_Mags = host_Mags - 0.21  # Transfer to Vega system
#host_LR = 10 ** (0.4*(4.43-host_Mags))
#Mstar = np.log10(host_LR * 0.54 * 0.684 * 1.4191)  
Mstar = load_host_p(ID, folder = '../../')[1]
host_flux_WFC3 = np.array(load_flux(ID, folder = '../../', flt = 'WFC3'))[:,0]
host_flux_ACS = []
for i in range(len(ID)):
    ifexit = glob.glob('../../analysis_ACS/{0}'.format(ID[i]))
    if ifexit!= []:
        host_flux_ACS.append(load_flux([ID[i]], folder = '../../', flt = 'ACS')[0][0])
    else:
        host_flux_ACS.append(-99)
host_flux_ACS = np.asarray(host_flux_ACS)

cl=plt.colorbar()          #cl take the inforamtion from the lastest plt
cl.set_label('filter flux ratio, WFC3 / ACS',rotation=270,size=30)
cl.ax.get_yaxis().labelpad=35     #the distance of the colorbar titel from bar
cl.ax.tick_params(labelsize=30)


if relation == 0:
    f1 ='../../M_BH_relation/data//Bennert+2011_local.txt'
    b11_l = np.loadtxt(f1)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
    b11_local_Reff = b11_l[:,-1]
    b11_local_mstar = b11_l[:,4]
    plt.scatter(b11_local_mstar,np.log10(b11_local_Reff),s=180, c ='black',
                marker="o",zorder=100, vmin=0, vmax=7, edgecolors='white')     
    plt.scatter(Mstar[host_flux_ACS>0],np.log10(ID_Reff_kpc)[host_flux_ACS>0],s=680, c =(host_flux_WFC3/host_flux_ACS)[host_flux_ACS>0],
                marker="*",zorder=100, vmin=0, vmax=7, edgecolors='white')    
    plt.scatter(Mstar[host_flux_ACS>0],np.log10(ID_Reff_kpc)[host_flux_ACS>0],s=680, c =(host_flux_WFC3/host_flux_ACS)[host_flux_ACS>0],
                marker="*",zorder=100, vmin=0, vmax=7, edgecolors='white')
    plt.errorbar(Mstar[host_flux_ACS>0],np.log10(ID_Reff_kpc)[host_flux_ACS>0],
             yerr= (np.log10(ID_Reff_kpc)-np.log10(ID_Reff_kpc-ID_Reff_kpc_e))[host_flux_ACS>0],
             color='k',ecolor='k', fmt='.',markersize=1, zorder = 100)  
    plt.scatter(Mstar[host_flux_ACS<0],np.log10(ID_Reff_kpc)[host_flux_ACS<0],s=680, c ='none',
                marker="*",zorder=101, vmin=0, vmax=7, edgecolors='black')
    plt.errorbar(Mstar[host_flux_ACS<0],np.log10(ID_Reff_kpc)[host_flux_ACS<0],
             yerr= (np.log10(ID_Reff_kpc)-np.log10(ID_Reff_kpc-ID_Reff_kpc_e))[host_flux_ACS<0],
             color='k',ecolor='k', fmt='.',markersize=1, zorder = 100)  
elif relation ==1:
    plt.scatter(Mstar[host_flux_ACS>0],(host_flux_WFC3/host_flux_ACS)[host_flux_ACS>0],s=680, c =(host_flux_WFC3/host_flux_ACS)[host_flux_ACS>0],
                marker="*",zorder=100, vmin=0, vmax=7, edgecolors='white')
    
elif relation ==2:
    plt.scatter(Mstar[host_flux_ACS>0],indexs[host_flux_ACS>0],s=680, c =(host_flux_WFC3/host_flux_ACS)[host_flux_ACS>0],
                marker="*",zorder=100, vmin=0, vmax=7, edgecolors='white')
    plt.scatter(Mstar[host_flux_ACS<0],indexs[host_flux_ACS<0],s=180, c ='black',
                marker="s",zorder=99, vmin=0, vmax=7, edgecolors='white') 
    
elif relation ==3:
    plt.scatter(np.log10(ID_Reff_kpc)[host_flux_ACS>0],indexs[host_flux_ACS>0],s=680, c =(host_flux_WFC3/host_flux_ACS)[host_flux_ACS>0],
                marker="*",zorder=100, vmin=0, vmax=7, edgecolors='white')
    plt.errorbar(np.log10(ID_Reff_kpc)[host_flux_ACS>0],indexs[host_flux_ACS>0],
             xerr= (np.log10(ID_Reff_kpc)-np.log10(ID_Reff_kpc-ID_Reff_kpc_e))[host_flux_ACS>0],
             color='k',ecolor='orange', fmt='.',markersize=1, zorder = 10)    
    plt.scatter(np.log10(ID_Reff_kpc)[host_flux_ACS<0],indexs[host_flux_ACS<0],s=180, c ='black',
                marker="s",zorder=99, vmin=0, vmax=7, edgecolors='white') 
    plt.errorbar(np.log10(ID_Reff_kpc)[host_flux_ACS<0],indexs[host_flux_ACS<0],
             xerr= (np.log10(ID_Reff_kpc)-np.log10(ID_Reff_kpc-ID_Reff_kpc_e))[host_flux_ACS<0],
             color='k',ecolor='orange', fmt='.',markersize=1, zorder = 10)  
elif relation ==4:
    plt.scatter(zs[host_flux_ACS>0], np.log10(ID_Reff_kpc)[host_flux_ACS>0],s=680, c =(host_flux_WFC3/host_flux_ACS)[host_flux_ACS>0],
                marker="*",zorder=100, vmin=0, vmax=7, edgecolors='white')
    plt.errorbar(zs[host_flux_ACS>0], np.log10(ID_Reff_kpc)[host_flux_ACS>0],
             yerr= (np.log10(ID_Reff_kpc)-np.log10(ID_Reff_kpc-ID_Reff_kpc_e))[host_flux_ACS>0],
             color='k',ecolor='orange', fmt='.',markersize=1, zorder = 10)    
    plt.scatter(zs[host_flux_ACS<0], np.log10(ID_Reff_kpc)[host_flux_ACS<0],s=180, c ='black',
                marker="s",zorder=99, vmin=0, vmax=7, edgecolors='white') 
    plt.errorbar(zs[host_flux_ACS<0], np.log10(ID_Reff_kpc)[host_flux_ACS<0],
             yerr= (np.log10(ID_Reff_kpc)-np.log10(ID_Reff_kpc-ID_Reff_kpc_e))[host_flux_ACS<0],
             color='k',ecolor='orange', fmt='.',markersize=1, zorder = 10)  
   

plt.xlim([9.5, 12.0])
plt.xlabel("log$(M_*/M_{\odot})$",fontsize=35)
plt.tick_params(labelsize=25)
if relation ==0:
    plt.ylabel("log$(Reff)$ (kpc)",fontsize=35)
    plt.title('$M_* - Reff$ relation, sample redshift range {0}'.format(z_range), fontsize = 25)
    plt.ylim([-0.5, 1.5])
#    plt.savefig('Mstar-Reff_z{0}-{1}_color.pdf'.format(z_range[0],z_range[1]))
elif relation ==1:
    plt.ylabel("filter flux ratio, WFC3 / ACS",fontsize=35)
    plt.title('$M_* -$ color relation, sample redshift range {0}'.format(z_range), fontsize = 25)
    plt.ylim([0, 20])
#    plt.savefig('Mstar-color{0}-{1}_Reff.pdf'.format(z_range[0],z_range[1]))
elif relation ==2:
    plt.ylabel("Sersic index",fontsize=35)
    plt.title('$M_* - $Sersic index relation, sample redshift range {0}'.format(z_range), fontsize = 25)
    plt.ylim([0, 7])
#    plt.savefig('Mstar-Sersic_n{0}-{1}_color.pdf'.format(z_range[0],z_range[1]))
elif relation ==3:
    plt.xlabel("log$(Reff)$ (kpc)",fontsize=35)
    plt.ylabel("Sersic index",fontsize=35)
    plt.title('$Reff - $Sersic index relation, sample redshift range {0}'.format(z_range), fontsize = 25)
    plt.xlim([-0.5, 1.5])
    plt.ylim([0, 8])
#    plt.savefig('Reff-Sersic_n{0}-{1}_color.pdf'.format(z_range[0],z_range[1]))
elif relation ==4:
    plt.xlabel("z",fontsize=35)
    plt.ylabel("log$(Reff)$ (kpc)",fontsize=35)
    plt.title('$z- Reff$ relation, sample redshift range {0}'.format(z_range), fontsize = 25)
    plt.xlim([1., 2])
    plt.ylim([-0.5, 1.5])
#    plt.savefig('Reff-z_{0}-{1}_color.pdf'.format(z_range[0],z_range[1]))

#cl.ax.tick_params(labelsize=15)   #the labe size
plt.show()
