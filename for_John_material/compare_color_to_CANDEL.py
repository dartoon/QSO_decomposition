#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 13:33:24 2019

@author: Dartoon

Compare the color
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

import sys
sys.path.insert(0,'../py_tools')
from load_result import load_zs, load_flux, load_host_p

ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281', 'CID255']
zs = np.asarray(load_zs(ID))

host_flux_WFC3 = np.array(load_flux(ID, folder = '../', flt = 'WFC3'))[:,0]
host_flux_ACS = []
for i in range(len(ID)):
    ifexit = glob.glob('../analysis_ACS/{0}'.format(ID[i]))
    if ifexit!= []:
        host_flux_ACS.append(load_flux([ID[i]], folder = '../', flt = 'ACS')[0][0])
    else:
        host_flux_ACS.append(-99)
        
host_flux_ACS = np.asarray(host_flux_ACS)
Mstar = load_host_p(ID, folder = '../')[1]

Mstar = Mstar[host_flux_ACS!=-99]
zs = zs[host_flux_ACS!=-99]
host_flux_WFC3 = host_flux_WFC3[host_flux_ACS!=-99]
host_flux_ACS = host_flux_ACS[host_flux_ACS!=-99]

#%% Load CANDEL information:
#The COSMOS files in http://www.mpia.de/homes/vdwel/candels.html by van der Wel et al. (2012).
#   NUMBER         RA        DEC          f        mag       dmag         re        dre          n         dn          q         dq         pa        dpa          sn
file_galfit_COSMOS =  np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/cos_2epoch_wfc3_f160w_060mas_v1.0_galfit.cat')
galfit_COSMOS_loc = file_galfit_COSMOS[:,[1,2]]  #RA, DEC
galfit_COSMOS = file_galfit_COSMOS[:,[4,6,8,3]] # mag, re, n, flag
##The data from 3D HST: https://3dhst.research.yale.edu/Data.php  (PHOTOMETRY)
file_stellar_COSMOS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/cosmos_3dhst.v4.1.cat')
stellar_COSMOS_loc = file_stellar_COSMOS[:,[3,4]]  # RA, DEC
stellar_COSMOS_flux_ap = file_stellar_COSMOS[:,[69,51,39]]  # flux, F140w, F125w, F814w.
stellar_COSMOS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/cosmos_3dhst.v4.1.fout')[:,[1,6,8]]  # redshift, stellar mass, star formation rate
color_COSMOS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/cosmos_3dhst.v4.1.cats/RF_colors/cosmos_3dhst.v4.1.master.RF')[:,[3, 7, 9]] #Johnson_U, Johnson_V, 2MASS/J,
## Check if they are in a same field.
#plt.plot(galfit_COSMOS_loc[:,0], galfit_COSMOS_loc[:,1],'.')
#plt.plot(stellar_COSMOS_loc[:,0], stellar_COSMOS_loc[:,1],'r.')
#plt.show()
#galfit_loc = galfit_COSMOS_loc
#galfit = galfit_COSMOS
#stellar_loc = stellar_COSMOS_loc
#stellar_flux_ap = stellar_COSMOS_flux_ap
#stellar = stellar_COSMOS
#color = color_COSMOS

####Extended on 2019/10/26 with the other 4 fields: AEGIS, Good_N, Good_S, UDS.
file_galfit_AEGIS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/egs_2epoch_wfc3_f160w_060mas_v0.8_galfit.cat')
galfit_AEGIS_loc = file_galfit_AEGIS[:,[1,2]]  #RA, DEC
galfit_AEGIS = file_galfit_AEGIS[:,[4,6,8,3]] # mag, re, n, flag
file_stellar_AEGIS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/aegis_3dhst.v4.1.cats/Catalog/aegis_3dhst.v4.1.cat')
stellar_AEGIS_loc = file_stellar_AEGIS[:,[3,4]]  # RA, DEC
stellar_AEGIS_flux_ap = file_stellar_AEGIS[:,[48,33,27]]  # flux, (F140w, F125w, F814w. filter has checked using aegis.readme)
stellar_AEGIS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/aegis_3dhst.v4.1.cats/Fast/aegis_3dhst.v4.1.fout')[:,[1,6,8]]  # redshift, stellar mass, star formation rate
color_AEGIS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/aegis_3dhst.v4.1.cats/RF_colors/aegis_3dhst.v4.1.master.RF')[:,[3, 7, 9]] #Johnson_U, Johnson_V, 2MASS/J,
#plt.plot(galfit_AEGIS_loc[:,0], galfit_AEGIS_loc[:,1],'.')
#plt.plot(stellar_AEGIS_loc[:,0], stellar_AEGIS_loc[:,1],'r.')
#plt.show()
#
#file_galfit_gn = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/gn_all_candels_wfc3_f160w_060mas_v0.8_galfit.cat')
#galfit_gn_loc = file_galfit_gn[:,[1,2]]  #RA, DEC
#galfit_gn = file_galfit_gn[:,[4,6,8,3]] # mag, re, n, flag
#file_stellar_gn = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodsn_3dhst.v4.1.cats/Catalog/goodsn_3dhst.v4.1.cat')
#stellar_gn_loc = file_stellar_gn[:,[3,4]]  # RA, DEC
#stellar_gn_flux_ap = file_stellar_gn[:,[69,51,39]]  # flux, (F140w: 54, F125w:48, F814w, missing)
#stellar_gn = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodsn_3dhst.v4.1.cats/Fast/goodsn_3dhst.v4.1.fout')[:,[1,6,8]]  # redshift, stellar mass, star formation rate
#color_gn = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodsn_3dhst.v4.1.cats/RF_colors/goodsn_3dhst.v4.1.master.RF')[:,[3, 7, 9]] #Johnson_U, Johnson_V, 2MASS/J,

##plt.plot(galfit_gn_loc[:,0], galfit_gn_loc[:,1],'.')
##plt.plot(stellar_gn_loc[:,0], stellar_gn_loc[:,1],'r.')
##plt.show()
#
file_galfit_gs = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/gs_all_candels_ers_udf_f160w_v0.5_galfit.cat')
galfit_gs_loc = file_galfit_gs[:,[1,2]]  #RA, DEC
galfit_gs = file_galfit_gs[:,[4,6,8,3]] # mag, re, n, flag
file_stellar_gs = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat')
stellar_gs_loc = file_stellar_gs[:,[3,4]]  # RA, DEC
stellar_gs_flux_ap = file_stellar_gs[:,[63,54,45]]  # flux, (F140w, F125w, F814w. filter has checked using *.readme)
stellar_gs = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodss_3dhst.v4.1.cats/Fast/goodss_3dhst.v4.1.fout')[:,[1,6,8]]  # redshift, stellar mass, star formation rate
color_gs = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodss_3dhst.v4.1.cats/RF_colors/goodss_3dhst.v4.1.master.RF')[:,[3, 7, 9]] #Johnson_U, Johnson_V, 2MASS/J,
#plt.plot(galfit_gs_loc[:,0], galfit_gs_loc[:,1],'.')
#plt.plot(stellar_gs_loc[:,0], stellar_gs_loc[:,1],'r.')
#plt.show()
#
file_galfit_uds = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/uds_2epoch_wfc3_f160w_060mas_v0.3_galfit.cat')
galfit_uds_loc = file_galfit_uds[:,[1,2]]  #RA, DEC
galfit_uds = file_galfit_uds[:,[4,6,8,3]] # mag, re, n, flag
file_stellar_uds = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/uds_3dhst.v4.2.cats/Catalog/uds_3dhst.v4.2.cat')
stellar_uds_loc = file_stellar_uds[:,[3,4]]  # RA, DEC
stellar_uds_flux_ap = file_stellar_uds[:,[42,36,30]]  # flux, (F140w, F125w, F814w. filter has checked from *.cat file)
stellar_uds = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/uds_3dhst.v4.2.cats/Fast/uds_3dhst.v4.2.fout')[:,[1,6,8]]  # redshift, stellar mass, star formation rate
color_uds = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/uds_3dhst.v4.2.cats/RF_colors/uds_3dhst.v4.2.master.RF')[:,[3, 7, 9]] #Johnson_U, Johnson_V, 2MASS/J,
#plt.plot(galfit_uds_loc[:,0], galfit_uds_loc[:,1],'.')
#plt.plot(stellar_uds_loc[:,0], stellar_uds_loc[:,1],'r.')
#plt.show()
#
galfit_loc = np.concatenate((galfit_COSMOS_loc, galfit_AEGIS_loc, galfit_gs_loc, galfit_uds_loc))
galfit = np.concatenate((galfit_COSMOS, galfit_AEGIS, galfit_gs, galfit_uds))
stellar_loc = np.concatenate((stellar_COSMOS_loc, stellar_AEGIS_loc, stellar_gs_loc, stellar_uds_loc))
stellar_flux_ap = np.concatenate((stellar_COSMOS_flux_ap, stellar_AEGIS_flux_ap, stellar_gs_flux_ap,stellar_uds_flux_ap))
stellar = np.concatenate((stellar_COSMOS, stellar_AEGIS, stellar_gs,stellar_uds))
color = np.concatenate((color_COSMOS, color_AEGIS, color_gs,color_uds))

#%% Summary the data and Combin them
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
    result9 = color[lists[i][1]]
    results.append([result0, result1, result2, result3, result4, result5, result6, result7, result8, result9])  #Redshift, Reff(arcsec), Sersic_n, Stellar_mass, flag, SFR
results = np.asarray(results)

#Clean up the sample
results = results[results[:,0] != -1] # Flag as good
results = results[results[:,1] != -999.0]
results = results[results[:,3] != 0] # Flag as good
results = results[np.nan_to_num(results[:,5]) != -99]
results = results[np.nan_to_num(results[:,5]) != 0]
results = results[np.nan_to_num(results[:,3]) != -1]
results = results[np.nan_to_num(results[:,3]) != 0]

results = results[(results[:,6]) != -99]
results = results[(results[:,8]) != -99]
results = results[(results[:,7]) != -99]
results = results[(results[:,6]) != 0]
results = results[(results[:,7]) != 0]
results = results[(results[:,8]) != 0]
#
results = results[(results[:,6]) >0]
results = results[(results[:,7]) >0]
results = results[(results[:,8]) >0]

z_range = [1.2, 1.44, 1.7]
z_cut_0 = ((results[:,0]>z_range[0]) * (results[:,0]<z_range[1]))
z_cut_1 = ((results[:,0]>z_range[1]) * (results[:,0]<z_range[2]))

mstar_cut = [(results[:,3]>9.5) * (results[:,3]<11.5)][0]
#mstar_cut = [(results[:,3]>10.4) * (results[:,3]<10.9)][0]
Reff_cut = [(results[:,1]>0.2)][0]

#%%
cosmos6, cosmos7, cosmos8 = results[:,6], results[:,7], results[:,8]
cosmos6, cosmos7, cosmos8 = cosmos6.astype('float64'), cosmos7.astype('float64'), cosmos8.astype('float64')


# =============================================================================
# Flux not make good sense since the zeropoint is different.
# =============================================================================
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
#
#plt.figure(figsize=(8,6))
#high0, high_x_0, _ = plt.hist((host_flux_ACS/host_flux_WFC3)[zs<z_range[1]], normed=True, histtype=u'step',
#         label=('high-z AGN, z<1.44'), linewidth = 2, color='blue')
#high0, high_x_0, _ = plt.hist((host_flux_ACS/host_flux_WFC3)[zs>z_range[1]], normed=True, histtype=u'step',
#         label=('high-z AGN, z>1.44'), linewidth = 2, color='red')
#candel0, candel_x0, _ = plt.hist((results[:,8]/results[:,7])[z_cut_0*mstar_cut*Reff_cut], normed=True, histtype=u'step',
#         label=('COSMOS galaxy, z<1.44'), linewidth = 2, color='c')
#candel1, candel_x1, _ = plt.hist((results[:,8]/results[:,6])[z_cut_1*mstar_cut*Reff_cut], normed=True, histtype=u'step',
#         label=('COSMOS galaxy, z>1.44'), linewidth = 2, color='lightcoral')
#plt.xlabel("ACS/WFC3 flux ratio",fontsize=27)
#plt.ylabel("Density",fontsize=27)
#plt.legend(prop={'size':20})
#plt.tick_params(labelsize=20)
#plt.show()
#
#plt.figure(figsize=(8,6))
#high0, high_x_0, _ = plt.hist((host_flux_WFC3)[zs<z_range[1]], normed=True, histtype=u'step',
#         label=('high-z AGN, z<1.44'), linewidth = 2, color='blue')
#high0, high_x_0, _ = plt.hist((host_flux_WFC3)[zs>z_range[1]], normed=True, histtype=u'step',
#         label=('high-z AGN, z>1.44'), linewidth = 2, color='red')
#candel0, candel_x0, _ = plt.hist((results[:,7])[z_cut_0*mstar_cut*Reff_cut], normed=True, histtype=u'step',
#         label=('COSMOS galaxy, z<1.44'), linewidth = 2, color='c')
#candel1, candel_x1, _ = plt.hist((results[:,6])[z_cut_1*mstar_cut*Reff_cut], normed=True, histtype=u'step',
#         label=('COSMOS galaxy, z>1.44'), linewidth = 2, color='lightcoral')
#plt.xlabel("WFC3 flux",fontsize=27)
#plt.ylabel("Density",fontsize=27)
#plt.legend(prop={'size':20})
#plt.tick_params(labelsize=20)
#plt.show()
#
#plt.figure(figsize=(8,6))
#high0, high_x_0, _ = plt.hist((host_flux_ACS)[zs<z_range[1]], normed=True, histtype=u'step',
#         label=('high-z AGN, z<1.44'), linewidth = 2, color='blue')
#high0, high_x_0, _ = plt.hist((host_flux_ACS)[zs>z_range[1]], normed=True, histtype=u'step',
#         label=('high-z AGN, z>1.44'), linewidth = 2, color='red')
#candel0, candel_x0, _ = plt.hist((results[:,8])[z_cut_0*mstar_cut*Reff_cut], normed=True, histtype=u'step',
#         label=('COSMOS galaxy, z<1.44'), linewidth = 2, color='c')
#candel1, candel_x1, _ = plt.hist((results[:,8])[z_cut_1*mstar_cut*Reff_cut], normed=True, histtype=u'step',
#         label=('COSMOS galaxy, z>1.44'), linewidth = 2, color='lightcoral')
#plt.xlabel("ACS flux",fontsize=27)
#plt.ylabel("Density",fontsize=27)
#plt.legend(prop={'size':20})
#plt.tick_params(labelsize=20)
#plt.show()

#%%

#        if filt == "F140w":
#            zp = 26.4524
#        elif filt == "F125w":
#            zp = 26.2303
#        elif filt == "F814w":
#            zp = 25.94333
def cal_mag(flux, zp):
    return -2.5*np.log10(flux) + zp

mstar_cut = [(results[:,3]>9.5) * (results[:,3]<11.5)][0]
plt.figure(figsize=(12,9))
plt.scatter(Mstar[zs<z_range[1]], (cal_mag(host_flux_ACS, 25.94333)- cal_mag(host_flux_WFC3,  26.2303))[zs<z_range[1]], edgecolors='black',
         label=('high-z AGN, z<1.44'), color='blue', zorder = 10, marker="*", s=440)
plt.scatter(Mstar[zs>z_range[1]], (cal_mag(host_flux_ACS, 25.94333)- cal_mag(host_flux_WFC3,  26.4524))[zs>z_range[1]], edgecolors='black',
         label=('high-z AGN, z>1.44'), color='red', zorder = 10, marker="*", s=440)
plt.scatter(results[:,3][z_cut_0*mstar_cut*Reff_cut], (cal_mag(cosmos8, 25.) - cal_mag(cosmos7, 25.))[z_cut_0*mstar_cut*Reff_cut],
         label=('CANDELS galaxy, z<1.44'),  color='c')
plt.scatter(results[:,3][z_cut_1*mstar_cut*Reff_cut], (cal_mag(cosmos8, 25.) - cal_mag(cosmos6, 25.))[z_cut_1*mstar_cut*Reff_cut],
         label=('CANDELS galaxy, z>1.44'), color='lightcoral')
plt.xlabel("log (M$_*$; units of M$_{\odot}$)",fontsize=35)
plt.ylabel(r"$\Delta$mag (ACS - WFC3)",fontsize=27)
plt.legend(prop={'size':15})
plt.ylim([0, 4])
plt.tick_params(labelsize=20)
plt.show()

mstar_cut = [(results[:,3]>10.4) * (results[:,3]<10.9)][0]
plt.figure(figsize=(8,6))
high0, high_x_0, _ = plt.hist( (cal_mag(host_flux_ACS, 25.94333)- cal_mag(host_flux_WFC3,  26.2303))[zs<z_range[1]], normed=True, histtype=u'step',
         label=('high-z AGN, z<1.44'), linewidth = 2, color='blue')
high0, high_x_0, _ = plt.hist( (cal_mag(host_flux_ACS, 25.94333)- cal_mag(host_flux_WFC3,  26.4524))[zs>z_range[1]], normed=True, histtype=u'step',
         label=('high-z AGN, z>1.44'), linewidth = 2, color='red')
candel0, candel_x0, _ = plt.hist( (cal_mag(cosmos8, 25.) - cal_mag(cosmos7, 25.))[z_cut_0*mstar_cut*Reff_cut], normed=True, histtype=u'step', linestyle=('dashed'),
         label=('CANDELS galaxy, z<1.44'), linewidth = 2, color='c')
candel1, candel_x1, _ = plt.hist( (cal_mag(cosmos8, 25.) - cal_mag(cosmos6, 25.))[z_cut_1*mstar_cut*Reff_cut], normed=True, histtype=u'step', linestyle=('dashed'),
         label=('CANDELS galaxy, z>1.44'), linewidth = 2, color='lightcoral')
plt.xlabel(r"$\Delta$mag (ACS - WFC3)",fontsize=27)
plt.ylabel("Density",fontsize=27)
plt.legend(prop={'size':15})
plt.tick_params(labelsize=20)
plt.show()
from scipy import stats
print "p-value: high_z VS CANDELS, z<1.44:", round(stats.ks_2samp((cal_mag(host_flux_ACS, 25.94333)- cal_mag(host_flux_WFC3,  26.2303))[zs<z_range[1]],
                                                                 (cal_mag(cosmos8, 25.) - cal_mag(cosmos7, 25.))[z_cut_0*mstar_cut*Reff_cut]).pvalue,3)
print "p-value: high_z VS CANDELS, z>1.44:", round(stats.ks_2samp((cal_mag(host_flux_ACS, 25.94333)- cal_mag(host_flux_WFC3,  26.4524))[zs>z_range[1]],
                                                                 (cal_mag(cosmos8, 25.) - cal_mag(cosmos6, 25.))[z_cut_1*mstar_cut*Reff_cut]).pvalue,3)
