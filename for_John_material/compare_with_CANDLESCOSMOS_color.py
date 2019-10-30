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

####Extended on 2019/10/26 with the other 4 fields: AEGIS, Good_N, Good_S, UDS.
file_galfit_AEGIS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/egs_2epoch_wfc3_f160w_060mas_v0.8_galfit.cat')
galfit_AEGIS_loc = file_galfit_AEGIS[:,[1,2]]  #RA, DEC
galfit_AEGIS = file_galfit_AEGIS[:,[4,6,8,3]] # mag, re, n, flag
file_stellar_AEGIS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/aegis_3dhst.v4.1.cats/Catalog/aegis_3dhst.v4.1.cat')
stellar_AEGIS_loc = file_stellar_AEGIS[:,[3,4]]  # RA, DEC
stellar_AEGIS_flux_ap = file_stellar_AEGIS[:,[69,51,39]]  # flux, (F140w, F125w, F814w. filter is not checked, copy from COSMOS line)
stellar_AEGIS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/aegis_3dhst.v4.1.cats/Fast/aegis_3dhst.v4.1.fout')[:,[1,6,8]]  # redshift, stellar mass, star formation rate
color_AEGIS = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/aegis_3dhst.v4.1.cats/RF_colors/aegis_3dhst.v4.1.master.RF')[:,[3, 7, 9]] #Johnson_U, Johnson_V, 2MASS/J,
#plt.plot(galfit_AEGIS_loc[:,0], galfit_AEGIS_loc[:,1],'.')
#plt.plot(stellar_AEGIS_loc[:,0], stellar_AEGIS_loc[:,1],'r.')
#plt.show()

file_galfit_gn = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/gn_all_candels_wfc3_f160w_060mas_v0.8_galfit.cat')
galfit_gn_loc = file_galfit_gn[:,[1,2]]  #RA, DEC
galfit_gn = file_galfit_gn[:,[4,6,8,3]] # mag, re, n, flag
file_stellar_gn = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodsn_3dhst.v4.1.cats/Catalog/goodsn_3dhst.v4.1.cat')
stellar_gn_loc = file_stellar_gn[:,[3,4]]  # RA, DEC
stellar_gn_flux_ap = file_stellar_gn[:,[69,51,39]]  # flux, (F140w, F125w, F814w. filter is not checked, copy from COSMOS line)
stellar_gn = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodsn_3dhst.v4.1.cats/Fast/goodsn_3dhst.v4.1.fout')[:,[1,6,8]]  # redshift, stellar mass, star formation rate
color_gn = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodsn_3dhst.v4.1.cats/RF_colors/goodsn_3dhst.v4.1.master.RF')[:,[3, 7, 9]] #Johnson_U, Johnson_V, 2MASS/J,

#plt.plot(galfit_gn_loc[:,0], galfit_gn_loc[:,1],'.')
#plt.plot(stellar_gn_loc[:,0], stellar_gn_loc[:,1],'r.')
#plt.show()

file_galfit_gs = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/gs_all_candels_ers_udf_f160w_v0.5_galfit.cat')
galfit_gs_loc = file_galfit_gs[:,[1,2]]  #RA, DEC
galfit_gs = file_galfit_gs[:,[4,6,8,3]] # mag, re, n, flag
file_stellar_gs = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat')
stellar_gs_loc = file_stellar_gs[:,[3,4]]  # RA, DEC
stellar_gs_flux_ap = file_stellar_gs[:,[69,51,39]]  # flux, (F140w, F125w, F814w. filter is not checked, copy from COSMOS line)
stellar_gs = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodss_3dhst.v4.1.cats/Fast/goodss_3dhst.v4.1.fout')[:,[1,6,8]]  # redshift, stellar mass, star formation rate
color_gs = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/goodss_3dhst.v4.1.cats/RF_colors/goodss_3dhst.v4.1.master.RF')[:,[3, 7, 9]] #Johnson_U, Johnson_V, 2MASS/J,
#plt.plot(galfit_gs_loc[:,0], galfit_gs_loc[:,1],'.')
#plt.plot(stellar_gs_loc[:,0], stellar_gs_loc[:,1],'r.')
#plt.show()

file_galfit_uds = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/uds_2epoch_wfc3_f160w_060mas_v0.3_galfit.cat')
galfit_uds_loc = file_galfit_uds[:,[1,2]]  #RA, DEC
galfit_uds = file_galfit_uds[:,[4,6,8,3]] # mag, re, n, flag
file_stellar_uds = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/uds_3dhst.v4.2.cats/Catalog/uds_3dhst.v4.2.cat')
stellar_uds_loc = file_stellar_uds[:,[3,4]]  # RA, DEC
stellar_uds_flux_ap = file_stellar_uds[:,[69,51,39]]  # flux, (F140w, F125w, F814w. filter is not checked, copy from COSMOS line)
stellar_uds = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/uds_3dhst.v4.2.cats/Fast/uds_3dhst.v4.2.fout')[:,[1,6,8]]  # redshift, stellar mass, star formation rate
color_uds = np.loadtxt('../Comparsion/CANDELS_catalog/CANDELS_data/uds_3dhst.v4.2.cats/RF_colors/uds_3dhst.v4.2.master.RF')[:,[3, 7, 9]] #Johnson_U, Johnson_V, 2MASS/J,
#plt.plot(galfit_uds_loc[:,0], galfit_uds_loc[:,1],'.')
#plt.plot(stellar_uds_loc[:,0], stellar_uds_loc[:,1],'r.')
#plt.show()

galfit_loc = np.concatenate((galfit_COSMOS_loc, galfit_AEGIS_loc, galfit_gn_loc, galfit_gs_loc, galfit_uds_loc))
galfit = np.concatenate((galfit_COSMOS, galfit_AEGIS, galfit_gn, galfit_gs, galfit_uds))
stellar_loc = np.concatenate((stellar_COSMOS_loc, stellar_AEGIS_loc, stellar_gn_loc, stellar_gs_loc, stellar_uds_loc))
stellar_flux_ap = np.concatenate((stellar_COSMOS_flux_ap, stellar_AEGIS_flux_ap, stellar_gn_flux_ap, stellar_gs_flux_ap,stellar_uds_flux_ap))
stellar = np.concatenate((stellar_COSMOS, stellar_AEGIS, stellar_gn, stellar_gs,stellar_uds))
color = np.concatenate((color_COSMOS, color_AEGIS, color_gn, color_gs,color_uds))

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

results = results[(results[:,8]) != -99]
results = results[(results[:,7]) != -99]
results = results[(results[:,8]) != 0]
results = results[(results[:,7]) != 0]

da_result = 1/(1+results[:,0])*c*vec_EE(results[:,0])/h0  #in Mpc

#%% Plot and input my sample

relation = 0  # 0 M*- Reff ; 1 M*-color; 2 M*-n; 3 Reff-n relation
z_range = [1.2,1.7]
z_cut = ((results[:,0]>z_range[0]) * (results[:,0]<z_range[1]))

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


#ssfr_break =-10.5
#blue_galaxy = ([results[:,5]>ssfr_break])[0]
#red_galaxy = ([results[:,5]<ssfr_break])[0]
all_galaxy = ([results[:,5]<100])[0]  #As all galaxy

#Define based on Color:
red_galaxy = []
blue_galaxy = []

all_color = []
for i in range(len(results)):
    all_color.append([results[:,9][i][0], results[:,9][i][1], results[:,9][i][2]])
all_color = np.asarray(all_color)
ColorUV = -(2.5* np.log10(all_color[:,0])-2.5* np.log10(all_color[:,1]))
ColorVJ = -(2.5* np.log10(all_color[:,1])-2.5* np.log10(all_color[:,2]))
for i in range(len(all_color)):
    if ColorVJ[i] < 0.8425:
        blue_galaxy.append((ColorUV[i] < 1.286))
        red_galaxy.append((ColorUV[i] > 1.286))
    else:
        k = 1.17
        b = 0.3
        line_p = k*ColorVJ[i]+b
        blue_galaxy.append((ColorUV[i] -line_p < 0))        
        red_galaxy.append((ColorUV[i] -line_p > 0))
blue_galaxy = np.asarray(blue_galaxy)
red_galaxy = np.asarray(red_galaxy)
#%%
#import scipy.optimize as opt
#def lfit(x,m,c):
#    return m*x+c

cmap_r = matplotlib.cm.get_cmap('RdBu_r')

fig, ax = plt.subplots(figsize=(14, 11))
Reff_kpc = da_result * 10 **3 * (results[:,1]/3600./180.*np.pi)
Reff_kpc = Reff_kpc.astype(float)

Mstar_candels = results[:,3]
Mstar_candels = Mstar_candels.astype(float)

mstar_cut = [(results[:,3]>9.5) * (results[:,3]<11.5)][0]
if relation == 0:
#    plt.scatter(results[:,3][z_cut* all_galaxy],np.log10(Reff_kpc[z_cut * all_galaxy]),
#                c=(results[:,7][z_cut* all_galaxy]/results[:,8][z_cut * all_galaxy]),s=280,marker=".",zorder=90, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)
#    mstar_line = np.linspace(10.5,11.5,20)
#    plt.plot(mstar_line, logR_mstar(mstar_line,logA=0.155 , alpha=0.76), 'r')
#    mstar_line = np.linspace(9,11.5,20)
#    plt.plot(mstar_line, logR_mstar(mstar_line,logA=0.675 , alpha=0.23), 'b')
    plt.scatter(Mstar_candels[z_cut* blue_galaxy],Reff_kpc[z_cut * blue_galaxy],
                c='lightskyblue',s=280,marker=".",zorder=-90, alpha=0.6, edgecolors='white', cmap=cmap_r, label='CANDELS galaxy, star-forming')
    plt.scatter(Mstar_candels[z_cut* red_galaxy],Reff_kpc[z_cut * red_galaxy],
                c='darksalmon',s=280,marker=".",zorder=-90, alpha=0.6, edgecolors='white', cmap=cmap_r, label='CANDELS galaxy, quiescent')
    
    mstar_line = np.linspace(10.5,11.5,20)
#    m_cut = [(results[:,3]>10) * (results[:,3]<11.5)][0]
#    fit_red = opt.curve_fit(lfit, results[:,3][z_cut* red_galaxy * m_cut],np.log10(Reff_kpc[z_cut * red_galaxy * m_cut]))
#    plt.plot(mstar_line, 10**(lfit(mstar_line, fit_red[0][0],fit_red[0][1])), 'r',linewidth=3)
    plt.plot(mstar_line, 10**(logR_mstar(mstar_line,logA=0.155 , alpha=0.76)), 'r',linewidth=3)
    
    mstar_line = np.linspace(9,11.5,20)
    plt.plot(mstar_line, 10**(logR_mstar(mstar_line,logA=0.675 , alpha=0.23)), 'b',linewidth=3)   
    
    mstar_line = np.linspace(9,11.5,20)
    plt.plot(mstar_line, 10**(0.54+ 0.57*(mstar_line-11.)), 'r--',linewidth=2,alpha=0.8)
    plt.text(9.56, 10**(-0.34), 'z = 0.06', color='red', fontsize=35)    
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
sys.path.insert(0,'../py_tools')
from load_result import load_zs, load_mag, load_re, load_n, load_flux

ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281', 'CID255']
zs = np.asarray(load_zs(ID))
mags = np.array(load_mag(ID, folder = '../')[0])
Reffs = np.array(load_re(ID, folder = '../'))[:,0]
Reffs_e = np.array(load_re(ID, folder = '../'))[:,1]
indexs = np.array(load_n(ID, folder = '../'))[:,0]

dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #in pc
da=1/(1+zs)*c*vec_EE(zs)/h0   #in Mpc
ID_Reff_kpc = da * 10 **3 * (Reffs/3600./180.*np.pi)
ID_Reff_kpc_e = da * 10 **3 * ((Reffs_e)/3600./180.*np.pi)

from load_result import load_host_p

Mstar = load_host_p(ID, folder = '../')[1]
host_flux_WFC3 = np.array(load_flux(ID, folder = '../', flt = 'WFC3'))[:,0]
host_flux_ACS = []
for i in range(len(ID)):
    ifexit = glob.glob('../analysis_ACS/{0}'.format(ID[i]))
    if ifexit!= []:
        host_flux_ACS.append(load_flux([ID[i]], folder = '../', flt = 'ACS')[0][0])
    else:
        host_flux_ACS.append(-99)
host_flux_ACS = np.asarray(host_flux_ACS)

#cl=plt.colorbar()          #cl take the inforamtion from the lastest plt
#cl.set_label('filter flux ratio, WFC3 / ACS',rotation=270,size=30)
#cl.ax.get_yaxis().labelpad=35     #the distance of the colorbar titel from bar
#cl.ax.tick_params(labelsize=30)
if relation == 0:
    f1 ='../M_BH_relation/data/Bennert+2011_local.txt'
    b11_l = np.loadtxt(f1)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
    b11_local_Reff = b11_l[:,-1]
    b11_local_mstar = b11_l[:,4]
    plt.scatter(b11_local_mstar,b11_local_Reff,s=180, c ='black',
                marker="o",zorder=100, vmin=0.5, vmax=5, edgecolors='white', label='local AGN (VB2011)')     
    plt.scatter(Mstar[host_flux_ACS>0],ID_Reff_kpc[host_flux_ACS>0],s=200, linewidth='2', c =indexs[host_flux_ACS>0],
                marker="D",zorder=101, vmin=0.5, vmax=5, edgecolors='black', label='our AGN sample, 1.2<z<1.7',cmap=cmap_r)    
    for i in range(len(Mstar)):
        if Reffs[i]-0.1 < 0.009:
            plt.arrow(Mstar[i], ID_Reff_kpc[i], 0, -0.3, length_includes_head=True,
                  head_width=0.08, head_length=0.05, zorder=102, color='black', linewidth='1.2')
    
    log_Rerr = (np.log10(ID_Reff_kpc)-np.log10(ID_Reff_kpc-ID_Reff_kpc_e))
    low_err = ID_Reff_kpc - 10**(np.log10(ID_Reff_kpc)-log_Rerr)
    up_err = 10**(np.log10(ID_Reff_kpc)+log_Rerr) - ID_Reff_kpc
    plt.errorbar(Mstar[host_flux_ACS>0],ID_Reff_kpc[host_flux_ACS>0],
                 yerr=  [low_err[host_flux_ACS>0],
                         up_err[host_flux_ACS>0]],
#                 yerr= 10**(np.log10(ID_Reff_kpc)-np.log10(ID_Reff_kpc-ID_Reff_kpc_e)) [host_flux_ACS>0],
                 color='k',ecolor='k', fmt='.',markersize=1, zorder = 99)  
    
    plt.scatter(Mstar[host_flux_ACS<0],ID_Reff_kpc[host_flux_ACS<0],s=200, linewidth='2', c =indexs[host_flux_ACS<0],
                marker="D",zorder=101, vmin=0.5, vmax=5, edgecolors='black',cmap=cmap_r)
    plt.errorbar(Mstar[host_flux_ACS<0],ID_Reff_kpc[host_flux_ACS<0],
                 yerr=  [low_err[host_flux_ACS<0],
                         up_err[host_flux_ACS<0]],                 
#             yerr= (np.log10(ID_Reff_kpc)-np.log10(ID_Reff_kpc-ID_Reff_kpc_e))[host_flux_ACS<0],
             color='k',ecolor='k', fmt='.',markersize=1, zorder = 99)  
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
plt.xlabel("log(M$_*$/M$_{\odot})$",fontsize=35)
plt.tick_params(labelsize=25)
plt.legend(loc='upper right',fontsize=21,numpoints=1)
if relation ==0:
    plt.ylabel(r"R$_{\rm eff}$ (kpc)",fontsize=35)
#    plt.title(r"M$_*$ - R$_{eff}$ relation, sample redshift range {0}".format(z_range), fontsize = 25)
#    plt.title(r"M$_*$ - R$_{\rm eff}$ relation"+', sample redshift range {0}'.format(z_range), fontsize = 25)
    plt.ylim([0.3, 31.5])
#    labels = [item.get_text() for item in ax.get_yticklabels()]
#    labels[1] = 'Testing'
#    ax.set_xticklabels(labels)
    plt.yscale('log')
    ax.tick_params(axis='both', which='major', length=10)
    ax.tick_params(axis='y', which='minor', length=5)
    ax.tick_params(axis='y', which='both', width=1.5)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=20)
    cbar.ax.set_ylabel('Sersic index', rotation=270, fontsize = 25, labelpad=25)
    plt.savefig('Mstar-Reff.pdf')
#    plt.savefig('Mstar-Reff_z{0}-{1}.pdf'.format(z_range[0],z_range[1]))
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
    plt.savefig('Reff-z_{0}-{1}_color.pdf'.format(z_range[0],z_range[1]))

#cl.ax.tick_params(labelsize=15)   #the labe size
plt.show()

#%%KS test for their relations.
Reff_highz = np.log10(ID_Reff_kpc)

Mstar_blue, Reff_blue = results[:,3][z_cut* blue_galaxy],np.log10(Reff_kpc[z_cut * blue_galaxy])
Mstar_red, Reff_red = results[:,3][z_cut* red_galaxy],np.log10(Reff_kpc[z_cut * red_galaxy])
Reff_blue = Reff_blue[(Mstar_blue>9.5) * (Mstar_blue<12)]
Reff_red = Reff_red[(Mstar_red>9.5) * (Mstar_red<12)]

fig, ax = plt.subplots(figsize=(9,7))
high0, x0, _ = plt.hist(Reff_highz,normed=True, histtype=u'step',
         label=('high-z galaxy'), linewidth = 2, color='firebrick')
high1, x1, _ = plt.hist(Reff_blue,normed=True, histtype=u'step',
         label=('star-forming galaxy'), linewidth = 2, color='lightskyblue')
high2, x2, _ = plt.hist(Reff_red,normed=True, histtype=u'step',
         label=('quiescent galaxy'), linewidth = 2, color='darksalmon')
x0_m = np.median(Reff_highz)
high_m0 = high0[np.where(abs(x0_m-x0) == abs(x0_m-x0).min())[0][0]]
x1_m = np.median(Reff_blue)
high_m1 = high1[np.where(abs(x1_m-x1) == abs(x1_m-x1).min())[0][0]]
x2_m = np.median(Reff_red)
high_m2 = high2[np.where(abs(x2_m-x2) == abs(x2_m-x2).min())[0][0]-1]

plt.plot(np.linspace(0,high_m0)*0+np.median(x0_m) , np.linspace(0,high_m0), linewidth = 4,color='firebrick')
plt.plot(np.linspace(0,high_m1)*0+np.median(x1_m) , np.linspace(0,high_m1), linewidth = 4, color='lightskyblue')
plt.plot(np.linspace(0,high_m2)*0+np.median(x2_m) , np.linspace(0,high_m2), linewidth = 4, color='darksalmon')
plt.text(np.median(x0_m)-0.1, high_m0*1.05, '{0}'.format(round(np.median(x0_m),2)), color='firebrick',fontsize=25)
plt.text(np.median(x1_m)-0.2, high_m1*1.05, '{0}'.format(round(np.median(x1_m),2)), color='lightskyblue',fontsize=25)
plt.text(np.median(x2_m)-0.2, high_m2*1.05, '{0}'.format(round(np.median(x2_m),2)), color='darksalmon',fontsize=25)
fig.canvas.draw()
labels = [item.get_text().encode('ascii', 'replace').replace('?','-') for item in ax.get_xticklabels()]
print labels
for i in range(len(labels)-2):
    labels[i+1] = '10$^{'+ labels[i+1] + '}$'
ax.set_xticklabels(labels)

plt.ylim([0,3])
plt.xlabel("log(Reff) kpc",fontsize=27)
plt.ylabel("Density",fontsize=27)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.show()

from scipy import stats
print "p-value: high_z VS star-forming:", round(stats.ks_2samp(Reff_highz, Reff_blue).pvalue,3)
print "p-value: high_z VS quiesent:", round(stats.ks_2samp(Reff_highz, Reff_red).pvalue,3)
