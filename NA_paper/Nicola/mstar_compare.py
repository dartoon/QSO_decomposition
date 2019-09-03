#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:58:05 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy, matplotlib
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
cmap = matplotlib.cm.get_cmap('viridis')
import matplotlib as mpl
mpl.rc('image', cmap='jet')

#%% Load Nicola's data
#filename = 'contour_all.dat'
filename = 'contour_eddvar.dat'
data = np.loadtxt(filename)
prop_host = data[:,0]
prop_mbh = data[:,1]
number_dens = data[:,2]
norm_number_dens = data[:,3]

'''
#%% Load the BH mass and the stellar mass
import sys
sys.path.insert(0,'../../py_tools')
from load_result import load_zs, load_mag, load_n
import re
ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
zs = np.asarray(load_zs(ID))
#lumi_s = 0.4*(4.61-host_mags)
f = open("../../M_BH_relation/fmos_MBH_table","r")
with f as g:
    lines = g.readlines()
porp_list = lines[0].replace('#','').split(' ')
samples = [lines[i].split(' ') for i in range(1,len(lines))]
ID_ser_dic =  {}

MB_ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'ECDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','LID1820','LID1622',\
'LID1878', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
for j in range(len(ID)):
    count = 0
    for i in range(len(samples)):
        if samples[i][1] == MB_ID[j]:
            ID_ser_dic.update({ID[j]:i})
            count += 1
    if count == 0:
        ID_ser_dic.update({ID[j]: -99})
MBs = []
CDFS_FWHMa = {'CDFS-1': 5449.4022,'CDFS-229': 2254.0105, 'CDFS-724': 3351.852239}
CDFS_logLHadr = {'CDFS-1': 43.08,'CDFS-229': 43.30, 'CDFS-724': 42.561413}
for tar_in in range(len(ID)):       
    t_name = ID[tar_in]
    ser = ID_ser_dic[t_name]
#    print ser
    if ser!=-99 and float(samples[ser][10]) != 0:
        FWMH_a = float(samples[ser][8])
        logLHadr = float(samples[ser][6])
        cal_logMa = 6.71+0.48*(logLHadr-42)+2.12*np.log10(FWMH_a/1000)  # as used in Andreas
#        cal_logMa = 6.459+0.55*(logLHadr-42)+2.*np.log10(FWMH_a/1000)  # as used in H0liCOW 7 and McGill
        mbh = cal_logMa
        MBs.append(mbh)
    if ser!=-99 and float(samples[ser][21]) != 0:
        print "use Hb for", ID[tar_in]
        FWMH_b = float(samples[ser][19])
        logL5100dr = float(samples[ser][16])
        cal_logMb = 6.91+0.5*(logL5100dr-44)+2.*np.log10(FWMH_b/1000)  # as used in Andreas
#        cal_logMb = 6.882+0.518*(logL5100dr-44)+2.*np.log10(FWMH_b/1000)        # calibrated in H0liCOW 7
        mbh = (cal_logMa + cal_logMb)/2
#        print mbh
        MBs[-1] = mbh
    if ser==-99 and t_name in ['CDFS-1','CDFS-229', 'CDFS-724']:
        FWMH_a = float(CDFS_FWHMa[t_name])
        logLHadr = float(CDFS_logLHadr[t_name])
        cal_logMa = 6.71+0.48*(logLHadr-42)+2.12*np.log10(FWMH_a/1000)  # as used in Andreas
        MBs.append(cal_logMa)
    elif ser==-99:
        MBs.append(-99)
#        print float(cal_logMa) - float(samples[ser][10])
MBs = np.asarray(MBs)
host_n = np.array(load_n(ID, folder = '../../'))[:,0]

# Load my sample and transfer to Stellar mass
from dmag import pass_dmag
from scipy.integrate import quad
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)
h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
dm = 0 #Not considering the passive evolution
mags = np.array(load_mag(ID, folder = '../../')[0])
from dmag import k_corr_R
from filter_info import filt_info
dm_k_R = []
for i in range(len(zs)):
    dm_k_R.append(k_corr_R(zs[i],filt_info[ID[i]], galaxy_age = '1Gyrs'))
dm_k_R = np.asarray(dm_k_R) # Get the k-correction for each target as an array
dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #in pc
host_Mags = mags -5*(np.log10(dl)-1) + dm_k_R + dm*pass_dmag(zs) # This is in AB system
host_Mags = host_Mags - 0.21  # Transfer to Vega system
host_LR = 10 ** (0.4*(4.61-host_Mags))
Mstar = np.log10(host_LR * 0.54 * 0.684 * 1.4191)  
# 0.54 is the ratio M*_tot/Lv, 
# 0.684 is the color correction Lv/Lr, 10**(-0.4*(V-R)) = 10**(-0.4*0.4119)
# 1.4191 is the Lr_sun/Lv_sun = 10**(-0.4*(Rsun-Vsun)) = 10**(-0.4*(4.43-4.81)

#SED_stellar = np.loadtxt('../SED_data/firstmasses.dat')
#SED_ID = np.array(SED_stellar[:,1], dtype=int)
#indexs = []
#goods, Mstar = [] , []  
#for tar_in in range(len(ID)):  
#    t_name = ID[tar_in]
#    numb = int(re.findall("\d+", t_name)[0])
#    if "SXDS" in t_name:
#        index = np.where(numb == SED_ID)[0]
#        index = index[index>24]
#    elif "CDFS" in t_name:
#        index = np.where(numb == SED_ID)[0]
#        index = index[index>19]
#    else:
#        index = np.where(numb == SED_ID)[0]
#        index = index[index<=19]
#    indexs.append(index)
#    goods.append(SED_stellar[index][0][0])
#    Mstar.append(SED_stellar[index][0][2])
#Mstar = np.asarray(Mstar)
#Mstar = np.log10(Mstar) + 10
#goods = np.asarray(goods)

#bools = ([MBs!=-99])
bools_b = ([MBs!=-99] and [host_n>3])[0]
bools_d = ([MBs!=-99] and [host_n<3])[0]
'''
#%% Plot the data
plt.figure(figsize=(13, 11))
vmax = 1.1 #1.1 #10**norm_number_dens.max()
vmin = -0.1 #10**norm_number_dens[norm_number_dens!=norm_number_dens.min()].min()
#my_cmap = copy.copy(matplotlib.cm.get_cmap('YlOrBr')) # copy the default cmap
my_cmap = copy.copy(matplotlib.cm.get_cmap('Oranges',6)) # copy the default cmap
#my_cmap.set_under('white')

for i in range(len(norm_number_dens)):
    if norm_number_dens[i]!=norm_number_dens.min() : 
        plt.scatter(prop_host[i], prop_mbh[i], c=10**norm_number_dens[i], s = 140, vmin = vmin, vmax = vmax,
                   marker='s', alpha=0.9, edgecolors='none', cmap = my_cmap)  # plot the Nicola
cl=plt.colorbar()
cl.set_label('Value in Col. 4', size=20)
#plt.scatter(Mstar[bools_d], MBs[bools_d],c=goods[bools_d],s=280,marker="s",zorder=100, vmin=0.8, vmax=3, edgecolors='k')
#plt.scatter(Mstar[bools_b], MBs[bools_b],c=goods[bools_b],s=280,marker="o",zorder=100, vmin=0.8, vmax=3, edgecolors='k')
#cl=plt.colorbar()
#cl.set_label('1 means good, 2 means bad', size=20)
#plt.clim(0.8,3)

#xi = prop_host[:50]
#yi = prop_mbh[::50]
#prop_host_cl = prop_host[norm_number_dens!=norm_number_dens.min()]
#prop_mbh_cl = prop_mbh[norm_number_dens!=norm_number_dens.min()]
#norm_number_dens_cl = norm_number_dens[norm_number_dens!=norm_number_dens.min()]
#plt.scatter(prop_host_cl, prop_mbh_cl, c=10**norm_number_dens_cl, s = 140, vmin = vmin, vmax = vmax,
#           marker='s', alpha=0.9, edgecolors='none', cmap = my_cmap)  # plot the Nicola
#import matplotlib.tri as tri
#triang = tri.Triangulation(prop_host_cl, prop_mbh_cl)
#interpolator = tri.LinearTriInterpolator(triang,10**norm_number_dens_cl)
#Xi, Yi = np.meshgrid(xi, yi)
#zi = interpolator(Xi, Yi)
#plt.contour(xi, yi, zi, 5,linewidths=0.5, colors='k')
#plt.contourf(xi, yi, zi, 5, cmap=my_cmap)

'''
plt.scatter(Mstar[bools_d], MBs[bools_d],c=zs[bools_d],s=280,marker="s",zorder=100, vmin=0.3, vmax=2, edgecolors='k')
plt.scatter(Mstar[bools_b], MBs[bools_b],c=zs[bools_b],s=280,marker="o",zorder=100, vmin=0.3, vmax=2, edgecolors='k')
'''
#cl=plt.colorbar()          #cl take the inforamtion from the lastest plt
#plt.clim(vmin=0.3, vmax=2)
#cl.set_label('Source redshift',rotation=270,size=20)

#x_aklant = np.linspace(10.5,11.8)
#y_aklant = np.linspace(7,9.2)
#plt.plot(x_aklant-0.25-0.15,y_aklant-0.15, x_aklant-0.15, y_aklant-0.15, x_aklant+0.25-0.15,y_aklant-0.15,color = 'gray')


#x_local = np.linspace(9,12.5)
#y_local = np.linspace(6, 9.9)
#plt.plot(x_local,y_local-0.4, x_local, y_local, x_local,y_local+0.4,color = 'black')


plt.xlabel("$log(M_*/M_{\odot})$",fontsize=35)
plt.ylabel("$log(M_{BH}/M_{\odot})$",fontsize=35)
plt.minorticks_on()
#plt.xlim([prop_host.min(), prop_host.max()])
#plt.ylim([prop_mbh.min(), prop_mbh.max()])
#plt.xlim([9, 12])
#plt.ylim([5, 10])
#from local_MMstar import *

if 'rmag' in filename:
    plt.gca().invert_xaxis()
plt.tick_params(labelsize=25)
plt.savefig('SAM_MMstar.pdf')
plt.show()

#for i in range(len(ID)):
#    print ID[i], round(MBs[i],3), round(Mstar[i],3)