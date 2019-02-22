#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 11:57:20 2019

@author: Dartoon

Compare with Aklant
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib as mat
import matplotlib.lines as mlines
from matplotlib import colors
mat.rcParams['font.family'] = 'STIXGeneral'
import matplotlib as mpl
mpl.rc('image', cmap='jet')

import sys
sys.path.insert(0,'../py_tools')
from dmag import pass_dmag

#==============================================================================
# My new inference
#==============================================================================
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

from load_result import load_zs, load_mag

ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
zs = np.asarray(load_zs(ID))
mags = np.array(load_mag(ID)[0])

from dmag import k_corr_R
import sys
sys.path.insert(0,'../py_tools')
from filter_info import filt_info
dm_k_R = []
for i in range(len(zs)):
    dm_k_R.append(k_corr_R(zs[i],filt_info[ID[i]], galaxy_age = '1Gyrs'))
dm_k_R = np.asarray(dm_k_R) # Get the k-correction for each target as an array
dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #in pc
host_mags=mags -5*(np.log10(dl)-1) + dm_k_R + dm*pass_dmag(zs)    # 0.7 is the k-correction value
#R=mags -5*(np.log10(dl)-1) +0.717      
print host_mags

lumi_s = 0.4*(4.61-host_mags)

f = open("fmos_MBH_table","r")
with f as g:
    lines = g.readlines()
porp_list = lines[0].replace('#','').split(' ')
samples = [lines[i].split(' ') for i in range(1,len(lines))]
ID_ser_dic =  {}

#XID2202 to LID1622
#XID2138 to LID1820
#XID2396 to LID1878
#CDFS321 to ECDFS321
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
for tar_in in range(len(ID)):       
    t_name = ID[tar_in]
    ser = ID_ser_dic[t_name]
#    print ser
    if ser!=-99 and float(samples[ser][10]) != 0:
        FWMH_a = float(samples[ser][8])
        logLHadr = float(samples[ser][6])
        cal_logMa = 6.71+0.48*(logLHadr-42)+2.12*np.log10(FWMH_a/1000)  # as used in Andreas
#        cal_logMa = 6.459+0.55*(logLHadr-42)+2.*np.log10(FWMH_a/1000)  # as used in H0liCOW 7 and McGill
        MBs.append(cal_logMa)
    elif ser!=-99 and float(samples[ser][21]) != 0:
        print "use Hb for", ID[tar_in]
        FWMH_b = float(samples[ser][19])
        logL5100dr = float(samples[ser][16])
        cal_logMb = 6.91+0.5*(logL5100dr-44)+2.*np.log10(FWMH_b/1000)  # as used in Andreas
#        cal_logMb = 6.882+0.518*(logL5100dr-44)+2.*np.log10(FWMH_b/1000)        # calibrated in H0liCOW 7
        MBs.append(cal_logMb)
    elif ser==-99:
        MBs.append(-99)
#        print float(cal_logMa) - float(samples[ser][10])
MBs = np.asarray(MBs)
#plt.scatter(lumi_s,MBs,c=zs,s=880,marker="p",zorder=100, vmin=0.3, vmax=5, edgecolors='k')

plt.subplots(figsize=(9.4,7.8))

m = np.linspace(7,9.2,30)
mr = np.linspace(-19.5,-22.6,30)
plt.plot(m,mr,'r')

plt.scatter(MBs[MBs!=-99],host_mags[MBs!=-99],c=zs[MBs!=-99],s=280,marker="o",zorder=100, vmin=0.3, vmax=5, edgecolors='k')
plt.gca().invert_yaxis()
plt.xlim(7,9.2)
plt.ylim(-17.5, -24.5)
plt.tick_params(labelsize=25)
cl=plt.colorbar()          #cl take the inforamtion from the lastest plt
plt.clim(1.2,1.8)
cl.set_label('Source redshift',rotation=270,size=20)
cl.ax.get_yaxis().labelpad=35     #the distance of the colorbar titel from bar
cl.ax.tick_params(labelsize=30)   #the
plt.xlabel("$log(M_{BH}/M_{\odot})$",fontsize=35)
plt.ylabel("$M_R$",fontsize=35)
plt.grid()
plt.savefig('/Users/Dartoon/Downloads/HST_sample_Mr-Mbh.pdf')
plt.show()