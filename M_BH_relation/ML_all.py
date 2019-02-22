#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 15:50:27 2018

@author: Dartoon
"""
import numpy as np
np.set_printoptions(precision=4)
import matplotlib.pyplot as plt
import matplotlib as mat
import matplotlib.lines as mlines
from matplotlib import colors
mat.rcParams['font.family'] = 'STIXGeneral'
plt.figure(figsize=(14.5,12))

import matplotlib as mpl
mpl.rc('image', cmap='jet')

import sys
sys.path.insert(0,'../py_tools')
from dmag import pass_dmag
########## input local data ####
#==============================================================================
# The seleting for dm and host_total and dmag are in this local
#==============================================================================
from local import *
#==============================================================================
#input Park's data 
#==============================================================================
################ bulge or total relaiton? #################
#host= input('with bulge or total relation??? 0 = bulge;   1= total:')
####### input Park data ####
#######in AB system, V band#######
if host == 0:
   f0 ='data/Park_b'
if host == 1:
   f0 ='data/Park_t'
pk = np.loadtxt(f0)  #0 redshift; 1 BH mass; 2 Lv,obs,solar=0.4*(4.83-M_v)
pk[:,2]=4.83-pk[:,2]/0.4  # tansfer from L to Magnitude 
pk[:,2]=pk[:,2]-0.46  # transfer from V to R; 0.46 is the mean value of the data from Taka
pk[:,2]=pk[:,2]+dm*pass_dmag(pk[:,0])  #evolution of stellar population makes the mag fainter.
pk[:,2]=0.4*(4.61-pk[:,2])
plt.scatter(pk[:,2],pk[:,1],c=pk[:,0],marker="^",s=80,zorder=100,vmin=0.3, vmax=5, edgecolors='k')

#==============================================================================
# input Peng's data
#==============================================================================
inp_peng = 1

from cal_peng import *

Mg_cali=vec_Mgcal(L_Mg,FWHM_Mg)
Mg[:,3]=Mg[:,3]+0.21+dm*pass_dmag(Mg[:,1])
Hb_cali=vec_Hcal(L_hb,FWHM_hb)
H[:,3]=H[:,3]+0.21+dm*pass_dmag(H[:,1])   #change Mg from Vega to AB
C_cali=vec_Ccal(L_c,FWHM_c)
C[:,3]=C[:,3]+0.21+dm*pass_dmag(C[:,1])   #change Mg from Vega to AB
##################plot#############################
###############plot the data in Peng################
#For lensed samples
C[:,3],Mg[:,3],H[:,3]=0.4*(4.61-C[:,3]),0.4*(4.61-Mg[:,3]),0.4*(4.61-H[:,3])  #transfer from Mag to Luminosity

if inp_peng == 1:
    plt.scatter(C[:21,3],C_cali[:21],c=C[:21,1],s=80,marker="o",zorder=100, vmin=0.3, vmax=5, edgecolors='k')
    plt.scatter(Mg[:18,3],Mg_cali[:18],c=Mg[:18,1],s=80,marker="o",zorder=100, vmin=0.3, vmax=5, edgecolors='k')
    plt.scatter(H[:,3],Hb_cali,c=H[:,1],s=80,marker="o",zorder=100, vmin=0.3, vmax=5, edgecolors='k')  #zorder is the position at z foreground.
    
    #For non-lensed samples
    plt.scatter(C[21:,3],C_cali[21:],c=C[21:,1],s=80,marker="s",zorder=100, vmin=0.3, vmax=5, edgecolors='k')
    plt.scatter(Mg[18:,3],Mg_cali[18:],c=Mg[18:,1],s=80,marker="s",zorder=100, vmin=0.3, vmax=5, edgecolors='k')

#==============================================================================
# Two sample in H0licow VII: i.e. HE0435 and RXJ1131
#==============================================================================
###########HE0435 information using Mg#############
z_0435,L_0435,FWHM_0435=1.69,46-np.log10(5.15),4930
M_0435=Mgcal(L_0435,FWHM_0435)
mag_0435=21.75 -1.2514 + dm*pass_dmag(z_0435) #to vega for K correction, -1.2514 is the filter at F160w
Mag_0435=mag_0435-(20.09-(-23.40))  #23.40 and 20.09 are in Vege system, the K-correction value is the same as AB.
Mag_0435=Mag_0435+0.21     #change Mg from Vega to AB in R band
L_0435=0.4*(4.61-Mag_0435)
#print "L0435=",Mag_0435
plt.scatter(L_0435,M_0435,c=1.69,s=880,marker="*",zorder=100, vmin=0.3, vmax=5, edgecolors='k')
HE0435= mlines.Line2D([], [], color='m', ls='', marker='*', markersize=28, markeredgecolor='k')

###########RXJ1131 information using Mg#############
z_1131,L_1131,FWHM_1131=0.65,45-np.log10(5.15),5630
M_Mg_1131=Mgcal(L_1131,FWHM_1131)
###########RXJ1131 information using Hb#############
z_1131,L_1131,FWHM_1131=0.65,45-np.log10(3.81),4545.
M_Hb_1131=Hcal(L_1131,FWHM_1131)
from RXJ_Mr import Mr_bulge, Mr_disk ###1=bulge; 2=disk, K_correction inferred by Takahiro. For the details see RXJ_Mr.py
Mag_1131_1= Mr_bulge   
Mag_1131_2= Mr_disk
#print "L_bulge_1131=",0.4*(4.61-Mag_1131_1),"L_disk_1131=",0.4*(4.61-Mag_1131_2) 
Mag_1131_tot=-2.5*np.log10(10**((Mag_1131_1+dm*pass_dmag(0.65))/-2.5)+10**((Mag_1131_2+dm*pass_dmag(0.65))/-2.5))  # with brightness correction
Mag_1131_b=Mag_1131_1+dm*pass_dmag(0.65)  #-21.85 for the bulge, -23.183 for the disk
M_1131=np.log10((10**M_Hb_1131+10**M_Mg_1131)/2)
if host==0:
    Mag_1131=Mag_1131_b
if host==1:
    Mag_1131=Mag_1131_tot
L_1131 = 0.4*(4.61-Mag_1131)
plt.scatter(L_1131,M_1131,c=0.65,s=880,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
RXJ1131= mlines.Line2D([], [], color='lightpink', ls='', marker='*', markersize=28, markeredgecolor='k')

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
plt.scatter(lumi_s[MBs!=-99],MBs[MBs!=-99],c=zs[MBs!=-99],s=880,marker="x",zorder=100, vmin=0.3, vmax=5, edgecolors='k')

#==============================================================================
# The colorbar label setting up
#==============================================================================
cl=plt.colorbar()          #cl take the inforamtion from the lastest plt
cl.set_label('Source redshift',rotation=270,size=20)
cl.ax.get_yaxis().labelpad=35     #the distance of the colorbar titel from bar
cl.ax.tick_params(labelsize=30)   #the labe size

if host == 0:
  plt.xlabel("$log(L_{R,bulge}/L_{\odot})$",fontsize=35)
if host == 1:
  plt.xlabel("$log(L_{R,total}/L_{\odot})$",fontsize=35)
plt.ylabel("$log(M_{BH}/M_{\odot})$",fontsize=35)
plt.xticks(np.arange(5,16.5,0.5))
plt.yticks(np.arange(6,16.5,0.5))
plt.axis([8.5,12.5,6.25,11])
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)


lens = mlines.Line2D([], [], color='cyan', ls='', marker='o', markersize=9)
unlens = mlines.Line2D([], [], color='cyan', ls='', marker='s', markersize=9)
park = mlines.Line2D([], [], color='blue', ls='', marker='^', markersize=9)
new_sample = mlines.Line2D([], [], color='royalblue', ls='', marker='x', markersize=20)


if inp_peng == 1:
    plt.legend([HE0435,RXJ1131,park,lens,unlens,Pkc, new_sample],['HE0435',\
    'RXJ1131',\
    "AGNs from P15",
    #'M$_{BH}$ calibrated from Mg$_{II}$',\
    #"M$_{BH}$ calibrated from H${\\beta}$",'M$_{BH}$ calibrated from C$_{IV}$'
    "lensed AGNs in P06","non-lensed AGNs in P06"
    ,"local AGNs",
    "our new samples"\
    ],scatterpoints=1,numpoints=1,loc=2,prop={'size':24},ncol=2)
else:
    plt.legend([HE0435,RXJ1131,park,Pkc, new_sample],['HE0435',\
    'RXJ1131',\
    "AGNs from P15",
    #'M$_{BH}$ calibrated from Mg$_{II}$',\
    #"M$_{BH}$ calibrated from H${\\beta}$",'M$_{BH}$ calibrated from C$_{IV}$'
    "local AGNs",
    "our new samples"\
    ],scatterpoints=1,numpoints=1,loc=2,prop={'size':24},ncol=2)
#plt.savefig("MBH-L.pdf")
plt.show()
