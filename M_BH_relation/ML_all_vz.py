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
host=plt.figure(figsize=(14.5,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label

import matplotlib as mpl
mpl.rc('image', cmap='jet')
import sys
sys.path.insert(0,'../py_tools')
from dmag import pass_dmag

inp_peng = 1

########## input local data ####
#==============================================================================
# The seleting for dm and host_total and dmag are in this local
#==============================================================================
from local_vz import *
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
park=plt.errorbar(np.log10(1+pk[:,0]),pk[:,1]-(m_mid*pk[:,2]+b_mid),yerr=(0.4**2+0.2**2)**0.5,fmt='.',color='darkseagreen')#,markersize=1)
SS13=plt.errorbar(np.log10(1+pk[63:,0]),pk[63:,1]-(m_mid*pk[63:,2]+b_mid),yerr=(0.4**2+0.2**2)**0.5,fmt='.',color='tomato')#,markersize=1)

#==============================================================================
# input Peng's data
#==============================================================================
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
    plt.errorbar(np.log10(1+Mg[:,1])[:14],(Mg_cali-(m_mid*Mg[:,3]+b_mid))[:14],yerr=(0.4**2+0.2**2)**0.5,fmt='o',color='goldenrod',markersize=9,zorder=250, mec='k')
    plt.errorbar(np.log10(1+H[:,1]),Hb_cali-(m_mid*H[:,3]+b_mid),yerr=(0.4**2+0.2**2)**0.5,fmt='>',color='goldenrod',markersize=9,zorder=200, mec='k')
    plt.errorbar(np.log10(1+C[:,1])[:18],(C_cali-(m_mid*C[:,3]+b_mid))[:18],yerr=(0.4**2+0.2**2)**0.5,fmt='s',color='goldenrod',markersize=9,zorder=250, mec='k')
    
    plt.errorbar(np.log10(1+Mg[:,1])[14:],(Mg_cali-(m_mid*Mg[:,3]+b_mid))[14:],yerr=(0.4**2+0.2**2)**0.5,fmt='o',color='dodgerblue',markersize=9,zorder=250, mec='k')
    plt.errorbar(np.log10(1+C[:,1])[18:],(C_cali-(m_mid*C[:,3]+b_mid))[18:],yerr=(0.4**2+0.2**2)**0.5,fmt='s',color='dodgerblue',markersize=9,zorder=250, mec='k')

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
plt.errorbar(np.log10(1+z_0435),M_0435-(m_mid*L_0435+b_mid),yerr=(0.4**2+0.2**2)**0.5,fmt='m*',markersize=28,zorder=250, mec='k')
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
plt.errorbar(np.log10(1+z_1131),M_1131-(m_mid*L_1131+b_mid),yerr=(0.4**2+0.2**2)**0.5,fmt='*',color='lightpink',markersize=28,zorder=250, mec='k')
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
from filter_info import filt_info
dm_k_R = []
for i in range(len(zs)):
    dm_k_R.append(k_corr_R(zs[i],filt_info[ID[i]], galaxy_age = '1Gyrs'))
dm_k_R = np.asarray(dm_k_R)   # Get the k-correction for each target as an array
dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #in pc
host_mags=mags-5*(np.log10(dl)-1) + dm_k_R + dm*pass_dmag(zs)  # 0.7 is the k-correction value
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
plt.scatter(np.log10(1+zs[MBs!=-99]),MBs[MBs!=-99]-(m_mid*lumi_s[MBs!=-99]+b_mid),c='royalblue',s=880,marker="x",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
#plt.errorbar(np.log10(1+zs[MBs!=-99]),MBs[MBs!=-99]-(m_mid*lumi_s[MBs!=-99]+b_mid),yerr=(0.4**2+0.2**2)**0.5,fmt='x',color='royalblue',markersize=28,zorder=250)#, mec='k')

#####fit the evolution##########
################################
#z_loc,y_loc,err_loc=0,0,0
z_pk,y_pk=pk[:,0],pk[:,1]-(m_mid*pk[:,2]+b_mid)
z_Mg,y_Mg=Mg[:,1],Mg_cali-(m_mid*Mg[:,3]+b_mid)
z_H,y_H=H[:,1],Hb_cali-(m_mid*H[:,3]+b_mid)
z_C,y_C=C[:,1],C_cali-(m_mid*C[:,3]+b_mid)
z_0435,y_0435=1.69,M_0435-(m_mid*L_0435+b_mid)
z_1131,y_1131=0.65,M_1131-(m_mid*L_1131+b_mid)
z_cosmos, y_cosmos = zs[MBs!=-99], MBs[MBs!=-99]-(m_mid*lumi_s[MBs!=-99]+b_mid)
                              
z_1131=0.65

if inp_peng == 1:
    z=np.concatenate((z_pk,z_Mg,z_H,z_C,np.array([z_0435,z_1131]), z_cosmos),axis=0)
    y=np.concatenate((y_pk,y_Mg,y_H,y_C,np.array([y_0435,y_1131]), y_cosmos),axis=0)
else:
    z=np.concatenate((z_pk,np.array([z_0435,z_1131]), z_cosmos),axis=0)
    y=np.concatenate((y_pk,np.array([y_0435,y_1131]), y_cosmos),axis=0)


#z=z_cosmos #np.concatenate((z_pk,z_Mg,z_H,z_C,np.array([z_0435,z_1131]), z_cosmos),axis=0)
#y=y_cosmos #np.concatenate((y_pk,y_Mg,y_H,y_C,np.array([y_0435,y_1131]), y_cosmos),axis=0)


#### fit with emcee ###############
x=np.log10(1+z)
y=y
yerr=y*0+(0.4**2+0.2**2)**0.5   # the error for the fitting
def lnlike(theta, x, y, yerr):
    b, sint= theta
    model = b*x
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [1.8, 0.3], args=(x, y, yerr))
b_ml,sint_ml= result["x"]
#print b_ml, sint_ml, "ka=",lnlike(theta=[b_ml, sint_ml],x=loc[:,0], y=loc[:,1], yerr=loc[:,2])

xp = np.array([5, 13])
#plt.plot(xp, m_ml*xp+b_ml, 'r-')
def lnprior(theta):
    b, sint	 = theta
    if -10 < b < 10.0 and 0 < sint < 10:
        return 0.0
    return -np.inf
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)
ndim, nwalkers = 2, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

b_mid, sint_mid =np.percentile(samples, 50,axis=0)
#print "lnlike=",lnlike(theta=[b_mid, sint_mid],x=x, y=y, yerr=yerr)
xl = np.linspace(0, 5, 100)
plt.plot(xl, xl*0+xl*b_mid, color="red", linewidth=4.0,zorder=0)

def find_n(array,value):           #get the corresponding b for a given m 
    idx= (np.abs(array-value)).argmin()
    return array[idx]
b=np.percentile(samples,50,axis=0)[0]
#print samples[:,1][samples[:,0]==find_n(samples[:,0],m)]
for i in range(100):
    posi=np.random.uniform(16,84)
    b=np.percentile(samples,posi,axis=0)[0]    
    #print b
    plt.plot(xl, xl*0+xl*b, color="lightgray", alpha=0.2,linewidth=7.0,zorder=-1)
value=round(b_mid,2)
#####################
value,sig=round(b_mid,2),round((np.percentile(samples,84,axis=0)[0]-np.percentile(samples,16,axis=0)[0])/2,2)
#print value,sig
#plt.text(0.15, -1.75, "$M_{BH} $VS$ L_{host}\propto (1+z)^{%s\pm%s}$"%(value,sig),color='blue',fontsize=25)
#"$M_{BH}-L_{host}\propto (1+z)

Mg = mlines.Line2D([], [], color='goldenrod', ls='', marker='o', markersize=9, markeredgecolor='k')
H = mlines.Line2D([], [], color='goldenrod', ls='', marker='>', markersize=9, markeredgecolor='k')
C = mlines.Line2D([], [], color='goldenrod', ls='', marker='s', markersize=9, markeredgecolor='k')
unlens = mlines.Line2D([], [], color='dodgerblue', ls='', marker='s', markersize=9, markeredgecolor='k')
new_sample = mlines.Line2D([], [], color='royalblue', ls='', marker='x', markersize=20)

#
plt.xlabel("log$(1+z)$",fontsize=35)
plt.ylabel("$\Delta$log$M_{BH}$",fontsize=35)
plt.text(0.15, -1.75, "$\Delta$log$M_{BH}$=$(%s\pm%s)$log$(1+z)$"%(value,sig),color='blue',fontsize=25)

plt.xticks(np.arange(-0.1,1,0.1))
plt.yticks(np.arange(-5.5,6,0.5))
xl=-0.01
xh=0.8
plt.axis([xl,xh,-2.0,3.5])
plt.grid()
plt.tick_params(labelsize=25)

ax2=ax.twiny()
tticks=np.array([10**xl-1,0.5,1,2,3,4,5,10**xh-1])
ax2.set_xticks([np.log(t+1) for t in tticks])  # for the entire scale
ax2.set_xticklabels([0,0.5,1,2,3,4,5])  # 0 actuall is corresponds to 10**-0.01-1
ax2.set_xlabel('$z$',fontsize=35)
plt.tick_params(labelsize=25)

if inp_peng == 1:
    plt.legend([HE0435,RXJ1131,park,SS13,Mg,H,C,unlens,Pkc, new_sample],['HE0435',\
    'RXJ1131',\
    "AGNs in Park 15",\
    "AGNs in Schramm et al. 2013",\
    'lensed AGNs in P06 using Mg$_{II}$',\
    "lensed AGNs in P06 using H${\\beta}$",\
    'lensed AGNs in P06 using C$_{IV}$',\
    "non-lensed AGNs in P06",\
    "local AGNs",\
    "our new samples"
    ],scatterpoints=1,numpoints=1,loc=2,prop={'size':23},ncol=2,handletextpad=0)
else:
    plt.legend([HE0435,RXJ1131,park,SS13,Pkc, new_sample],['HE0435',\
    'RXJ1131',\
    "AGNs from Park 15",\
    "AGNs from Schramm et al. 2013",\
    "local AGNs",\
    "our new samples"
    ],scatterpoints=1,numpoints=1,loc=2,prop={'size':23},ncol=2,handletextpad=0)
#plt.savefig("MBHL-z.pdf")
plt.show()
