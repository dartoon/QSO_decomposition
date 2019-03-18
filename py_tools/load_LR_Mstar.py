#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 19:09:55 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from load_result import load_zs, load_mag
from dmag import pass_dmag
import copy

h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]

from scipy.integrate import quad
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)

def load_host_p(ID, temp='1Gyr', dm = 0):
    zs = np.asarray(load_zs(ID))
    mags = np.array(load_mag(ID, folder = '../')[0])
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
    return np.log10(host_LR), Mstar

def load_MBH(ID, MB_ID):
    #ID change:
    #XID2202 to LID1622
    #XID2138 to LID1820
    #XID2396 to LID1878
    #CDFS321 to ECDFS321
    f = open("../M_BH_relation/fmos_MBH_table","r")
    with f as g:
        lines = g.readlines()
#    porp_list = lines[0].replace('#','').split(' ')
    samples = [lines[i].split(' ') for i in range(1,len(lines))]
    ID_ser_dic =  {}
    
    for j in range(len(MB_ID)):
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
    return MBs

#ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
#'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
#MB_ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'ECDFS-321', 'CID1174',\
#'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
#'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
#'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','LID1820','LID1622',\
#'LID1878', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281']
#MBs = load_MBH(ID, MB_ID)