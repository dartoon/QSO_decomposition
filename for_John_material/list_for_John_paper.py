#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 13:45:03 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../py_tools')

from load_result import load_host_p, load_MBH, load_err, load_re, load_zs

ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']
MB_ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'ECDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','LID1820','LID1622',\
'LID1878', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']

MBs = load_MBH(ID,MB_ID,if_reportHb=0)
Mstar = load_host_p(ID)[1]
Mstar_err = load_err(prop = 'Mstar', ID=ID)
zs = np.asarray(load_zs(ID))

Reffs = np.array(load_re(ID))

tab_list = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 'CID3570',
            'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 'CID70', 'LID1273',
            'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396', 'CDFS-1', 'CDFS-229',
            'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']

from scipy.integrate import quad
h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)
dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #in pc
da=1/(1+zs)*c*vec_EE(zs)/h0   #in Mpc
Reff_kpc = da * 10 **3 * (Reffs[:,0]/3600./180.*np.pi)

for target in tab_list:
    i = [i for i in range(len(ID)) if target == ID[i]][0]
#    print M_BH, M_BH_error(dex), M_stellar, M_stellar_error, R_e(arcsec), R_e_err, R_e(kpc), R_e_err, ID
    print round(MBs[i],2), 0.4, round(Mstar[i],2), round(Mstar_err[i][0],2), '+'+repr(round(Mstar_err[i][1],2)),\
    round(Reffs[i][0],2), round(Reffs[i][1],2), round(Reff_kpc[i],2), round(Reff_kpc[i]*(Reffs[i][1]/Reffs[i][0]),2), "#"+ID[i]
    
#%%
import sys
sys.path.insert(0,'../py_tools')
#==============================================================================
# My new inference
#==============================================================================
from load_result import load_host_p, load_MBH, load_err
from load_result import load_zs, load_n
ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']
MB_ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'ECDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','LID1820','LID1622',\
'LID1878', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']
zs = np.asarray(load_zs(ID))
host_n = np.array(load_n(ID, folder = '../'))[:,0]
Mstar = load_host_p(ID)[1]
import pickle
n_BT_relation = pickle.load(open('../Comparsion/CANDELS_catalog/bulge_disk_fit/n_BT_relation.pkl','rb'))
n_list,BTR_smooth_mean,BTR_smooth_median = n_BT_relation
idx = [np.where(abs(host_n[i] - n_list) == abs(host_n[i] - n_list).min())[0][0] for i in range(len(host_n))]
host_BT = np.array([BTR_smooth_mean[idx[i]] for i in range(len(host_n))])
for target in tab_list:
    i = [i for i in range(len(ID)) if target == ID[i]][0]
#    print M_BH, M_BH_error(dex), M_stellar, M_stellar_error, R_e(arcsec), R_e_err, R_e(kpc), R_e_err, ID
    print round(host_BT[i],3), "#"+ID[i]
   

