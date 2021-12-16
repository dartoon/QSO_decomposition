#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 17:53:10 2019

@author: Dartoon

compare my inference to Sun's SED fitting.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'../../py_tools')
from load_result import load_host_p, load_MBH

sun_data = np.loadtxt('Sun_data.txt')
sun_pos = sun_data[:,1:3]
sun_z = sun_data[:,3]
sun_MBH = sun_data[:,-2]
sun_Mstar = sun_data[:,-1]

my_data = np.loadtxt('my_table1.txt')
my_z = my_data[:,0]
my_pos = my_data[:,1:3]

IDs = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 
       'CID3570', 'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 
       'CID70', 'LID1273', 'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396', 
       'CDFS-1', 'CDFS-229', 'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 
       'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']
MB_IDs = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 
       'CID3570', 'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 
       'CID70', 'LID1273', 'LID1538', 'LID360', 'LID1820', 'LID1622', 'LID1878', 
       'CDFS-1', 'CDFS-229', 'ECDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 
       'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']

Mstar = load_host_p(IDs, folder='../../')[1]
MBs = load_MBH(IDs,MB_IDs,if_reportHb=0,folder='../../')
plt.figure(figsize=(8,8))
for i in range(len(my_pos)):
    dis = my_pos[i] - sun_pos
    diff = np.sqrt(dis[:,0]**2 + dis[:,1]**2)
    if diff.min()*3600/0.0642 < 30: #diff.min()*3600/0.0642 is in the pixel scale in 0.0642
        print(my_z[i], sun_z[diff==diff.min()])
        print(MBs[i], sun_MBH[diff==diff.min()])
        print(Mstar[i], sun_Mstar[diff==diff.min()])
        print(IDs[i])
        plt.scatter(Mstar[i], sun_Mstar[diff==diff.min()])
#        plt.scatter(MBs[i], Mstar[i], marker="*", c='k')
#        plt.scatter(sun_MBH[diff==diff.min()], sun_Mstar[diff==diff.min()],marker="o", c='r')
x_line = np.linspace(10,12,10)
plt.plot(x_line, x_line)
plt.show()
        
        

        
    
    

