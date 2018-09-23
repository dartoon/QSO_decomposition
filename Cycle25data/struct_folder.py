#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 23 10:04:02 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

id_list= ['idnl28', 'idnl29', 'idnl32', 'idnl09', 'idnl02', 'idnl03', 'idnl11', 'idnl12', 'idnl06', 'idnl17', 'idnl01', 'idnl08', 'idnl14', 'idnl15', 'idnl30', 'idnl20', 'idnl16', 'idnl13', 'idnl31', 'idnl26', 'idnl18', 'idnl23', 'idnl19', 'idnl21', 'idnl04', 'idnl22', 'idnl24', 'idnl25', 'idnl27']
name_list = ['CDFS-1', 'CDFS-229', 'CDFS-724', 'CID1174', 'CID206', 'CID216', 'CID3242', 'CID3570', 'CID452', 'CID454', 'CID50', 'CID607', 'CID70', 'LID1273', 'LID360', 'XID2138', 'XID2202', 'XID2396', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X735', 'CID543', 'LID1538', 'CID237', 'SXDS-X717', 'SXDS-X763', 'SXDS-X969', 'CDFS-321']


import os
path = os.path.dirname(os.path.realpath(__file__))

import subprocess 

for i in range(len(id_list)):
    name = path+'/'+name_list[i]
    os.mkdir(name)
    files_name = id_list[i]+'*'
    subprocess.call("mv {0} {1}".format(files_name,name), shell=True)
#    shutil.move(id_list[i]+'*', name)