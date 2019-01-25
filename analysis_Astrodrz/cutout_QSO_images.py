#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 09:11:54 2019

@author: Dartoon

Cutout the QSO from the drizzed HST image and share with Federica.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../py_tools')
import os
path = os.getcwd()
from cut_image import cut_image, cut_center_bright, save_loc_png, grab_pos
import copy
import glob
import re

IDs = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID1281', 'CID597'\
]

for i in range(len(IDs)):
    ID = IDs[i]
    filename= '{0}/analysis/{0}.reg'.format(ID)
    if_file = glob.glob(filename)
    if if_file != []:
        c_psf_list, QSO_loc = grab_pos(filename,reg_ty = 'astrodrz_06', QSO_reg_return=True)
        center_QSO = c_psf_list[QSO_loc]
    else:
        filename= '{0}/analysis/cut_PSFs.py'.format(ID)
        with open(filename) as f:
            content = f.readlines()
        lines = [x.strip() for x in content] 
        infos = [lines[i] for i in range(len(lines)) if 'center_QSO =' in lines[i]]
        if len(infos) > 1:
            print "Warning the center_QSO is not clear"
        info = infos[0].split( )
        info = [info[i] for i in range(len(info)) if 'np.array' in info[i]]
        center_QSO = np.asarray(re.findall(r'\d+', info[0]), dtype=int)
    fitsFile = pyfits.open('{0}/astrodrz/final_drz.fits'.format(ID))
    print ID
    img = fitsFile[1].data # check the back grounp
    
    QSO_sm, cut_center = cut_center_bright(image=img, center=center_QSO, radius=60, return_center=True, plot=False)
    QSO_outer = cut_image(image=img, center=cut_center, radius=200)
    
    file_header = copy.deepcopy(fitsFile[1].header)
    file_header['CRPIX1'] = file_header['CRPIX1']-center_QSO[0]+len(QSO_outer)/2
    file_header['CRPIX2'] = file_header['CRPIX2']-center_QSO[1]+len(QSO_outer)/2
    pyfits.PrimaryHDU(QSO_outer, header=file_header).writeto('QSOs_cutout/{0}_cutout.fits'.format(ID),overwrite=True)