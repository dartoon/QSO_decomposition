#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 20:38:51 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from subprocess import call
#print call('ls', shell=True)

filt_info = {'CDFS-1': 'F140w', 'CDFS-229': 'F125w', 'CDFS-724': 'F125w',\
'CID1174': 'F140w', 'CID206': 'F140w', 'CID216': 'F140w', 'CID3242': 'F140w',\
'CID3570': 'F125w', 'CID452': 'F125w', 'CID454': 'F140w', 'CID50': 'F125w',\
'CID607': 'F125w', 'CID70': 'F140w', 'LID1273': 'F140w', 'LID360': 'F140w',\
'XID2138': 'F140w', 'XID2202': 'F140w', 'XID2396': 'F140w', 'ECDFS-358': 'F140w',\
'SXDS-X1136': 'F125w', 'SXDS-X50': 'F125w', 'SXDS-X735': 'F140w',\
'CID543': 'F125w', 'LID1538': 'F140w', 'CID237': 'F140w', 'SXDS-X717': 'F125w',\
'SXDS-X763': 'F125w', 'SXDS-X969': 'F140w', 'CDFS-321': 'F140w'}

#for key in filt_info.keys():
#    ID = key
#    print key
#    print call("mkdir {0}".format(ID), shell=True)
#    print call("cp -r ../template/temp_analysis/Drz06_temp/* {0}".format(ID), shell=True)
#    print call("cp ../Cycle25data/{0}/*_flt.fits {0}/astrodrz".format(ID), shell=True)

#for key in filt_info.keys():
#    ID = key
#    print ID
#    try:
#        runfile('/Users/Dartoon/Astro/QSO_decomposition/analysis_Astrodrz/{0}/astrodrz/img_bkg_sub.py'.format(ID), wdir='/Users/Dartoon/Astro/QSO_decomposition/analysis_Astrodrz/{0}/astrodrz'.format(ID))
#    except:
#        print(" ")

#import glob
#for key in filt_info.keys():
#    ID = key
##    in_fits = glob.glob("{0}/astrodrz/i*flt.fits".format(ID))
#    in_fits = glob.glob("{0}/astrodrz/final_drz.fits".format(ID))
#    hdul = pyfits.open(in_fits[0])
##    print hdul.info()
#    print ID,":", hdul[1].header['ORIENTAT']

#for key in filt_info.keys():
#    ID = key
#    if ID not in ['XID2202', 'CID50', 'CID1174', 'CID206']:
#        print key
#        print call("cp ../template/temp_analysis/Drz06_temp/analysis/cut_PSFs.py {0}/analysis/".format(ID), shell=True)

#for key in filt_info.keys():
#    ID = key
#    if ID not in ['CDFS-1']:
#        print ID
#        print call("cp {0}/analysis/4_analysis.py {0}/analysis/best_fit_plot_for_paper.py".format(ID), shell=True)
    
