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

import sys
sys.path.insert(0,'../py_tools')
from filter_info import filt_info, redshift_info

IDs = filt_info.keys()
IDs.sort()

# =============================================================================
# For re-fit the sample
# =============================================================================
for key in IDs:
    ID = key
    print key
    print call("mkdir {0}/deep_analysis".format(ID), shell=True)
    print call("mkdir {0}/deep_analysis/fit_result_each/".format(ID), shell=True)
    print call("cp {0}/analysis/{0}*.fits {0}/analysis/wht_err.fits {0}/analysis/*.pdf {0}/deep_analysis/".format(ID), shell=True)
    print call("cp {0}/analysis/4_analysis.py {0}/deep_analysis/".format(ID), shell=True)
    with open("{0}/deep_analysis/4_analysis.py".format(ID)) as f:
            contents = f.readlines()
    for i in range(len(contents)):
        if 'deep_seed' in contents[i] and 'False' in contents[i]:
            print ID, "is deep_seed = False", "line:", i

#    print call("cp ../Cycle25data/{0}/*_flt.fits {0}/astrodrz".format(ID), shell=True)
# =============================================================================
# Before
# =============================================================================
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
    
