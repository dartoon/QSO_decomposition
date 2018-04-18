#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 22:52:55 2018

@author: dartoon

PSF compare
"""

import astropy.io.fits as pyfits
import sys
sys.path.insert(0,'../../py_tools')
from flux_profile import profiles_compare
import matplotlib.pyplot as plt

PSF1 = pyfits.getdata('PSF1.fits')
PSF1_swarp = pyfits.getdata('../../analysis_SWarp/CID70/analysis/PSF1.fits')

prf_list = [PSF1,PSF1_swarp]
scal_list = [1,1]
prf_name_list = ['PSF1_by_HST_drz', 'PSF1_by_SWarp']
fig_pro_compare = profiles_compare(prf_list, scal_list, prf_name_list=prf_name_list, gridspace = 'log',if_annuli=True)
fig_pro_compare.savefig('PSFavd_vs_QSO_xlin_annu1.pdf')
plt.show()
