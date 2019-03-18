#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 13:48:52 2018

@author: Dartoon

Test if the Lenstronmy can derive the right answer of the simulated QSO.
"""

# import of standard python libraries
import numpy as np
import os
import time
import copy
import corner
import astropy.io.fits as pyfits

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#==============================================================================
# To do the fitting:
#==============================================================================
image_noisy = pyfits.getdata('sim_image.fits')
QSO_std = pyfits.getdata('error_map.fits')
#fix_n = 4.
fixed_source, kwargs_source_init,kwargs_source_sigma,kwargs_lower_source, kwargs_upper_source= [], [], [], [], []
fixed_source.append({'n_sersic': fix_n,'R_sersic': .3})
#fixed_source.append({})
kwargs_source_init.append({'R_sersic': .3, 'n_sersic': fix_n, 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
kwargs_source_sigma.append({'n_sersic_sigma': 1., 'R_sersic_sigma': 0.5, 'e1_sigma': 0.1, 'e2_sigma': 0.1, 'center_x_sigma': 0.1, 'center_y_sigma': 0.1})
kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -10, 'center_y': -10})
kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': 10, 'center_y': 10})
source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]

import sys
sys.path.insert(0,'../../../py_tools/')
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result
PSF1 = pyfits.getdata('PSF1.fits')
#PSF1 /= PSF1.sum()
fixcenter = False
#tag = 'test'
tag = 'fix_n_fix_re_{0}'.format(int(fix_n))
ID = 'Sim'
exp_time=2400.
background_rms = 0.0076
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im=image_noisy, psf_ave=PSF1, source_params=source_params,pix_sz = 'drz06',
                                                                 background_rms=background_rms, image_plot = True, corner_plot=False,
                                                                 flux_ratio_plot=False, deep_seed = False, fixcenter= fixcenter,
                                                                 exp_time = exp_time,tag=tag, no_MCMC=True, QSO_std=QSO_std)

result = transfer_to_result(data=image_noisy,
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=QSO_std**2,
                            cut=0, filt='F140w', fixcenter=fixcenter,ID=ID,
                            QSO_msk = '', plot_compare= True, tag=tag) 
import glob
filename = 'result_fix_n_fix_re_.txt'
if_file = glob.glob(filename)   
if if_file == []:
    f =  open(filename,'w') 
elif if_file is not []:
    f = open(filename,"r+")
    f.read()
f.write("\n#========================================")    
#f.write("\n#fix_n_to: "+repr(fix_n))
del result['QSO_amp'], result['center_x'], result['center_y'] , result['host_amp']
del result['phi_G'], result['q'], result['qso_x'], result['qso_y']
del result['host_mag'], result['amp']
f.write("\n#fit result:\n"+repr(result))
f.close()
f=open(filename,"r+")
f.read()
