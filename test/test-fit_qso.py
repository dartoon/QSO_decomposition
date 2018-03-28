#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 20:30:07 2018

@author: Dartoon

Test fit_qso
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pylab as plt
import glob
import sys
sys.path.insert(0,'../py_tools')
from psfs_average import psf_ave
from matplotlib.colors import LogNorm

ID = 'CID1174'
folder = '../fits_image/{0}/analysis/'.format(ID)

psf_name_list = glob.glob(folder+"PSF*.fits")   # Read *.reg files in a list.

psf_list = []
for i in range(len(psf_name_list)):
    psf_get = pyfits.getdata(folder+'PSF{0}.fits'.format(i))
    psf_list.append(psf_get)

mask_list = glob.glob(folder+"PSF*.reg")   # Read *.reg files in a list.
psf_ave, psf_std=psf_ave(psf_list,mode = 'CI', not_count=(1,4,5),
                  mask_list=mask_list)
plt.imshow(psf_ave, origin='low', norm = LogNorm())
plt.show()
QSO_im = pyfits.getdata(folder+'{0}_cutout.fits'.format(ID))

from fit_qso import fit_qso
## here are the options for the host galaxy fitting
#fixed_source = []
#kwargs_source_init = []
#kwargs_source_sigma = []
#kwargs_lower_source = []
#kwargs_upper_source = []
#
## Disk component, as modelled by an elliptical Sersic profile
#fixed_source.append({})  # we fix the Sersic index to n=1 (exponential)
#kwargs_source_init.append({'R_sersic': 1., 'n_sersic': 1, 'e1': 0, 'e2': 0, 'center_x': 0, 'center_y': 0})
#kwargs_source_sigma.append({'n_sersic_sigma': 0.5, 'R_sersic_sigma': 0.5, 'ellipse_sigma': 0.1, 'center_x_sigma': 0.1, 'center_y_sigma': 0.1})
#kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.001, 'n_sersic': .5, 'center_x': -10, 'center_y': -10})
#kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10, 'n_sersic': 5., 'center_x': 10, 'center_y': 10})
#
### Buldge component, as modelled by a spherical Sersic profile
##fixed_source.append({'n_sersic': 4})  # we fix the Sersic index to n=4 (buldgy)
##kwargs_source_init.append({'R_sersic': .5, 'n_sersic': 4, 'center_x': 0, 'center_y': 0})
##kwargs_source_sigma.append({'n_sersic_sigma': 0.5, 'R_sersic_sigma': 0.3, 'center_x_sigma': 0.1, 'center_y_sigma': 0.1})
##kwargs_lower_source.append({'R_sersic': 0.001, 'n_sersic': .5, 'center_x': -10, 'center_y': -10})
##kwargs_upper_source.append({'R_sersic': 10, 'n_sersic': 5., 'center_x': 10, 'center_y': 10})
#source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]

source_result, ps_result, image_ps, image_host=fit_qso(QSO_im, psf_ave, psf_std=psf_std,
                                                       source_params=None, image_plot = False, corner_plot=False, flux_ratio_plot=True)

#==============================================================================
#Plot the images for adopting in the paper
#==============================================================================
from flux_profile import total_compare
data = QSO_im
QSO = image_ps
host = image_host
flux_list = [data, QSO, host]
label = ['data', 'QSO', 'host', 'model', 'residual']
import glob
mask_list = glob.glob("files/QSO*.reg")   # Read *.reg files in a list.
total_compare(label_list = label, flux_list = flux_list, target_ID = ID, data_mask_list = mask_list)

