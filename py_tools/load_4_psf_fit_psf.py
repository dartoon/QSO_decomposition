#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 11:48:07 2018

@author: Dartoon
"""

import sys
sys.path.insert(0,'../py_tools')
import astropy.io.fits as pyfits
from psfs_average import psf_ave
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from flux_profile import SB_profile, QSO_psfs_compare

# data specifics need to set up based on the data situation
background_rms = 0.04  #  background noise per pixel (Gaussian)
exp_time = 2400.  #  exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
numPix = len(psf_init_ave)  #  cutout pixel size
deltaPix = 1  #  pixel size in arcsec (area per pixel = deltaPix**2)
fwhm = 0.1  # full width half max of PSF (only valid when psf_type='gaussian')
psf_type = 'PIXEL'  # 'gaussian', 'pixel', 'NONE'
kernel_size = len(psf_init_ave)
kernel = psf_init_ave

# lens model choicers (lenstronomy requires the instances of them, but we can keep them empty)
fixed_lens = [{}]
kwargs_lens_init = [{}]
kwargs_lens_sigma = [{}]
kwargs_lower_lens = [{}]
kwargs_upper_lens = [{}]

# lens light model choices (lenstronomy requires the instances of them, but we can keep them empty)
fixed_lens_light = [{}]
kwargs_lens_light_init = [{}]
kwargs_lens_light_sigma = [{}]
kwargs_lower_lens_light = [{}]
kwargs_upper_lens_light = [{}]

# here are the options for the host galaxy fitting
fixed_source = []
kwargs_source_init = []
kwargs_source_sigma = []
kwargs_lower_source = []
kwargs_upper_source = []

source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
center_x = 0.0
center_y = 0.0
point_amp = fitted_PSF.sum()

fixed_ps = [{}]
kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]
kwargs_ps_init = kwargs_ps
kwargs_ps_sigma = [{'pos_sigma': 0.01, 'pos_sigma': 0.01}]
kwargs_lower_ps = [{'ra_image': [-10], 'dec_image': [-10]}]
kwargs_upper_ps = [{'ra_image': [10], 'dec_image': [10]}]
ps_param = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]

kwargs_params = {'source_model': source_params,
               'point_source_model': ps_param}

#==============================================================================
#Doing the QSO fitting 
#==============================================================================
from lenstronomy.SimulationAPI.simulations import Simulation
SimAPI = Simulation()
data_class = SimAPI.data_configure(numPix, deltaPix, exp_time, background_rms)
psf_class = SimAPI.psf_configure(psf_type=psf_type, fwhm=fwhm, kernelsize=kernel_size, deltaPix=deltaPix, kernel=kernel)
data_class.update_data(fitted_PSF)

from lenstronomy.PointSource.point_source import PointSource
point_source_list = ['UNLENSED']
pointSource = PointSource(point_source_type_list=point_source_list)

### Make simulation:
from lenstronomy.ImSim.image_model import ImageModel
kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
imageModel = ImageModel(data_class, psf_class, #source_model_class=lightModel,
                                point_source_class=pointSource, kwargs_numerics=kwargs_numerics)

kwargs_model = { #'source_light_model_list': light_model_list,
                'point_source_model_list': point_source_list
                 }

kwargs_constraints = {'joint_center_source_light': True,  # if set to True, all the components in the host galaxy will have a shared center
                      'fix_to_point_source_list': [True, True],  # this results in a shared center of the host galaxy with the point source (quasar)
                      'num_point_source_list': [1]
                     }

kwargs_likelihood = {'check_bounds': True,  #Set the bonds, if exceed, reutrn "penalty"
                     'source_marg': False,  #In likelihood_module.LikelihoodModule -- whether to fully invert the covariance matrix for marginalization
                             }
kwargs_data = data_class.constructor_kwargs() # The "dec_at_xy_0" means the dec at the (0,0) point.
kwargs_psf = psf_class.constructor_kwargs()
kwargs_psf['psf_error_map'] = psf_std

image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
multi_band_list = [image_band]

from lenstronomy.Workflow.fitting_sequence import FittingSequence

mpi = False  

fitting_seq = FittingSequence(multi_band_list, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)

fitting_kwargs_list = [
        {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 50,
         'n_iterations': 50},
#        {'fitting_routine': 'MCMC', 'n_burn': 10, 'n_run': 20, 'walkerRatio': 50, 'mpi': False,   ##Inputs  to CosmoHammer:
#            #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
#         'sigma_scale': .1}
]