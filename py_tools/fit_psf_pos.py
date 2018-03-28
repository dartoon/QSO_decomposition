#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 12:28:29 2018

@author: Dartoon

fit psf and return the ra and dec position.
"""
import sys
import astropy.io.fits as pyfits
from psfs_average import psf_ave
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from flux_profile import SB_profile, QSO_psfs_compare, profiles_compare
import time
from lenstronomy.Util.kernel_util import de_shift_kernel

def fit_psf_pos(fitted_PSF, ave_psf, psf_std = None):
    '''
    Fit the fitted_PSF and with ave_psf and return the ra and dec.
    Note that this fitting is in pixel sacle
    '''
    # data specifics need to set up based on the data situation
    background_rms = 0.04  #  background noise per pixel (Gaussian)
    exp_time = 2400.  #  exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
    numPix = len(ave_psf)  #  cutout pixel size
    deltaPix = 1  #  pixel size in arcsec (area per pixel = deltaPix**2)
    fwhm = 0.1  # full width half max of PSF (only valid when psf_type='gaussian')
    psf_type = 'PIXEL'  # 'gaussian', 'pixel', 'NONE'
    kernel_size = len(ave_psf)
    kernel = ave_psf
     
#    # lens model choicers (lenstronomy requires the instances of them, but we can keep them empty)
#    fixed_lens = [{}]
#    kwargs_lens_init = [{}]
#    kwargs_lens_sigma = [{}]
#    kwargs_lower_lens = [{}]
#    kwargs_upper_lens = [{}]
#    
#    # lens light model choices (lenstronomy requires the instances of them, but we can keep them empty)
#    fixed_lens_light = [{}]
#    kwargs_lens_light_init = [{}]
#    kwargs_lens_light_sigma = [{}]
#    kwargs_lower_lens_light = [{}]
#    kwargs_upper_lens_light = [{}]
    
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
    
    point_source_list = ['UNLENSED']
    
    ### Make simulation:
    kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
    if psf_std is not None:
        kwargs_numerics.get('psf_error_map', False)     #Turn on the PSF error map
     
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
    if psf_std is not None:
        kwargs_psf['psf_error_map'] = psf_std
    
    image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
    multi_band_list = [image_band]
    
    from lenstronomy.Workflow.fitting_sequence import FittingSequence
    fitting_seq = FittingSequence(multi_band_list, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
    
    fitting_kwargs_list = [
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 50,
             'n_iterations': 50},
    #       {'fitting_routine': 'MCMC', 'n_burn': 10, 'n_run': 20, 'walkerRatio': 50, 'mpi': False,   ##Inputs  to CosmoHammer:
    #          #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
    #        'sigma_scale': .1}
    ]
    
    lens_result, source_result, lens_light_result, ps_result, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(fitting_kwargs_list)
#    print(end_time - start_time, 'total time needed for computation')
#   print('============ CONGRATULATION, YOUR JOB WAS SUCCESSFUL ================ ')
    return ps_result[0]['ra_image'][0] , ps_result[0]['dec_image'][0]

