#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 20:23:36 2018

@author: Dartoon

A func for fitting the QSO with a input PSF. The PSF std is optional.
"""
from matplotlib.pylab import plt
import numpy as np
import corner

def fit_qso(QSO_im, psf_ave, psf_std=None, source_params=None, background_rms=0.04, pix_sz = 'drz06',
            exp_time = 2400., fix_n=None, image_plot = True, corner_plot=True,
            flux_ratio_plot=True, deep_seed = False, fixcenter = False, QSO_msk=None, QSO_std=None,
            tag = None, no_MCMC= False):
    '''
    A quick fit for the QSO image with (so far) single sersice + one PSF. The input psf noise is optional.
    
    Parameter
    --------
        QSO_im: An array of the QSO image.
        psf_ave: The psf image.
        psf_std: The psf noise, optional.
        source_params: The prior for the source. Default is given.
        background_rms: default as 0.04
        exp_time: default at 2400.
        deep_seed: if Ture, more mcmc steps will be performed.
        tag: The name tag for save the plot
            
    Return
    --------
        Will output the fitted image (Set image_plot = True), the corner_plot and the flux_ratio_plot.
        source_result, ps_result, image_ps, image_host
    
    To do
    --------
        
    '''
    # data specifics need to set up based on the data situation
    background_rms = background_rms  #  background noise per pixel (Gaussian)
    exp_time = exp_time  #  exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
    numPix = len(QSO_im)  #  cutout pixel size
    if pix_sz == 'swarp' :
        deltaPix = 0.127985  #  pixel size in arcsec (area per pixel = deltaPix**2)
    elif pix_sz == 'drz06':
        deltaPix = 0.0642
    elif pix_sz == 'acs':
        deltaPix = 0.03
    fwhm = 0.1  # full width half max of PSF (only valid when psf_type='gaussian')
    psf_type = 'PIXEL'  # 'gaussian', 'pixel', 'NONE'
    kernel_size = len(psf_ave)
    kernel = psf_ave
    
    if psf_std is not None:
        kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False, 'psf_error_map': True}     #Turn on the PSF error map
    else: 
        kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
    
    if source_params is None:
        # here are the options for the host galaxy fitting
        fixed_source = []
        kwargs_source_init = []
        kwargs_source_sigma = []
        kwargs_lower_source = []
        kwargs_upper_source = []
        
        # Disk component, as modelled by an elliptical Sersic profile
        if fix_n == None:
            fixed_source.append({})  # we fix the Sersic index to n=1 (exponential)
            kwargs_source_init.append({'R_sersic': 0.3, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
            kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
            kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -10, 'center_y': -10})
            kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': 10, 'center_y': 10})
        elif fix_n is not None:
            fixed_source.append({'n_sersic': fix_n})
            kwargs_source_init.append({'R_sersic': 0.3, 'n_sersic': fix_n, 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
            kwargs_source_sigma.append({'n_sersic': 0.001, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
            kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': fix_n, 'center_x': -10, 'center_y': -10})
            kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3, 'n_sersic': fix_n, 'center_x': 10, 'center_y': 10})
        source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
    else:
        source_params = source_params

    center_x = 0.0
    center_y = 0.0
    point_amp = QSO_im.sum()/2.
    
    fixed_ps = [{}]
    kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]
    kwargs_ps_init = kwargs_ps
    kwargs_ps_sigma = [{'ra_image': [0.01], 'dec_image': [0.01]}]
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
    data_class.update_data(QSO_im)
    
    from lenstronomy.LightModel.light_model import LightModel
    light_model_list = ['SERSIC_ELLIPSE']
    lightModel = LightModel(light_model_list=light_model_list)
    from lenstronomy.PointSource.point_source import PointSource
    point_source_list = ['UNLENSED']
    pointSource = PointSource(point_source_type_list=point_source_list)
    
    from lenstronomy.ImSim.image_model import ImageModel
    if QSO_msk is not None:
        kwargs_numerics['mask'] = QSO_msk
    
    kwargs_model = { 'source_light_model_list': light_model_list,
                    'point_source_model_list': point_source_list
                    }
    # numerical options and fitting sequences
    
    kwargs_constraints = {'joint_center_source_light': fixcenter,  # if set to True, all the components in the host galaxy will have a shared center
                          'fix_to_point_source_list': [fixcenter],  # this results in a shared center of the host galaxy with the point source (quasar)
                          'num_point_source_list': [1]
                          }
    
    kwargs_likelihood = {'check_bounds': True,  #Set the bonds, if exceed, reutrn "penalty"
                         'source_marg': False,  #In likelihood_module.LikelihoodModule -- whether to fully invert the covariance matrix for marginalization
                          'check_positive_flux': True,        
                         }
    kwargs_data = data_class.constructor_kwargs() # The "dec_at_xy_0" means the dec at the (0,0) point.
    if QSO_std is not None:
        kwargs_data['noise_map'] = QSO_std
    
    kwargs_psf = psf_class.constructor_kwargs()
    if psf_std is not None:
        kwargs_psf['psf_error_map'] = psf_std
    
    from lenstronomy.Data.psf import PSF
    psf_class = PSF(kwargs_psf)
    
    
    imageModel = ImageModel(data_class, psf_class, source_model_class=lightModel,
                                    point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
                  
    image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
    multi_band_list = [image_band]
    
    from lenstronomy.Workflow.fitting_sequence import FittingSequence
#    mpi = False  # MPI possible, but not supported through that notebook.
    # The Params for the fitting. kwargs_init: initial input. kwargs_sigma: The parameter uncertainty. kwargs_fixed: fixed parameters;
    #kwargs_lower,kwargs_upper: Lower and upper limits.

    fitting_seq = FittingSequence(multi_band_list, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
    
    if deep_seed == False:
        fitting_kwargs_list = [
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 0.8, 'n_particles': 100,
             'n_iterations': 100},
            {'fitting_routine': 'MCMC', 'n_burn': 10, 'n_run': 20, 'walkerRatio': 50, 'mpi': False,   ##Inputs  to CosmoHammer:
               #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
            'sigma_scale': .1}
            ]
    elif deep_seed == True:
         fitting_kwargs_list = [
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 150,
             'n_iterations': 150},
            {'fitting_routine': 'MCMC', 'n_burn': 50, 'n_run': 100, 'walkerRatio': 50, 'mpi': False,   ##Inputs  to CosmoHammer:
               #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
            'sigma_scale': .1}
            ]
    if no_MCMC == True:
        fitting_kwargs_list = [fitting_kwargs_list[0],
                               ]        
    
    import time
    start_time = time.time()
    lens_result, source_result, lens_light_result, ps_result, cosmo_temp, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(fitting_kwargs_list)
    end_time = time.time()
    print(end_time - start_time, 'total time needed for computation')
    print('============ CONGRATULATION, YOUR JOB WAS SUCCESSFUL ================ ')
    
    # this is the linear inversion. The kwargs will be updated afterwards
    image_reconstructed, error_map, _, _ = imageModel.image_linear_solve(kwargs_source=source_result, kwargs_ps=ps_result)
    image_ps = imageModel.point_source(ps_result)
    image_host = imageModel.source_surface_brightness(source_result)
    # let's plot the output of the PSO minimizer
    from lenstronomy.Plots.output_plots import LensModelPlot
    lensPlot = LensModelPlot(kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, lens_result, source_result,
                                 lens_light_result, ps_result, arrow_size=0.02, cmap_string="gist_heat")
    
    if image_plot:
        f, axes = plt.subplots(3, 3, figsize=(16, 16), sharex=False, sharey=False)
        lensPlot.data_plot(ax=axes[0,0])
        lensPlot.model_plot(ax=axes[0,1])
        lensPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
        
        lensPlot.decomposition_plot(ax=axes[1,0], text='Host galaxy', source_add=True, unconvolved=True)
        lensPlot.decomposition_plot(ax=axes[1,1], text='Host galaxy convolved', source_add=True)
        lensPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
        
        lensPlot.subtract_from_data_plot(ax=axes[2,0], text='Data - Point Source', point_source_add=True)
        lensPlot.subtract_from_data_plot(ax=axes[2,1], text='Data - host galaxy', source_add=True)
        lensPlot.subtract_from_data_plot(ax=axes[2,2], text='Data - host galaxy - Point Source', source_add=True, point_source_add=True)
        
        f.tight_layout()
        #f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
        if tag is not None:
            f.savefig('{0}_fitted_image.pdf'.format(tag))
        plt.show()
        
    if corner_plot==True and no_MCMC==False:
        # here the (non-converged) MCMC chain of the non-linear parameters
        if not samples_mcmc == []:
           n, num_param = np.shape(samples_mcmc)
           plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
           if tag is not None:
               plot.savefig('{0}_para_corner.pdf'.format(tag))
           plt.show()
           
        
    if flux_ratio_plot==True and no_MCMC==False:
        from lenstronomy.Workflow.parameters import Param
        param = Param(kwargs_model, kwargs_constraints, kwargs_fixed_source=source_params[2], kwargs_fixed_ps=fixed_ps)
        mcmc_new_list = []
        labels_new = [r"Quasar flux", r"host_flux", r"source_x", r"source_y"]
        
        # transform the parameter position of the MCMC chain in a lenstronomy convention with keyword arguments #
        for i in range(len(samples_mcmc)/10):
            kwargs_lens_out, kwargs_light_source_out, kwargs_light_lens_out, kwargs_ps_out, kwargs_cosmo = param.getParams(samples_mcmc[i+ len(samples_mcmc)/10*9])
            image_reconstructed, _, _, _ = imageModel.image_linear_solve(kwargs_source=kwargs_light_source_out, kwargs_ps=kwargs_ps_out)
            
            image_ps = imageModel.point_source(kwargs_ps_out)
            flux_quasar = np.sum(image_ps)
            image_disk = imageModel.source_surface_brightness(kwargs_light_source_out, k=0)
            flux_disk = np.sum(image_disk)
            source_x = kwargs_ps_out[0]['ra_image']
            source_y = kwargs_ps_out[0]['dec_image']
            #    image_buldge = imageModel.source_surface_brightness(kwargs_light_source_out, k=1)
            #    flux_buldge = np.sum(image_buldge)
#            kwargs_ps_out
            if flux_disk>0:
                mcmc_new_list.append([flux_quasar, flux_disk, source_x, source_y])
        plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
        if tag is not None:
            plot.savefig('{0}_HOSTvsQSO_corner.pdf'.format(tag))
        plt.show()
    if QSO_std is None:
        return source_result, ps_result, image_ps, image_host, np.sqrt(data_class.C_D+np.abs(error_map))
    else:
        return source_result, ps_result, image_ps, image_host, np.sqrt(QSO_std**2+np.abs(error_map))

def fit_ps(QSO_im, psf_ave, psf_std=None, background_rms=0.04, source_params=None, pix_sz = None,
            exp_time = 2400., fix_n=None, image_plot = True, corner_plot=True,
            flux_ratio_plot=True, deep_seed = False, fixcenter = True, QSO_msk=None, QSO_std=None,
            tag = None, no_MCMC= False):
    '''
    Fit the fitted_PSF and with ave_psf and return the ra and dec.
    Note that this fitting is in pixel sacle
    '''
    # data specifics need to set up based on the data situation
    numPix = len(QSO_im)  #  cutout pixel size
    if pix_sz is not 'acs':
        deltaPix = 0.0642  #  pixel size
    else:
        deltaPix = 0.03
    fwhm = 0.1  # full width half max of PSF (only valid when psf_type='gaussian')
    psf_type = 'PIXEL'  # 'gaussian', 'pixel', 'NONE'
    kernel_size = len(psf_ave)
    kernel = psf_ave
     
    # here are the options for the host galaxy fitting
    fixed_source = []
    kwargs_source_init = []
    kwargs_source_sigma = []
    kwargs_lower_source = []
    kwargs_upper_source = []
    
    source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
    center_x = 0.0
    center_y = 0.0
    point_amp = QSO_im.sum()
    
    fixed_ps = [{}]
    kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]
    kwargs_ps_init = kwargs_ps
    kwargs_ps_sigma = [{'ra_image': [0.01], 'dec_image': [0.01]}]
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
    data_class.update_data(QSO_im)
    from lenstronomy.PointSource.point_source import PointSource
    point_source_list = ['UNLENSED']
    pointSource = PointSource(point_source_type_list=point_source_list)
    ### Make simulation:
    if psf_std is not None:
        kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False, 'psf_error_map': True}     #Turn on the PSF error map
    else: 
        kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
    if QSO_msk is not None:
        kwargs_numerics['mask'] = QSO_msk

    from lenstronomy.ImSim.image_model import ImageModel
    imageModel = ImageModel(data_class, psf_class,
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
    if psf_std is not None:
        kwargs_psf['psf_error_map'] = psf_std
    from lenstronomy.Data.psf import PSF
    psf_class = PSF(kwargs_psf)
    
    image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
    multi_band_list = [image_band]
    
    from lenstronomy.Workflow.fitting_sequence import FittingSequence
    fitting_seq = FittingSequence(multi_band_list, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
    
    fitting_kwargs_list = [
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 50,
             'n_iterations': 50},
    ]
    
    import time
    start_time = time.time()
    lens_result, source_result, lens_light_result, ps_result, cosmo_temp, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc  = fitting_seq.fit_sequence(fitting_kwargs_list)
    end_time = time.time()
    print(end_time - start_time, 'total time needed for computation')
    print('============ CONGRATULATION, YOUR JOB WAS SUCCESSFUL ================ ')
    # this is the linear inversion. The kwargs will be updated afterwards
    image_reconstructed, _, _, _ = imageModel.image_linear_solve(kwargs_source=source_result, kwargs_ps=ps_result)
    image_ps = imageModel.point_source(ps_result)
    # let's plot the output of the PSO minimizer
    if image_plot:
        from lenstronomy.Plots.output_plots import LensModelPlot
        lensPlot = LensModelPlot(kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, lens_result, source_result,
                                 lens_light_result, ps_result, arrow_size=0.02, cmap_string="gist_heat")
    
        f, axes = plt.subplots(2, 2, figsize=(9, 9), sharex=False, sharey=False)
        lensPlot.data_plot(ax=axes[0,0])
        lensPlot.model_plot(ax=axes[0,1])
        lensPlot.normalized_residual_plot(ax=axes[1,0], v_min=-6, v_max=6)
        lensPlot.subtract_from_data_plot(ax=axes[1,1], text='Data - Point Source', point_source_add=True)
        f.tight_layout()
        #f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
        if tag is not None:
            f.savefig('{0}_fitted_image.pdf'.format(tag))
        plt.show()    
#    return source_result, ps_result, image_ps, QSO_im-image_ps, data_class.C_D
#
    if QSO_std is None:
        return source_result, ps_result, image_ps, QSO_im-image_ps, np.sqrt(data_class.C_D)
    else:
        return source_result, ps_result, image_ps, QSO_im-image_ps, QSO_std


"""
def fit_qso_disk_buldge(QSO_im, psf_ave, psf_std=None, source_params=None, background_rms=0.04, pix_sz = 'swarp',
            exp_time = 2400., fix_n=None, image_plot = True, corner_plot=True,
            flux_ratio_plot=True, deep_seed = False, fixcenter = True, QSO_msk=None,QSO_std=None):
    '''
    A quick fit for the QSO image with (so far) single sersice + one PSF. The input psf noise is optional.
    
    Parameter
    --------
        QSO_im: An array of the QSO image.
        psf_ave: The psf image.
        psf_std: The psf noise, optional.
        source_params: The prior for the source. Default is given.
        background_rms: default as 0.04
        exp_time: default at 2400.
        deep_seed: if Ture, more mcmc steps will be performed.
            
    Return
    --------
        Will output the fitted image (Set image_plot = True), the corner_plot and the flux_ratio_plot.
        source_result, ps_result, image_ps, image_host
    
    To do
    --------
        
    '''
    # data specifics need to set up based on the data situation
    background_rms = background_rms  #  background noise per pixel (Gaussian)
    exp_time = exp_time  #  exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
    numPix = len(QSO_im)  #  cutout pixel size
    if pix_sz == 'swarp' :
        deltaPix = 0.127985  #  pixel size in arcsec (area per pixel = deltaPix**2)
    elif pix_sz == 'drz06':
        deltaPix = 0.0642
    elif pix_sz == 'acs':
        deltaPix = 0.03
    fwhm = 0.1  # full width half max of PSF (only valid when psf_type='gaussian')
    psf_type = 'PIXEL'  # 'gaussian', 'pixel', 'NONE'
    kernel_size = len(psf_ave)
    kernel = psf_ave
    
    kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
    if psf_std is not None:
        kwargs_numerics.get('psf_error_map', True)     #Turn on the PSF error map
    
    if source_params is None:
        # here are the options for the host galaxy fitting
        fixed_source = []
        kwargs_source_init = []
        kwargs_source_sigma = []
        kwargs_lower_source = []
        kwargs_upper_source = []
        
#         Disk component, as modelled by an elliptical Sersic profile
        fixed_source.append({'n_sersic': 1})  # we fix the Sersic index to n=1 (exponential)
        kwargs_source_init.append({'R_sersic': 1., 'n_sersic': 1, 'e1': 0, 'e2': 0, 'center_x': 0, 'center_y': 0})
        kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': .3, 'center_x': -10, 'center_y': -10})
        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3, 'n_sersic': 5., 'center_x': 10, 'center_y': 10})
#         Buldge component, as modelled by a spherical Sersic profile
        fixed_source.append({'n_sersic': 4})  # we fix the Sersic index to n=4 (buldgy)
        kwargs_source_init.append({'R_sersic': .3, 'n_sersic': 4, 'e1': 0, 'e2': 0, 'center_x': 0, 'center_y': 0})
        kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.3, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': .5, 'center_x': -10, 'center_y': -10})
        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3, 'n_sersic': 5., 'center_x': 10, 'center_y': 10})
        source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
    else:
        source_params = source_params

    center_x = 0.0
    center_y = 0.0
    point_amp = QSO_im.sum()/2.
    
    fixed_ps = [{}]
    kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]
    kwargs_ps_init = kwargs_ps
    kwargs_ps_sigma = [{'ra_image': [0.01], 'dec_image': [0.01]}]
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
    data_class.update_data(QSO_im)
    
    from lenstronomy.LightModel.light_model import LightModel
    light_model_list = ['SERSIC_ELLIPSE','SERSIC_ELLIPSE']
    lightModel = LightModel(light_model_list=light_model_list)
    from lenstronomy.PointSource.point_source import PointSource
    point_source_list = ['UNLENSED']
    pointSource = PointSource(point_source_type_list=point_source_list)
    
    ### Make simulation:
    from lenstronomy.ImSim.image_model import ImageModel
    if QSO_msk is None:
        kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
    else:
        kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False, 'mask': QSO_msk}
    imageModel = ImageModel(data_class, psf_class, source_model_class=lightModel,
                                    point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
    
    kwargs_model = { 'source_light_model_list': light_model_list,
                    'point_source_model_list': point_source_list
                    }
    # numerical options and fitting sequences
    
    kwargs_constraints = {'joint_center_source_light': fixcenter,  # if set to True, all the components in the host galaxy will have a shared center
                          'fix_to_point_source_list': [fixcenter,fixcenter],  # this results in a shared center of the host galaxy with the point source (quasar)
                          'num_point_source_list': [1]
                          }
    
    kwargs_likelihood = {'check_bounds': True,  #Set the bonds, if exceed, reutrn "penalty"
                         'source_marg': False,  #In likelihood_module.LikelihoodModule -- whether to fully invert the covariance matrix for marginalization
                                 }
    kwargs_data = data_class.constructor_kwargs() # The "dec_at_xy_0" means the dec at the (0,0) point.
    if QSO_std is not None:
        kwargs_data['noise_map'] = QSO_std
        
    kwargs_psf = psf_class.constructor_kwargs()
    if psf_std is not None:
        kwargs_psf['psf_error_map'] = psf_std
                  
    image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
    multi_band_list = [image_band]
    
    from lenstronomy.Workflow.fitting_sequence import FittingSequence
#    mpi = False  # MPI possible, but not supported through that notebook.
    # The Params for the fitting. kwargs_init: initial input. kwargs_sigma: The parameter uncertainty. kwargs_fixed: fixed parameters;
    #kwargs_lower,kwargs_upper: Lower and upper limits.

    fitting_seq = FittingSequence(multi_band_list, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
    
    if deep_seed == False:
        fitting_kwargs_list = [
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 50,
             'n_iterations': 50},
            {'fitting_routine': 'MCMC', 'n_burn': 10, 'n_run': 20, 'walkerRatio': 50, 'mpi': False,   ##Inputs  to CosmoHammer:
               #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
            'sigma_scale': .1}
            ]
    elif deep_seed == True:
         fitting_kwargs_list = [
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 50,
             'n_iterations': 50},
            {'fitting_routine': 'MCMC', 'n_burn': 30, 'n_run': 30, 'walkerRatio': 100, 'mpi': False,   ##Inputs  to CosmoHammer:
               #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
            'sigma_scale': .1}
            ]
    
    import time
    start_time = time.time()
    lens_result, source_result, lens_light_result, ps_result, cosmo_temp, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(fitting_kwargs_list)
    end_time = time.time()
    print(end_time - start_time, 'total time needed for computation')
    print('============ CONGRATULATION, YOUR JOB WAS SUCCESSFUL ================ ')
    
    # this is the linear inversion. The kwargs will be updated afterwards
    image_reconstructed, _, _, _ = imageModel.image_linear_solve(kwargs_source=source_result, kwargs_ps=ps_result)
    image_ps = imageModel.point_source(ps_result)
    image_host = imageModel.source_surface_brightness(source_result)
    # let's plot the output of the PSO minimizer
    if image_plot:
        from lenstronomy.Plots.output_plots import LensModelPlot
        lensPlot = LensModelPlot(kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, lens_result, source_result,
                                 lens_light_result, ps_result, arrow_size=0.02, cmap_string="gist_heat")
    
        f, axes = plt.subplots(3, 3, figsize=(16, 16), sharex=False, sharey=False)
        lensPlot.data_plot(ax=axes[0,0])
        lensPlot.model_plot(ax=axes[0,1])
        lensPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
        
        lensPlot.decomposition_plot(ax=axes[1,0], text='Source light', source_add=True, unconvolved=True)
        lensPlot.decomposition_plot(ax=axes[1,1], text='Source light convolved', source_add=True)
        lensPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
        
        lensPlot.subtract_from_data_plot(ax=axes[2,0], text='Data - Point Source', point_source_add=True)
        lensPlot.subtract_from_data_plot(ax=axes[2,1], text='Data - host galaxy', source_add=True)
        lensPlot.subtract_from_data_plot(ax=axes[2,2], text='Data - host galaxy - Point Source', source_add=True, point_source_add=True)
        
        f.tight_layout()
        #f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
        plt.show()
        
    if corner_plot:
        # here the (non-converged) MCMC chain of the non-linear parameters
        if not samples_mcmc == []:
           n, num_param = np.shape(samples_mcmc)
           plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
        
    if flux_ratio_plot:
        from lenstronomy.Workflow.parameters import Param
        param = Param(kwargs_model, kwargs_constraints, kwargs_fixed_source=fixed_source, kwargs_fixed_ps=fixed_ps)
        mcmc_new_list = []
        labels_new = [r"Quasar flux", r"host_disk",r"host_bluge", r"source_x", r"source_y"]
        
        # transform the parameter position of the MCMC chain in a lenstronomy convention with keyword arguments #
        for i in range(len(samples_mcmc)):
            kwargs_lens_out, kwargs_light_source_out, kwargs_light_lens_out, kwargs_ps_out, kwargs_cosmo = param.getParams(samples_mcmc[i])
            image_reconstructed, _, _, _ = imageModel.image_linear_solve(kwargs_source=kwargs_light_source_out, kwargs_ps=kwargs_ps_out)
            
            image_ps = imageModel.point_source(kwargs_ps_out)
            flux_quasar = np.sum(image_ps)
            image_disk = imageModel.source_surface_brightness(kwargs_light_source_out, k=0)
            flux_disk = np.sum(image_disk)
            image_buldge = imageModel.source_surface_brightness(kwargs_light_source_out, k=1)
            flux_buldge = np.sum(image_buldge)
#            kwargs_ps_out
            source_x = kwargs_ps_out[0]['ra_image']
            source_y = kwargs_ps_out[0]['dec_image']
            if flux_disk>0:
                mcmc_new_list.append([flux_quasar, flux_disk, flux_buldge, source_x, source_y])
        plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
    plt.show()
    return source_result, ps_result, image_ps, image_host, data_class.C_D
"""


def fit_single_host(QSO_im, psf_ave, background_rms=0.04, pix_sz = 'swarp',
            exp_time = 2400., QSO_msk=None, QSO_std=None, tag=None):
    '''
    A quick fit for the QSO image with (so far) single sersice + one PSF. The input psf noise is optional.
    
    Parameter
    --------
        QSO_im: An array of the QSO image.
        psf_ave: The psf image.
        psf_std: The psf noise, optional.
        source_params: The prior for the source. Default is given.
        background_rms: default as 0.04
        exp_time: default at 2400.
        deep_seed: if Ture, more mcmc steps will be performed.
            
    Return
    --------
        Will output the fitted image (Set image_plot = True), the corner_plot and the flux_ratio_plot.
        source_result, ps_result, image_ps, image_host
    
    To do
    --------
        
    '''
    # data specifics need to set up based on the data situation
    background_rms = background_rms  #  background noise per pixel (Gaussian)
    exp_time = exp_time  #  exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
    numPix = len(QSO_im)  #  cutout pixel size
    if pix_sz == 'swarp' :
        deltaPix = 0.127985  #  pixel size in arcsec (area per pixel = deltaPix**2)
    elif pix_sz == 'drz06':
        deltaPix = 0.0642
    elif pix_sz == 'acs':
        deltaPix = 0.03
    fwhm = 0.1  # full width half max of PSF (only valid when psf_type='gaussian')
    psf_type = 'PIXEL'  # 'gaussian', 'pixel', 'NONE'
    kernel_size = len(psf_ave)
    kernel = psf_ave
    
    kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
    
    #==============================================================================
    #Doing the QSO fitting 
    #==============================================================================
    from lenstronomy.SimulationAPI.simulations import Simulation
    SimAPI = Simulation()
    data_class = SimAPI.data_configure(numPix, deltaPix, exp_time, background_rms)
    psf_class = SimAPI.psf_configure(psf_type=psf_type, fwhm=fwhm, kernelsize=kernel_size, deltaPix=deltaPix, kernel=kernel)
    
    data_class.update_data(QSO_im)
    
    # lens light model choices
    fixed_lens_light = []
    kwargs_lens_light_init = []
    kwargs_lens_light_sigma = []
    kwargs_lower_lens_light = []
    kwargs_upper_lens_light = []
    
    fixed_lens_light.append({})
    kwargs_lens_light_init.append({'R_sersic': 1., 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
    kwargs_lens_light_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
    kwargs_lower_lens_light.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -10, 'center_y': -10})
    kwargs_upper_lens_light.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': 10, 'center_y': 10})
        
    lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma,
                         fixed_lens_light, kwargs_lower_lens_light,
                         kwargs_upper_lens_light]
    kwargs_params = {'lens_light_model': lens_light_params}
    
    ### Make simulation:
    if QSO_msk is None:
        kwargs_numerics = {'subgrid_res': 1}
    else:
        kwargs_numerics = {'subgrid_res': 1, 'mask': QSO_msk}
        
    
    lens_light_model_list = ['SERSIC_ELLIPSE']
    
    kwargs_model = {'lens_light_model_list': lens_light_model_list}
    kwargs_constraints = {'joint_center_lens_light': False}
    kwargs_numerics_galfit = {'subgrid_res': 1}
    kwargs_likelihood = {'check_bounds': True}
    
    kwargs_data = data_class.constructor_kwargs()
    kwargs_psf = psf_class.constructor_kwargs()

    from lenstronomy.LightModel.light_model import LightModel
    lightModel = LightModel(lens_light_model_list)
    from lenstronomy.ImSim.image_model import ImageModel
    imageModel = ImageModel(data_class, psf_class, lens_light_model_class=lightModel, kwargs_numerics=kwargs_numerics)
    
    
    if QSO_std is not None:
        kwargs_data['noise_map'] = QSO_std
    
    image_band = [kwargs_data, kwargs_psf, kwargs_numerics_galfit]
    multi_band_list = [image_band]
    
    from lenstronomy.Workflow.fitting_sequence import FittingSequence
    fitting_seq = FittingSequence(multi_band_list, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
    
    fitting_kwargs_list = [
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 50,
             'n_iterations': 50},
#            {'fitting_routine': 'MCMC', 'n_burn': 10, 'n_run': 20, 'walkerRatio': 50, 'mpi': False,   ##Inputs  to CosmoHammer:
#               #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
#            'sigma_scale': .1}
            ]
    lens_result, source_result, lens_light_result, ps_result, cosmo_result, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(fitting_kwargs_list)
    
    from lenstronomy.Plots.output_plots import LensModelPlot
    import lenstronomy.Plots.output_plots as out_plot
    lensPlot = LensModelPlot(kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, lens_result, source_result,
                             lens_light_result, ps_result, arrow_size=0.02, cmap_string="gist_heat")
    
    for i in range(len(chain_list)):
        if len(param_list[i]) > 0:
            f, axes = out_plot.plot_chain(chain_list[i], param_list[i])
    
    f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    lensPlot.data_plot(ax=axes[0,0])
    lensPlot.model_plot(ax=axes[0,1])
    lensPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
    #f.tight_layout()
    #f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    #plt.show()
    #f, axes = plt.subplots(1, 3, figsize=(16, 8), sharex=False, sharey=False)
    
    lensPlot.decomposition_plot(ax=axes[1,0], text='Lens light', lens_light_add=True, unconvolved=True)
    lensPlot.decomposition_plot(ax=axes[1,1], text='Lens light convolved', lens_light_add=True)
    lensPlot.subtract_from_data_plot(ax=axes[1,2], text='Data - Lens Light', lens_light_add=True)
    f.tight_layout()
    #f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    plt.show()
    if tag is not None:
            f.savefig('{0}_fitted_image.pdf'.format(tag))
    image_reconstructed, _, _, _ = imageModel.image_linear_solve(kwargs_lens_light=lens_light_result)
    
    image_host = imageModel.lens_surface_brightness(lens_light_result)
#    print lens_light_result
    if QSO_std is None:
        return lens_light_result, image_host, np.sqrt(data_class.C_D)
    else:
        return lens_light_result, image_host, QSO_std
