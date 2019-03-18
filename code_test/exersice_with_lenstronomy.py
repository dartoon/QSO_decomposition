#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 16:56:38 2018

@author: Dartoon

Exercise with Simon's notebook
"""

import numpy as np
import os
import time
import copy
import corner
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm

from lenstronomy.SimulationAPI.simulations import Simulation
SimAPI = Simulation()

# import PSF file
path = os.getcwd()
kernel = pyfits.getdata('psf.fits')

plt.matshow(kernel, norm=LogNorm(), origin='lower')
plt.show()

# data specifics
background_rms = 0.1  #  background noise per pixel (Gaussian)
exp_time = 100.  #  exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
numPix = 81  #  cutout pixel size
deltaPix = 0.05  #  pixel size in arcsec (area per pixel = deltaPix**2)
fwhm = 0.1  # full width half max of PSF (only valid when psf_type='gaussian')
psf_type = 'PIXEL'  # 'gaussian', 'pixel', 'NONE'
kernel_size = 91

# initial input simulation
# generate the coordinate grid and image properties
data_class = SimAPI.data_configure(numPix, deltaPix, exp_time, background_rms)
# generate the psf variables
psf_class = SimAPI.psf_configure(psf_type=psf_type, fwhm=fwhm, kernelsize=kernel_size, deltaPix=deltaPix, kernel=kernel)

# quasar center (we chose a off-centered position)
center_x = 0.11
center_y = 0.01

# quasar brightness (as measured as the sum of pixel values)
point_amp = 10000 
from lenstronomy.PointSource.point_source import PointSource
point_source_list = ['UNLENSED']
pointSource = PointSource(point_source_type_list=point_source_list)
kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]

from lenstronomy.LightModel.light_model import LightModel
light_model_list = ['SERSIC_ELLIPSE', 'SERSIC']
lightModel = LightModel(light_model_list=light_model_list)
import lenstronomy.Util.param_util as param_util
e1, e2 = param_util.phi_q2_ellipticity(phi=0.3, q=0.6)
kwargs_disk = {'I0_sersic': 1, 'n_sersic': 1, 'R_sersic': 0.7, 'e1': e1, 'e2': e2, 'center_x': center_x, 'center_y': center_y}
kwargs_buldge = {'I0_sersic': 1, 'n_sersic': 4, 'R_sersic': 0.3, 'center_x': center_x, 'center_y': center_y}
kwargs_host = [kwargs_disk, kwargs_buldge]


### Make simulation:
from lenstronomy.ImSim.image_model import ImageModel
kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
imageModel = ImageModel(data_class, psf_class, source_model_class=lightModel,
                                point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
# simulate image with the parameters we have defined above #
image = imageModel.image(kwargs_source=kwargs_host, kwargs_ps=kwargs_ps)
# we can also add noise #
import lenstronomy.Util.image_util as image_util
poisson = image_util.add_poisson(image, exp_time=exp_time)
bkg = image_util.add_background(image, sigma_bkd=background_rms)
image_noisy = image + bkg + poisson
plt.matshow(image_noisy, norm=LogNorm(), origin='lower')
plt.colorbar()
plt.show()
data_class.update_data(image_noisy)

'''
# we can also simulate the different components separately
imageModel_ps = ImageModel(data_class, psf_class, point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
image_ps = imageModel_ps.image(kwargs_ps=kwargs_ps)
plt.matshow(np.log10(image_ps), origin='lower')
plt.show()

imageModel_host = ImageModel(data_class, psf_class, source_model_class=lightModel, kwargs_numerics=kwargs_numerics)
image_host = imageModel_host.image(kwargs_source=kwargs_host)
plt.matshow(np.log10(image_host), origin='lower')
plt.show()
'''
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

# Disk component, as modelled by an elliptical Sersic profile
fixed_source.append({'n_sersic': 1})  # we fix the Sersic index to n=1 (exponential)
kwargs_source_init.append({'R_sersic': 1., 'n_sersic': 1, 'e1': 0, 'e2': 0, 'center_x': 0, 'center_y': 0})
kwargs_source_sigma.append({'n_sersic_sigma': 0.5, 'R_sersic_sigma': 0.5, 'ellipse_sigma': 0.1, 'center_x_sigma': 0.1, 'center_y_sigma': 0.1})
kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.001, 'n_sersic': .5, 'center_x': -10, 'center_y': -10})
kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10, 'n_sersic': 5., 'center_x': 10, 'center_y': 10})

# Buldge component, as modelled by a spherical Sersic profile
fixed_source.append({'n_sersic': 4})  # we fix the Sersic index to n=4 (buldgy)
kwargs_source_init.append({'R_sersic': .5, 'n_sersic': 4, 'center_x': 0, 'center_y': 0})
kwargs_source_sigma.append({'n_sersic_sigma': 0.5, 'R_sersic_sigma': 0.3, 'center_x_sigma': 0.1, 'center_y_sigma': 0.1})
kwargs_lower_source.append({'R_sersic': 0.001, 'n_sersic': .5, 'center_x': -10, 'center_y': -10})
kwargs_upper_source.append({'R_sersic': 10, 'n_sersic': 5., 'center_x': 10, 'center_y': 10})

source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]

fixed_ps = [{}]
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
kwargs_model = { 'source_light_model_list': light_model_list,
                'point_source_model_list': point_source_list
                 }

# numerical options and fitting sequences

kwargs_constraints = {'joint_center_source_light': True,  # if set to True, all the components in the host galaxy will have a shared center
                      'fix_to_point_source_list': [True, True],  # this results in a shared center of the host galaxy with the point source (quasar)
                      'num_point_source_list': [1]
                     }

kwargs_likelihood = {'check_bounds': True,  #Set the bonds, if exceed, reutrn "penalty"
                     'source_marg': False,  #In likelihood_module.LikelihoodModule -- whether to fully invert the covariance matrix for marginalization
                             }
kwargs_data = data_class.constructor_kwargs() # The "dec_at_xy_0" means the dec at the (0,0) point.
kwargs_psf = psf_class.constructor_kwargs()
image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
multi_band_list = [image_band]

from lenstronomy.Workflow.fitting_sequence import FittingSequence

mpi = False  # MPI possible, but not supported through that notebook.

# The Params for the fitting. kwargs_init: initial input. kwargs_sigma: The parameter uncertainty. kwargs_fixed: fixed parameters;
#kwargs_lower,kwargs_upper: Lower and upper limits.

fitting_seq = FittingSequence(multi_band_list, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
'''
Set up the FittingSequence:
    1. multi_band_list: A 'two-folder' list including the image information dicts:
        [kwargs_data, kwargs_psf, kwargs_numerics]
    2. kwargs_model: model list including {'point_source_model_list', 'source_light_model_list'}.
    3. kwargs_constraints: a constraints prior of the source position.
    4. kwargs_likelihood. As noted above
    5. kwargs_params:
        A list of: [kwargs_init, kwargs_sigma, kwargs_fixed, kwargs_lower, kwargs_upper]. 

'''


fitting_kwargs_list = [
        {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 50,
         'n_iterations': 50},
        {'fitting_routine': 'MCMC', 'n_burn': 10, 'n_run': 10, 'walkerRatio': 10, 'mpi': False,   ##Inputs  to CosmoHammer:
            #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
         'sigma_scale': .1}
]

start_time = time.time()
lens_result, source_result, lens_light_result, ps_result, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(fitting_kwargs_list)
end_time = time.time()
print(end_time - start_time, 'total time needed for computation')
print('============ CONGRATULATION, YOUR JOB WAS SUCCESSFUL ================ ')

# let's plot the output of the PSO minimizer
from lenstronomy.Plots.output_plots import LensModelPlot
lensPlot = LensModelPlot(kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, lens_result, source_result,
                             lens_light_result, ps_result, arrow_size=0.02, cmap_string="gist_heat", high_res=5)
    
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

# here the (non-converged) MCMC chain of the non-linear parameters
if not samples_mcmc == []:
        n, num_param = np.shape(samples_mcmc)
        plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
        
# this is the linear inversion. The kwargs will be updated afterwards
image_reconstructed, _, _, _ = imageModel.image_linear_solve(kwargs_source=source_result, kwargs_ps=ps_result)

# flux count in point source
image_ps = imageModel.point_source(ps_result)
print np.sum(image_ps)
print ps_result
# for point sources, the fluxes in 'point_amp' are equivalent to the flux counts in the image.
# The only difference is the smaller cutout size in the image

# flux count in host galaxy
image_host = imageModel.source_surface_brightness(source_result)
print np.sum(image_host)

# if we only want the first component (disk in our case), we can do that
image_disk = imageModel.source_surface_brightness(source_result, k=0)
print np.sum(image_disk)

# and if we only want the second component (buldge in our case)
image_buldge = imageModel.source_surface_brightness(source_result, k=1)
print np.sum(image_buldge)

# to summarize
print("quasar-to-host galaxy ratio: ", np.sum(image_ps)/np.sum(image_host))
print("buldge-to-disk ratio:", np.sum(image_buldge)/np.sum(image_disk))

from lenstronomy.Workflow.parameters import Param
param = Param(kwargs_model, kwargs_constraints, kwargs_fixed_source=fixed_source, kwargs_fixed_ps=fixed_ps)

mcmc_new_list = []
labels_new = [r"Quasar flux", r"Buldge", r"source_x", r"source_y"]
for i in range(len(samples_mcmc)):
    # transform the parameter position of the MCMC chain in a lenstronomy convention with keyword arguments #
    kwargs_lens_out, kwargs_light_source_out, kwargs_light_lens_out, kwargs_ps_out = param.getParams(samples_mcmc[i])
    image_reconstructed, _, _, _ = imageModel.image_linear_solve(kwargs_source=kwargs_light_source_out, kwargs_ps=kwargs_ps_out)

    image_ps = imageModel.point_source(kwargs_ps_out)
    flux_quasar = np.sum(image_ps)
    image_disk = imageModel.source_surface_brightness(kwargs_light_source_out, k=0)
    flux_disk = np.sum(image_disk)
    source_x = kwargs_ps_out[0]['ra_image']
    source_y = kwargs_ps_out[0]['dec_image']
#    image_buldge = imageModel.source_surface_brightness(kwargs_light_source_out, k=1)
#    flux_buldge = np.sum(image_buldge)
    kwargs_ps_out
    mcmc_new_list.append([flux_quasar, flux_disk, source_x, source_y])

plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)