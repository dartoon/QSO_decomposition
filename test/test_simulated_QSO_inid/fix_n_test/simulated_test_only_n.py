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

## import main simulation class of lenstronomy
#from lenstronomy.SimulationAPI.simulations import Simulation
#SimAPI = Simulation()
#
#numPix = 61  #  cutout pixel size
#deltaPix = 0.13  #  pixel size in arcsec (area per pixel = deltaPix**2)
#psf_type = 'PIXEL'  # 'gaussian', 'pixel', 'NONE'
#PSF0 = pyfits.getdata('CID50_PSF.fits')
#PSF0 /= PSF0.sum()
#
#kernel_size = len(PSF0)
#kernel = PSF0
#
## initial input simulation
## generate the coordinate grid and image properties
#data_class = SimAPI.data_configure(numPix, deltaPix, 0, 0)
## generate the psf variables
#psf_class = SimAPI.psf_configure(psf_type=psf_type, fwhm=0.1, kernelsize=kernel_size, deltaPix=deltaPix, kernel=kernel)
#
## quasar center (we chose a off-centered position)
#center_x = 0.03
#center_y = 0.02
## quasar brightness (as measured as the sum of pixel values)
#total_flux = 400.
#
#sim_n,ini_n, host_ratio = 1., 2., 0.3
#
#point_amp = total_flux * (1-host_ratio)
#from lenstronomy.PointSource.point_source import PointSource
#point_source_list = ['UNLENSED']
#pointSource = PointSource(point_source_type_list=point_source_list)
#kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]
#
#from lenstronomy.LightModel.light_model import LightModel
#light_model_list = ['SERSIC_ELLIPSE']
#lightModel = LightModel(light_model_list=light_model_list)
#import lenstronomy.Util.param_util as param_util
#e1, e2 = param_util.phi_q2_ellipticity(phi=0.3, q=0.8)
#from flux_to_Sersic_amp import getAmp
#
#kwargs_disk = {'n_sersic': sim_n, 'R_sersic': 0.3, 'e1': e1, 'e2': e2, 'center_x': center_x, 'center_y': center_y}
#flux = total_flux *  host_ratio
#amp = getAmp(kwargs_disk,flux * 6.23/5.,deltaPix=deltaPix)
#kwargs_disk['amp'] = amp  #
#kwargs_host = [kwargs_disk]
##==============================================================================
##Make simulation 
##==============================================================================
#from lenstronomy.ImSim.image_model import ImageModel
#
#kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
#imageModel = ImageModel(data_class, psf_class, source_model_class=lightModel,
#                                point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
## simulate image with the parameters we have defined above #
#image = imageModel.image(kwargs_source=kwargs_host, kwargs_ps=kwargs_ps)
##plt.matshow(np.log10(image), origin='lower')
##plt.show()
## we can also add noise #
#import lenstronomy.Util.image_util as image_util
#poisson = image_util.add_poisson(image, exp_time=400.)
background_rms = 0.09
#bkg = image_util.add_background(image, sigma_bkd=background_rms)
#noise = bkg + poisson
#image_noisy = image + noise
##plt.matshow(np.log10(image_noisy), origin='lower')
##plt.show()
#
## we can also simulate the different components separately
#imageModel_ps = ImageModel(data_class, psf_class, point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
#image_ps_sim = imageModel_ps.image(kwargs_ps=kwargs_ps)
##plt.matshow(np.log10(image_ps), origin='lower')
##plt.show()
#imageModel_host = ImageModel(data_class, psf_class, source_model_class=lightModel, kwargs_numerics=kwargs_numerics)
#image_host_sim = imageModel_host.image(kwargs_source=kwargs_host)
##plt.matshow(np.log10(image_host), origin='lower')
##plt.show()
##print 'total_flux:', image_noisy.sum()
##print 'host_total_flux: ',image_host_sim.sum(),'\npoint source flux: ',image_ps_sim.sum()
##host_ratio = image_host_sim.sum()/(image_host_sim.sum()+image_ps_sim.sum())
##print 'host/total ratio:', host_ratio
#==============================================================================
# To do the fitting:
#==============================================================================

image_noisy=pyfits.getdata('sim_data.fits')
QSO_std = pyfits.getdata('noise_map.fits')
background_rms = 0.09
#fix_n =1
fixed_source, kwargs_source_init,kwargs_source_sigma,kwargs_lower_source, kwargs_upper_source= [], [], [], [], []
fixed_source.append({'n_sersic': fix_n,'R_sersic': 0.3})#, 'center_x': -0.027,'center_y': -0.003, 'qso_x': 0.009,'qso_y': -0.017,'e1': -0.062592224365596166, 'e2': 0.080528893354289297})
kwargs_source_init.append({'R_sersic': .3, 'n_sersic': 1, 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
kwargs_source_sigma.append({'n_sersic_sigma': 1, 'R_sersic_sigma': 0.5, 'e1_sigma': 0.1, 'e2_sigma': 0.1, 'center_x_sigma': 0.1, 'center_y_sigma': 0.1})
kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -10, 'center_y': -10})
kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': 10, 'center_y': 10})
source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]

import sys
sys.path.insert(0,'../../../py_tools/')
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result
#QSO_std = pyfits.getdata('{0}_Err.fits'.format(ID))
PSF1 = pyfits.getdata('PSF1.fits')
PSF1 /= PSF1.sum()

#from flux_profile import QSO_psfs_compare, profiles_compare
#fig_pro_compare = profiles_compare([PSF0,PSF1], [1,1], prf_name_list=['PSF0','PSF1'],
#                                   gridspace = 'log',if_annuli=True,norm_pix=2.0)
#plt.show()

fixcenter = False
#tag = 'test'
tag = 'fix_n{0}_only_n'.format(int(fix_n))
ID = 'sim_test'
source_result, ps_result, image_ps, image_host, data_C_D=fit_qso(QSO_im=image_noisy, psf_ave=PSF1, source_params=source_params,
                                                                 background_rms=background_rms, image_plot = True, corner_plot=False,
                                                                 flux_ratio_plot=False, deep_seed = False, fixcenter= fixcenter,
                                                                 exp_time = 400.,tag=tag, no_MCMC=True, QSO_std=QSO_std)

result = transfer_to_result(data=image_noisy,
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, data_C_D=QSO_std**2,
                            cut=0, filt='F140w', fixcenter=fixcenter,ID=ID,
                            QSO_msk = 'QSO_msk*.reg', plot_compare= True, tag=tag) 
import glob
filename = 'fix_n_test_only_n.txt'
if_file = glob.glob(filename)   
if if_file == []:
    f =  open(filename,'w') 
elif if_file is not []:
    f = open(filename,"r+")
    f.read()
f.write("\n#========================================")    
f.write("\n#fix_n: "+repr(fix_n)+ ' red_Chisq: '+repr(result['redu_Chisq'])+ ' host_ratio: '+repr(round(result['host_flux_ratio_percent'],1))+'%')
f.close()
f=open(filename,"r+")
f.read()
