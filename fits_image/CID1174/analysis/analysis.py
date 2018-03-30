#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 20:17:41 2018

@author: Dartoon

On the analysis of CID1174
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pylab as plt
import glob
import sys
sys.path.insert(0,'../../../py_tools')
from psfs_average import psf_ave
from flux_profile import QSO_psfs_compare, profiles_compare
from matplotlib.colors import LogNorm

ID = 'CID1174'

# =============================================================================
# Read PSF and QSO image
# =============================================================================
psf_name_list = glob.glob("PSF*.fits")   # Read *.reg files in a list.
psf_list = []
for i in range(len(psf_name_list)):
    psf_get = pyfits.getdata('PSF{0}.fits'.format(i))
    psf_list.append(psf_get)
mask_list = glob.glob("PSF*.reg")   # Read *.reg files in a list.
QSO_im = pyfits.getdata('{0}_cutout.fits'.format(ID))

#==============================================================================
# Compare the profile and derive the Average image
#==============================================================================
cut = 20      #cut_range
fig = QSO_psfs_compare(QSO=QSO_im[cut:-cut,cut:-cut], psfs=psf_list,
#                 plt_which_PSF=(0,1,2,3,4,5,6),
                 mask_list=mask_list,
                 include_QSO=False, radius=len(psf_list[0])/2, grids=20,
                 gridspace= 'log')

psf_ave_dirt, psf_std_dirt=psf_ave(psf_list,mode = 'direct', not_count=(5,6),
                  mask_list=mask_list)

psf_ave_wght, psf_std_wght=psf_ave(psf_list,mode = 'CI', not_count=(5,6),
                  mask_list=mask_list)

prf_list = [QSO_im,psf_ave_dirt, psf_ave_wght]
scal_list = [1,1,1]
prf_name_list = ['QSO', 'PSF_ave_direct', 'PSF_ave_by_wght']
profiles_compare(prf_list, scal_list, prf_name_list=prf_name_list, gridspace = 'log')

from fit_qso import fit_qso
#print "by psf_ave_dirt"
#source_result, ps_result, image_ps, image_host=fit_qso(QSO_im[20:-20,20:-20], psf_ave=psf_ave_dirt,
#                                                       source_params=None, image_plot = True, corner_plot=True, flux_ratio_plot=True)

print "by psf_ave_wght"
source_result, ps_result, image_ps, image_host=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_ave_wght, psf_std = psf_std_wght,
                                                       source_params=None, image_plot = True, corner_plot=True, flux_ratio_plot=True,
                                                       deep_seed = False)
plt.show()
#==============================================================================
# Translate the e1, e2 to phi_G and q
#==============================================================================
import lenstronomy.Util.param_util as param_util
source_result[0]['phi_G'], source_result[0]['q'] = param_util.ellipticity2phi_q(source_result[0]['e1'], source_result[0]['e2'])

#==============================================================================
# Save the result
#==============================================================================
from roundme import roundme
import copy
result = copy.deepcopy(source_result[0])
del result['e1']
del result['e2']
result['QSO_amp'] = ps_result[0]['point_amp'][0]
result['host_amp'] = image_host.sum()
result['host_flux_ratio_percent']= image_host.sum()/(image_ps.sum() + image_host.sum())*100
zp = 26.4524
result['host_mag'] = - 2.5 * np.log10(result['host_amp']) + zp 
result=roundme(result)
#print "The host flux is ~:", image_host.sum()/(image_ps.sum() + image_host.sum())

##==============================================================================
##Plot the images for adopting in the paper
##==============================================================================
from flux_profile import total_compare
data = QSO_im[cut:-cut,cut:-cut]
QSO = image_ps
host = image_host
flux_list = [data, QSO, host]
label = ['data', 'QSO', 'host', 'model', 'residual']
import glob
mask_list = glob.glob("QSO*.reg")   # Read *.reg files in a list.
total_compare(label_list = label, flux_list = flux_list, target_ID = ID,
              data_mask_list = mask_list, data_cut = cut, facility = 'F140w')
fig.savefig("SB_profile_{0}.pdf".format(ID))
