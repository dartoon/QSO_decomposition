#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 20:17:41 2018

@author: Dartoon

On the analysis of CID216
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

ID = 'CID216'
filt = 'F140w'

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
#                 plt_which_PSF=(0,1,2,3,4,5,6,7,8),
                 mask_list=mask_list,
                 include_QSO=True, grids=30)#, gridspace= 'log')
psf_ave_pa, psf_std_pa=psf_ave(psf_list,mode = 'CI', not_count=(1,2,3,8),
                  mask_list=mask_list)

psf_ave_pb, psf_std_pb=psf_ave(psf_list,mode = 'CI', not_count=(1,2,3,7,8),
                  mask_list=mask_list)


#pyfits.PrimaryHDU(psf_ave_pa).writeto('../../PSF_legacy/{0}_PSF.fits'.format(ID),overwrite=True)
#pyfits.PrimaryHDU(psf_std_pa).writeto('../../PSF_legacy/{0}_PSF_std.fits'.format(ID),overwrite=True)

PSF_1174 = pyfits.getdata("../../PSF_legacy/CID1174_PSF.fits")
PSF_1174_std = pyfits.getdata("../../PSF_legacy/CID1174_PSF_std.fits")
PSF_452= pyfits.getdata("../../PSF_legacy/CID452_PSF.fits")
PSF_452_std = pyfits.getdata("../../PSF_legacy/CID452_PSF_std.fits")

prf_list = [QSO_im,psf_ave_pa, PSF_1174, PSF_452]
scal_list = [1,1,1,1]
prf_name_list = ['QSO_CID216', 'Plan a', 'PSF_CID1174', 'PSF_CID452']
profiles_compare(prf_list, scal_list, prf_name_list=prf_name_list, gridspace = 'log')

from fit_qso import fit_qso
fixcenter = False

print "Plan a"
source_result, ps_result, image_ps, image_host,data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_ave_pa, background_rms=0.041, psf_std = psf_std_pa,
                                                       source_params=None, image_plot = True, corner_plot=True, flux_ratio_plot=True,
                                                       fixcenter= fixcenter)

#print "Plan b"
#source_result, ps_result, image_ps, image_host,data_C_D=fit_qso(QSO_im[cut:-cut,cut:-cut], psf_ave=psf_ave_pb, background_rms=0.041, psf_std = psf_std_pb,
#                                                       source_params=None, image_plot = True, corner_plot=True, flux_ratio_plot=True,
#                                                       fixcenter= fixcenter)
plt.show()
#=============================================================================
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
if filt == 'F140w':
    zp = 26.4524
elif filt == 'F125w':
    zp = 26.2303
result['host_mag'] = - 2.5 * np.log10(result['host_amp']) + zp 
result=roundme(result)
#print "The host flux is ~:", image_host.sum()/(image_ps.sum() + image_host.sum())

# =============================================================================
# Save QSO position to result if not fix center
# =============================================================================
if fixcenter == False:
    result['qso_x'] = ps_result[0]['ra_image'][0]
    result['qso_y'] = ps_result[0]['dec_image'][0]

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
fig = total_compare(label_list = label, flux_list = flux_list, target_ID = ID,
              data_mask_list = mask_list, data_cut = cut, facility = 'F140w')
#fig.savefig("SB_profile_{0}.pdf".format(ID))

# =============================================================================
# Calculate reduced Chisq and save to result
# =============================================================================
from flux_profile import cr_mask_img
QSO_mask = cr_mask_img(QSO_im[cut:-cut,cut:-cut], mask_list, mask_reg_cut = cut)
chiq_map = ((QSO_im[cut:-cut,cut:-cut]-image_ps-image_host)/np.sqrt(data_C_D))**2 * QSO_mask
pixels=len(data_C_D)**2 - (1-QSO_mask).sum()
reduced_Chisq = chiq_map.sum()/pixels
result['redu_Chisq'] = reduced_Chisq

result=roundme(result)
print result
