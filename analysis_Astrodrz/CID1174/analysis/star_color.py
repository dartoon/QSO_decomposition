#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 29 13:55:23 2018

@author: dartoon

Infer the color information of the stars.

'drz06': deltaPix = 0.0642
'acs':   deltaPix = 0.03

Change part: set up ID, psf_list_IR, psf_list_acs
"""
import numpy as np
import sys
sys.path.insert(0,'../../../py_tools')
from cut_image import cut_image, cut_center_bright, QSO_star_mag, grab_pos, QSO_star_color
from flux_profile import pix_region,flux_in_region
import copy
import astropy.io.fits as pyfits

ID = 'CID1174'
# =============================================================================
# ### Information for the IR band
# =============================================================================
filt = 'F140w'
if filt == 'F140w':
    zp_IR = 26.4524
elif filt == 'F125w':
    zp_IR = 26.2303
deltaPix_IR = 0.0642
filename_0= 'stars_and_QSO.reg'
IR_psf_list = grab_pos(filename_0,reg_ty = 'astrodrz_06')

fitsFile = pyfits.open('../astrodrz/final_drz.fits')
img_IR = fitsFile[1].data # check the back grounp
QSO_pos_IR = IR_psf_list[-1]
count=0
psf_list_IR = copy.deepcopy(IR_psf_list[:-1])
psf_list_IR = psf_list_IR[psf_list_IR[:,1].argsort()]
flux_IR = np.zeros(len(psf_list_IR))
flux_QSO_IR = 146.409   # Use the value as taking by plan b
for i in range(len(psf_list_IR)):
    PSF, center = cut_center_bright(image=img_IR, center=psf_list_IR[i], radius=60, return_center=True)
    region = pix_region(center=center, radius=12)    # take the radius = 12
    flux_IR[i] = flux_in_region(img_IR,region=region,mode='exact')
    count += 1
mag_QSO_IR = -2.5 * np.log10(flux_QSO_IR) + zp_IR
mag_IR = -2.5 * np.log10(flux_IR) + zp_IR
QSO_star_mag(img_IR, QSO_pos_IR, mag_QSO_IR, psf_list_IR, ID=ID,mag=mag_IR, reg_ty='astrodrz_06', ifsave=False)

# =============================================================================
# #### Information for the ACS band.
# =============================================================================
zp_acs = 25.94333
deltaPix_acs = 0.03
filename_1= '../../../analysis_ACS/{0}/{0}.reg'.format(ID)
acs_psf_list = grab_pos(filename_1,reg_ty = 'acs')
fitsFile = pyfits.open('../../../Cycle25data/ACS_data/{0}_acs_I_mosaic_180mas_sci.fits'.format(ID))
img_acs = fitsFile[0].data  #- (-0.003)  # check the back grounp
QSO_pos_acs = acs_psf_list[-1]
count=0
psf_list_acs = copy.deepcopy(acs_psf_list[:-1])
psf_list_acs = psf_list_acs[psf_list_acs[:,0].argsort()]
flux_acs = np.zeros(len(psf_list_acs))
flux_QSO_acs = 78.51 # Use the value as taking by plan a, this is realiable for CID1174 ACS fits.
for i in range(len(psf_list_acs)):
    PSF, center = cut_center_bright(image=img_acs, center=psf_list_acs[i], radius=60, return_center=True)
    region = pix_region(center=center, radius=16)  # take the radius = 16
    flux_acs[i] = flux_in_region(img_acs,region=region,mode='exact')
    count += 1
mag_QSO_acs = -2.5 * np.log10(flux_QSO_acs) + zp_acs
mag_acs = -2.5 * np.log10(flux_acs) + zp_acs
QSO_star_mag(img_acs, QSO_pos_acs, mag_QSO_acs, psf_list_acs, ID=ID,mag=mag_acs, reg_ty='acs', ifsave=False)

# =============================================================================
# #### Compare the color
# =============================================================================
PSF2QSO_IR = (psf_list_IR-QSO_pos_IR)*deltaPix_IR
PSF2QSO_acs = (psf_list_acs-QSO_pos_acs)*deltaPix_acs
#Take IR, compare to acs
mag_acs_corrs = np.zeros_like(mag_IR)
mag_diff=np.zeros_like(mag_IR)
for i in range(len(PSF2QSO_IR)):
    pos_diff = PSF2QSO_acs - PSF2QSO_IR[i]
    pos_diff_sq = pos_diff[:,0]**2 + pos_diff[:,1]**2
    if pos_diff_sq.min() > 1.:
        mag_diff[i]= np.nan
        mag_acs_corrs[i] = np.nan
    else:
        PSF_corrsp=np.where(pos_diff_sq==pos_diff_sq.min())[0][0]
        mag_diff[i] = mag_IR[i] - mag_acs[PSF_corrsp]
        mag_acs_corrs[i] = mag_acs[PSF_corrsp]
QSO_mags = np.array([mag_QSO_IR,mag_QSO_acs])
QSO_star_color(img_IR, QSO_pos_IR, QSO_mags, psf_list_IR, mag_IR, mag_acs_corrs,
               mag_diff, ID=ID, reg_ty='astrodrz_06')
