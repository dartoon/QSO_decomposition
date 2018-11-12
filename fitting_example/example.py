#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 20:17:41 2018

@author: Dartoon

On the analysis of xxx
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pylab as plt
import glob
import sys
sys.path.insert(0,'../py_tools')

import os
path = os.getcwd()
ID = path.split('/')[-2]

filt = 'F140w'

# =============================================================================
# Read PSF and QSO image
# =============================================================================
QSO_im = pyfits.getdata('data.fits'.format(ID))
frame_size = 61   # The final frame size to fit the QSO image
#frame = '{0}'.format(frame_size)
QSO_fm = len(QSO_im)
ct = (QSO_fm-frame_size)/2     # If want to cut to 61, i.e. (121-61)/2=30

# =============================================================================
# Doing the fitting
# =============================================================================
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result
QSO_im = QSO_im[ct:-ct,ct:-ct]
QSO_msk = QSO_im*0 +1    # This means is no mask is added in the image
QSO_std = pyfits.getdata('error_map.fits')[ct:-ct,ct:-ct]
##############################Fit
fixcenter = False
filename = 'fit_result.txt'
if_file = glob.glob(filename)   
if if_file == []:
    fit_result =  open(filename,'w') 
elif if_file is not []:
    fit_result = open(filename,"r+")
    fit_result.read()

tag = 'result'
psf = pyfits.getdata('psf.fits')[ct:-ct,ct:-ct]
source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_im, psf_ave=psf, psf_std = None,
                                                                 background_rms=.0047221,
                                                                 source_params=None, QSO_msk = QSO_msk,
                                                                 fixcenter=fixcenter,
                                                                 pix_sz = 'drz06', no_MCMC =True,
                                                                 QSO_std =QSO_std, tag=tag)
result = transfer_to_result(data=QSO_im, pix_sz = 'drz06',
        source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, error_map=error_map,
        filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk =QSO_msk, tag=tag)

fit_result.write("#QSO_img intensity: {0} \n".format(round(np.sum(QSO_im*QSO_msk),2)))
fit_result.write('PSF_intensity:{0} \n'.format(round(np.sum(psf),2)))
fit_result.write(repr(result) + "\n")
fit_result.close()

import os
os.system('say "your program has finished"')
