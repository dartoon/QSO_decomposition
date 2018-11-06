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
sys.path.insert(0,'../../py_tools')
from psfs_average import psf_ave
from flux_profile import QSO_psfs_compare, profiles_compare
from matplotlib.colors import LogNorm

import os
path = os.getcwd()
ID = path.split('/')[-1]

IDs = ['CID1174','CID255','CID50','CID70','XID2138','CID1281','CID3242','CID526',\
'LID1273','XID2202','CID206','CID3570','CID543','LID1538','XID2396','CID216','CID452',\
'CID597','LID360','CID237','CID454','CID607']

filt = 'acs'
# =============================================================================
# Read PSF and QSO image
# =============================================================================
QSO_bkg_value= 0.
QSO_im = pyfits.getdata('{0}_cutout.fits'.format(ID)) - QSO_bkg_value
QSO_msk = pyfits.getdata('{0}_msk.fits'.format(ID))
frame_size = 61
#frame = '{0}'.format(frame_size)
QSO_fm = len(QSO_im)
ct = (QSO_fm-frame_size)/2     # If want to cut to 61, i.e. (121-61)/2=30
        
import pickle
PSFs_dict = {}
QSOs_dict = {}
for key in IDs:
    PSFs, _=pickle.load(open('../{0}/{0}_PSFs_QSO'.format(key),'rb'))
    PSFs_dict.update({'{0}'.format(key):PSFs})
#    QSOs_dict.update({'{0}'.format(key):QSOs})

PSF_list = []
PSF_id = []
for key in IDs:
    psfs_dict = PSFs_dict[key]
    psfs = [psfs_dict[i] for i in range(len(psfs_dict))]
    PSF_list += psfs
    name_id = [key+"_"+str(i) for i in range(len(psfs_dict))]
    PSF_id = PSF_id + name_id
    if len(PSF_list) != len(PSF_id):
        raise ValueError("The PSF_list is not consistent with PSF_id")
psf_list = [PSF_list[i][0] for i in range(len(PSF_list))]
PSF_mask_img_list = [PSF_list[i][3] for i in range(len(PSF_list))]
psf_name_list = PSF_id
# =============================================================================
# Doing the fitting
# =============================================================================
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result
#from flux_profile import cr_mask_img
QSO_outer = pyfits.getdata('{0}_cutout_outer.fits'.format(ID))
from photutils import make_source_mask
mask = make_source_mask(QSO_outer, snr=2, npixels=5, dilate_size=11)
plt.imshow(QSO_outer* (1-mask*1), origin='low')
plt.close()
background_rms = np.std(QSO_outer* (1-mask*1))
print "background_rms: ", background_rms
QSO_msk = QSO_msk[ct:-ct,ct:-ct]
QSO_im = QSO_im[ct:-ct,ct:-ct]
QSO_msk = QSO_msk*0 +1    # This means no mask is added
QSO_std = pyfits.getdata('wht_err.fits')[ct:-ct,ct:-ct]
##############################Fit
#==============================================================================
# Activate the following part if their is no objects around the QSO.
#==============================================================================
from mask_objects import detect_obj
objs, Q_index = detect_obj(QSO_im)
obj = [objs[i] for i in range(len(objs)) if i != Q_index]
pix_s = 0.03

fixed_source = []
kwargs_source_init = []
kwargs_source_sigma = []
kwargs_lower_source = []
kwargs_upper_source = []      

fix_re, fix_n = 0.49, 1.293
fixed_source.append({'R_sersic': fix_re,'n_sersic': fix_n})  
kwargs_source_init.append({'R_sersic': fix_re, 'n_sersic': fix_n, 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
kwargs_source_sigma.append({'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'center_x': -0.5, 'center_y': -0.5})
kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'center_x': 0.5, 'center_y': 0.5})

## If want to add obj manually
#pos = np.array([15, 46])
#pos -= np.array([len(QSO_im)/2, len(QSO_im)/2])
#pos = pos * pix_s
#fixed_source.append({})  
#kwargs_source_init.append({'R_sersic': 0.1, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -pos[0], 'center_y': pos[1]})
#kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
#kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.01, 'n_sersic': 0.3, 'center_x': -pos[0]-0.5, 'center_y': pos[1]-0.5})
#kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': -pos[0]+0.5, 'center_y': pos[1]+0.5})

#if len(obj) >= 1:
#    for i in range(len(obj)):
#        fixed_source.append({})  
#        kwargs_source_init.append({'R_sersic': obj[i][2] * pix_s, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_s, 'center_y': obj[i][0][1]*pix_s})
#        kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
#        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': obj[i][2] * pix_s/5, 'n_sersic': 0.3, 'center_x': -obj[i][0][0]*pix_s-10, 'center_y': obj[i][0][1]*pix_s-10})
#        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': -obj[i][0][0]*pix_s+10, 'center_y': obj[i][0][1]*pix_s+10})

source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]

if os.path.exists('fit_result_each_fix')==False:
    os.mkdir('fit_result_each_fix')

import time
t1 = time.time()
fixcenter = False
filename = 'fit_result_each_fix/each_PSF_fit_qso.txt'
if_file = glob.glob(filename)   
if if_file == []:
    fit_result =  open(filename,'w') 
elif if_file is not []:
    fit_result = open(filename,"r+")
    fit_result.read()
count = 0
for i in np.array(range(len(psf_name_list))):
    print "{2} by PSF: {0}, PSF NO. {1}".format(psf_name_list[i], i, ID)
    tag = 'fit_result_each_fix/qso_fit_PSF{0}'.format(i)
    psf_i = psf_list[i] * PSF_mask_img_list[i]
    psf_i = psf_i[ct:-ct,ct:-ct]
    source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_im, psf_ave=psf_i, psf_std = None,
                                                                     background_rms=background_rms,
                                                                     source_params=source_params, QSO_msk = QSO_msk, fixcenter=fixcenter,
                                                                     pix_sz = 'acs', no_MCMC =True,
                                                                     QSO_std =QSO_std, tag=tag)
    result = transfer_to_result(data=QSO_im, pix_sz = 'acs',
            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, error_map=error_map,
            filt=filt, fixcenter=fixcenter,ID=ID,QSO_msk =QSO_msk, tag=tag)
    if count == 0:
        fit_result.write("#QSO_img intensity: {0} \n".format(round(np.sum(QSO_im*QSO_msk),2)))
    fit_result.write("#fit by PSF{0}: \n".format(psf_name_list[i]))
    fit_result.write('PSF_intensity:{0} \n'.format(round(np.sum(psf_i),2)))
    fit_result.write(repr(result) + "\n")
    count += 1
    t_i = time.time()
    print "assuming time need for this target:", (t_i-t1)/(i+1)*len(psf_name_list)/60, 'mins'
    print "time remaining:", "{0}%".format(round(float(i+1)/len(psf_name_list)*100,2)),'mins left:',round((t_i-t1)/(i+1)*len(psf_name_list)/60*(float(len(psf_name_list)-i-1)/len(psf_name_list)),2)
fit_result.close()
t2 = time.time()
t_tot= (t2-t1)/60
print "total time:", t_tot, "mins"

import os
os.system('say "your program has finished"')
