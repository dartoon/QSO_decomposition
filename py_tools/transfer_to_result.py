#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 15:59:47 2018

@author: dartoon

Transfer the fit to a result dict
"""
import numpy as np
def transfer_to_result(data, source_result, ps_result, image_ps, image_host,
                       data_C_D, cut, filt, fixcenter,ID, plot_compare=False,
                       QSO_msk = "QSO_msk*.reg",pix_sz = 'swarp', QSO_msk_image = None,
                       tag=None):
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
    if filt == 'F140w':
        zp = 26.4524
    elif filt == 'F125w':
        zp = 26.2303
    elif filt == 'acs':
        zp = 25.94333  # The AB zp for F814w, after 2006/7/4 
    result['host_mag'] = - 2.5 * np.log10(result['host_amp']) + zp 
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
    data = data
    QSO = image_ps
    host = image_host
    flux_list = [data, QSO, host]
    label = ['data', 'QSO', 'host', 'model', 'residual']
    import glob
    mask_list = glob.glob(QSO_msk)   # Read *.reg files in a list.
#    print "mask_list,muhahah", mask_list
    fig = total_compare(label_list = label, flux_list = flux_list, target_ID = ID, pix_sz=pix_sz,
                  data_mask_list = mask_list, data_cut = cut, facility = filt,
                  plot_compare = plot_compare, msk_image = QSO_msk_image)
    if tag is not None:
        fig.savefig("{0}_SB_profile.pdf".format(tag))
    fig.show()
    
    # =============================================================================
    # Calculate reduced Chisq and save to result
    # =============================================================================
    from flux_profile import cr_mask_img
    if QSO_msk_image is None:
        QSO_mask = cr_mask_img(data, mask_list, mask_reg_cut = cut)
    else:
        QSO_mask = QSO_msk_image
    chiq_map = ((data-image_ps-image_host)/np.sqrt(data_C_D))**2 * QSO_mask
    pixels=len(data_C_D)**2 - (1-QSO_mask).sum()
    reduced_Chisq = chiq_map.sum()/pixels
    result['redu_Chisq'] = reduced_Chisq
    
    result=roundme(result)
    return result
