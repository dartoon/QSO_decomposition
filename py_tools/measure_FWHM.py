#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 11:41:11 2018

@author: Dartoon
"""
import numpy as np

def measure_FWHM(image):
    center= len(image)/2
    line_range= (center-10,center+10)
    seed = range(line_range[0],line_range[1])
    frm = len(image)
    q_frm = frm/4
    center = np.where(image == image[q_frm:-q_frm,q_frm:-q_frm].max())[0][0]
    x_n = np.asarray([image[i][center] for i in seed]) # The x value, vertcial 
    y_n = np.asarray([image[center][i] for i in seed]) # The y value, horizontal 
    #popt, pcov = curve_fit(func, x, yn)
    from astropy.modeling import models, fitting
    g_init = models.Gaussian1D(amplitude=1., mean=60, stddev=1.)
    fit_g = fitting.LevMarLSQFitter()
    g_x = fit_g(g_init, seed, x_n)
    g_y = fit_g(g_init, seed, y_n)
    FWHM_ver = g_x.stddev.value * 2.355  # The FWHM = 2*np.sqrt(2*np.log(2)) * stdd = 2.355*stdd
    FWHM_hor = g_y.stddev.value * 2.355
    if (FWHM_hor-FWHM_ver)/FWHM_ver > 0.20:
        print "Warning, this image have inconsistent FWHM", FWHM_ver, FWHM_hor
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    seed = np.linspace(seed.min(),seed.max(),50)
#    ax.plot(seed, g_x(seed), c='r', label='Gaussian')
#    ax.legend()
#    ax.scatter(x, sample)
#    plt.show()
    return (FWHM_ver+FWHM_hor)/2., FWHM_ver, FWHM_hor

