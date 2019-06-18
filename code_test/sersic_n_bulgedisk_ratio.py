#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:54:56 2019

@author: Dartoon

The sersic n and B/D ratio
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

def n_bulge_ratio(n):
#    n = 5., log(B/D) = 0
#    n = 0.5, log(B/D) = - 2.4
    b = -1.67752801
    a = -b/np.log10(5.)
    logbd = a * np.log10(n) + b
    return 10**logbd
    
    
print n_bulge_ratio(5)
print n_bulge_ratio(0.5)
    