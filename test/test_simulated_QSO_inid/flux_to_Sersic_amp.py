#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 20:36:34 2018

@author: Dartoon
"""
from scipy.special import gamma
from math import exp,pi

def getAmp(Sersic,flux,deltaPix=None):
#    mag=SERSIC_in_mag['mag_sersic']
    n=Sersic['n_sersic']
    re = Sersic['R_sersic']/deltaPix
    k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
    flux= flux
#    print 'conts_getamp', cnts
    amp= flux/((re**2)*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi)
    return amp