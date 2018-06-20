#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 09:58:24 2018

@author: dartoon

Transfer the dict to the calculated Magnitude
"""
from math import exp,pi,log10
from scipy.special import gamma
import copy
import numpy as np

def getMag(Sersic,filt,deltaPix=0.0624):
    if filt == 140:
        zp = 26.4524
    elif filt == 125:
        zp = 26.2303
    elif filt == 814:
        zp = 25.94333
    n=Sersic['n_sersic']
    re = Sersic['R_sersic']/deltaPix
    k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
    cnts = (re**2)*Sersic['I0_sersic']*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi
#    print 'conts_getmag', cnts
    mag=-2.5*log10(cnts) + zp
    return mag

result = input('data save in result.txt:\n')
filt = input('which filt? 140? 125? 814? :\n')

Sersic = copy.deepcopy(result)
Sersic['R_sersic'] *= np.sqrt(Sersic['q'])

print getMag(Sersic, filt=filt, deltaPix=0.0642)