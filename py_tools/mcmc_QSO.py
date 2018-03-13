#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 21:57:56 2018

@author: Dartoon

A class for modeling the QSO image
"""

class QSO_decomp():
    def __init__(self,img=None, psf=None):
        self.img = img
        self.psf = psf
    
        