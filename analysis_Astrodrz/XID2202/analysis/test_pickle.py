#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 21:20:49 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

ID = 'XID2202'
import pickle
datafile = open('{0}_PSFs_QSO'.format(ID),'rb')
PSFs, QSO=pickle.load(open('XID2202_PSFs_QSO','rb'))
datafile.close()

