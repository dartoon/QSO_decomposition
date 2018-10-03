#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 11:41:01 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

filt_info = {'CDFS-1': 'F140w', 'CDFS-229': 'F125w', 'CDFS-724': 'F125w',\
'CID1174': 'F140w', 'CID206': 'F140w', 'CID216': 'F140w', 'CID3242': 'F140w',\
'CID3570': 'F125w', 'CID452': 'F125w', 'CID454': 'F140w', 'CID50': 'F125w',\
'CID607': 'F125w', 'CID70': 'F140w', 'LID1273': 'F140w', 'LID360': 'F140w',\
'XID2138': 'F140w', 'XID2202': 'F140w', 'XID2396': 'F140w', 'ECDFS-358': 'F140w',\
'SXDS-X1136': 'F125w', 'SXDS-X50': 'F125w', 'SXDS-X735': 'F140w',\
'CID543': 'F125w', 'LID1538': 'F140w', 'CID237': 'F140w', 'SXDS-X717': 'F125w',\
'SXDS-X763': 'F125w', 'SXDS-X969': 'F140w', 'CDFS-321': 'F140w'}