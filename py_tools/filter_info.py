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
'SXDS-X763': 'F125w', 'SXDS-X969': 'F140w', 'CDFS-321': 'F140w',
'CID1281': 'F140w', 'CID597':'F125w', 'CID255':'F140w'\
}

redshift_info = {'CID50': 1.239, 'CID206': 1.483, 'CID216': 1.567, 'CID237': 1.618,\
'CID255': 1.664, 'CID452': 1.407, 'CID597': 1.272, 'CID607': 1.294, 'CID1174': 1.552,\
'CID1281': 1.445, 'CID3242': 1.532, 'CID3570': 1.244, 'XID2396': 1.600, 'CID70': 1.667,\
'LID1273': 1.617, 'XID2202': 1.516, 'CID454': 1.478, 'CID543': 1.301, 'XID2138': 1.551,
'LID1538': 1.527, 'CDFS-1': 1.630, 'CDFS-724': 1.337, 'LID360': 1.579, 'CDFS-321': 1.570,\
'CDFS-229': 1.326, 'ECDFS-358': 1.626, 'SXDS-X50': 1.411, 'SXDS-X717': 1.276,\
'SXDS-X735': 1.447, 'SXDS-X763': 1.412, 'SXDS-X969': 1.585, 'SXDS-X1136': 1.325}