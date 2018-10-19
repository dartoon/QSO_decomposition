#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 21:47:29 2018

@author: Dartoon

auto run the test
"""

#import os 
#dir_path = os.path.dirname(os.path.realpath(__file__))

file_list=['CID237', 'CID70', 'LID360', 'SXDS-X763', 'CDFS-1', 'CID3242', 'SXDS-X969',\
           'CDFS-229', 'CID3570', 'XID2138', 'CDFS-321', 'CID452', 'XID2202', 'CID454',\
           'SXDS-X1136', 'XID2396', 'CID1174', 'CID50', 'SXDS-X50', 'CID206', 'CID543',\
           'LID1273', 'SXDS-X717', 'CID216', 'CID607', 'LID1538', 'SXDS-X735']

for i in range(len(file_list)):
    print 'Run for: ',file_list[i]
    runfile('/Users/Dartoon/Astro/QSO_decomposition/analysis_Astrodrz/{0}/analysis/analysis.py'.format(file_list[i]),
            wdir='/Users/Dartoon/Astro/QSO_decomposition/analysis_Astrodrz/{0}/analysis'.format(file_list[i]))