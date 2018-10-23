#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 21:47:29 2018

@author: Dartoon

auto run the test
"""

#import os 
#dir_path = os.path.dirname(os.path.realpath(__file__))

#Done:
#'CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',
#'CID206','CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',

#file_list= ['CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136']
#file_list= ['SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202', 'XID2396']


for i in range(len(file_list)):
    print 'Run for: ',file_list[i]
    runfile('/lhome/dxh/QSO_decomposition/analysis_Astrodrz/{0}/analysis/4_analysis.py'.format(file_list[i]),
            wdir='/lhome/dxh/QSO_decomposition/analysis_Astrodrz/{0}/analysis'.format(file_list[i]))