#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 21:47:29 2018

@author: Dartoon

auto run the test
"""

import sys
sys.path.insert(0,'../py_tools')
from filter_info import filt_info

import matplotlib
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
cmap = matplotlib.cm.get_cmap('viridis')

import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

file_list = ['CDFS-1',  'CDFS-229',  'CDFS-321',  'CDFS-724',  'CID1174',  'CID1281', 'CID206',  'CID216']
file_list = ['CID237',  'CID3242',  'CID3570',  'CID452', 'CID454',  'CID50',  'CID543']
file_list = ['CID597',  'CID607',  'CID70',  'ECDFS-358', 'LID1273',  'LID1538',  'LID360',  'SXDS-X1136']
file_list = ['SXDS-X50',  'SXDS-X717', 'SXDS-X735',  'SXDS-X763',  'SXDS-X969',  'XID2138',  'XID2202',  'XID2396']

for i in range(len(file_list)):
    print 'Run for: ',file_list[i]
    runfile('{1}/{0}/deep_analysis/4_analysis.py'.format(file_list[i],dir_path),
            wdir='{1}/{0}/deep_analysis'.format(file_list[i],dir_path))
    
#
##import os 
##dir_path = os.path.dirname(os.path.realpath(__file__))
#
##Done:
##'CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',
##'CID206','CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',
#
##file_list= ['CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136']
##file_list= ['SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202', 'XID2396']
##
##file_list = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
##'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
##'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
##'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
##'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724'\
##]
##import matplotlib as matt
##matt.rcParams['font.family'] = 'STIXGeneral'
##
##for i in range(len(file_list)):
##    print 'Run for: ',file_list[i]
##    runfile('/Users/Dartoon/Astro/QSO_decomposition/analysis_Astrodrz/{0}/analysis/best_fit_plot_for_paper.py'.format(file_list[i]),
##            wdir='/Users/Dartoon/Astro/QSO_decomposition/analysis_Astrodrz/{0}/analysis'.format(file_list[i]))
#    
#
#'''
##
##file_list = ['CID1174', 'CID206', 'CID216', 'CID237', 'CID3242', 'CID3570', 'CID452',\
##'CID454', 'CID50', 'CID543', 'CID607', 'CID70', 'XID2138', 'XID2202', 'XID2396',\
##'LID1273', 'LID1538', 'LID360', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735',\
##'SXDS-X763', 'SXDS-X969', 'CDFS-1', 'CDFS-229', 'CDFS-321', 'CDFS-724', 'ECDFS-358']
##for i in range(len(file_list)):
##    print "{0}_SB_profile.pdf}}".format(file_list[i])
#    
#name_list = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 'CID3570', 'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 'CID70', 'LID1273', 'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396', 'CDFS-1', 'CDFS-229', 'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']
##for i in range(len(name_list)):
###    print name_list[i]
##    if name_list[i] in file_list:
##        print filt_info[name_list[i]]
##    else:
##        print '*'
#ACS_IDs = ['CID1174','CID216', 'CID50','CID70','XID2138','CID3242',\
#'LID1273','XID2202','CID206','CID3570','CID543','LID1538','XID2396','CID452',\
#'LID360','CID237','CID454','CID607']        
#for i in range(len(name_list)):
##    print name_list[i]
#    if name_list[i] in ACS_IDs:
#        print name_list[i], "+ACS/F814W"
#    else:
#        print '\n'
#'''