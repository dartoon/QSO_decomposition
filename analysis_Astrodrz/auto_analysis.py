#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 21:47:29 2018

@author: Dartoon

auto run the test
"""

#import os 
#dir_path = os.path.dirname(os.path.realpath(__file__))

file_list=['454', '50', '70']

for i in range(len(file_list)):
    print 'Run for: CID',file_list[i]
    runfile('/Users/Dartoon/Astro/QSO_decomposition/analysis_Astrodrz/CID{0}/analysis/analysis.py'.format(file_list[i]),
            wdir='/Users/Dartoon/Astro/QSO_decomposition/analysis_Astrodrz/CID{0}/analysis'.format(file_list[i]))