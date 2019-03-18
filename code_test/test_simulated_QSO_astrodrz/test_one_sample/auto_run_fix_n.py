#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 21:47:29 2018

@author: Dartoon

auto run the test
"""
import numpy as np
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

fix_n_list = np.arange(1.65,1.8,0.01)
#fix_n_list = np.arange(1,2.51,0.1)

for i in range(len(fix_n_list)):
    fix_n = fix_n_list[i]
    print 'fix_n',fix_n
    runfile('{0}/analysis_sim_fix_n.py'.format(dir_path), wdir=dir_path)
            
