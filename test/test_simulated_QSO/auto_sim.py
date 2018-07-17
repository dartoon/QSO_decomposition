#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 21:47:29 2018

@author: Dartoon

auto run the test
"""

ini_n_sq = [2,4,6]
host_ratio_sq = [0.3, 0.6, 0.9]

for i in range(len(ini_n_sq)):
    for j in range(len(host_ratio_sq)):
        sim_n,ini_n, host_ratio = 1., ini_n_sq[i], host_ratio_sq[j]
        for k in range(5):
            print 'Sn_{0}_In_{1}_Rato_{2}.txt'.format(int(sim_n),int(ini_n),int(host_ratio*10))
            runfile('/Users/Dartoon/Astro/QSO_decomposition/test/test_simulated_QSO/simulated_test.py', wdir='/Users/Dartoon/Astro/QSO_decomposition/test/test_simulated_QSO')