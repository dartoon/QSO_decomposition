#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 19:59:43 2018

@author: dartoon

Transfer the dict to a list for out put
"""

import numpy as np

result = input('data save in result.txt:\n')
if 'qso_x' not in result:
    print  result['host_mag'], '\t', result['n_sersic'],'\t', result['R_sersic'],'\t', result['host_flux_ratio_percent'],'\t', \
result['host_flux_ratio_percent'], '\t', result['redu_Chisq']
elif 'qso_x' in result:
    dx = result['center_x']-result['qso_x']
    dy = result['center_y']-result['qso_y']
    print result['host_mag'], '\t', result['n_sersic'],'\t', result['R_sersic'],'\t', result['host_flux_ratio_percent'],'\t', \
result['host_flux_ratio_percent'], '\t', result['redu_Chisq'],'\t', dx, ',', dy