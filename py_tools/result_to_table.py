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
#    print  result['host_mag'], '\t', result['n_sersic'],'\t', result['R_sersic'],'\t', result['host_flux_ratio_percent'],'%','\t', \
#, result['redu_Chisq']
    print "{0}\t{1}\t{2}\t{3}%\t{4}".format(result['host_mag'],result['n_sersic'],result['R_sersic'],result['host_flux_ratio_percent'],
           result['redu_Chisq'])
elif 'qso_x' in result:
    dx = result['center_x']-result['qso_x']
    dy = result['center_y']-result['qso_y']
#    print result['host_mag'], '\t', result['n_sersic'],'\t', result['R_sersic'],'\t', result['host_flux_ratio_percent'],'%','\t', \
#result['redu_Chisq'],'\t', dx, ',', dy
    print "{0}\t{1}\t{2}\t{3}%\t{4}\t{5},{6}".format(result['host_mag'],result['n_sersic'],result['R_sersic'],result['host_flux_ratio_percent'],
           result['redu_Chisq'],dx,dy)