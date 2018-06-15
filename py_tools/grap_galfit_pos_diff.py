#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 14:44:18 2018

@author: Dartoon

Read the difference between QSO and host from Galfit.01
"""

import numpy as np

filt = 'F140w'
pix_sz = 'drz06'

if pix_sz == 'drz06':
    deltaPix = 0.0642
elif pix_sz == 'acs':
    deltaPix = 0.03

#ID = 'CID70' #raw_input('Target ID')
#==============================================================================
# Grab file
#==============================================================================
folder = '../analysis_Astrodrz/XID2396/for_galfit/'
filename = 'galfit.02' #'result_QSO_02' 

fit_out = open(folder+filename,'r')
lines = fit_out.readlines()
host_x, host_y = float(lines[40][4:12]),  float(lines[40][13:20])
QSO_x, QSO_y = float(lines[53][4:12]),  float(lines[53][13:20])
dx_arc = (host_x-QSO_x)*deltaPix
dy_arc =  (host_y-QSO_y)*deltaPix
print round(dx_arc,3),',',round(dy_arc,3)