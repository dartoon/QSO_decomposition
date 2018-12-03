#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 14:15:19 2018

@author: Dartoon

Test the MBH calibration by Andreas
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

f = open("fmos_MBH_table","r")
with f as g:
    lines = g.readlines()
porp_list = lines[0].replace('#','').split(' ')
samples = [lines[i].split(' ') for i in range(1,len(lines))]

two_line_sample = []
for i in range(len(samples)):
    if float(samples[i][10])!=0 and float(samples[i][21])!=0:
        two_line_sample.append(samples[i])

MB_info_a, MB_info_b = [], []
for i in range(len(two_line_sample)):       
    t_name = two_line_sample[i][1]
    FWMH_a = float(two_line_sample[i][8])
    logLHadr = float(two_line_sample[i][6])
    cal_logMa_And = 6.71+0.48*(logLHadr-42)+2.12*np.log10(FWMH_a/1000)  # as used in Andreas
    cal_logMa_my = 6.459+0.55*(logLHadr-42)+2.*np.log10(FWMH_a/1000)  # as used in H0liCOW 7 and McGill
    MB_info_a.append([t_name, FWMH_a, logLHadr, cal_logMa_And, float(two_line_sample[i][10]), cal_logMa_my])
    
    FWMH_b = float(two_line_sample[i][19])
    logL5100dr = float(two_line_sample[i][16])
    cal_logMb_And = 6.91+0.5*(logL5100dr-44)+2.*np.log10(FWMH_b/1000)  # as used in Andreas
    cal_logMb_my = 6.882+0.518*(logL5100dr-44)+2.*np.log10(FWMH_b/1000)        # calibrated in H0liCOW 7
    MB_info_b.append([t_name, FWMH_b, logL5100dr, cal_logMb_And, float(two_line_sample[i][21]), cal_logMb_my])


 # compare using Andreas's recipe to his result
#listm = []
#for i in range(len(two_line_sample)):
#    listm.append(MB_info_b[i][3]- MB_info_b[i][4])
#listm = np.asarray(listm)
#print listm.max()

plt.figure(figsize=(10, 10))
diff_And, diff_my = [], []
for i in range(len(MB_info_a)):
    plt.plot(MB_info_a[i][3], MB_info_b[i][3],'bo')
    diff_And.append([MB_info_a[i][3]-MB_info_b[i][3]])
    plt.plot(MB_info_a[i][5], MB_info_b[i][5],'ks')
    diff_my.append([MB_info_a[i][5]-MB_info_b[i][5]])
diff_And = np.asarray(diff_And)
diff_my = np.asarray(diff_my)
#texts = []
x=np.linspace(6,12,20)
y = x
plt.plot(x,y,'y')
#plt.title('', fontsize=30)
plt.xlabel('Calibrated by Ha',fontsize=25)
plt.ylabel('Calibrated by Hb', fontsize=25)
#val_min, val_max = np.min([mag_k_corrected_UV, mag_k_corrected_IR]), np.max([mag_k_corrected_UV, mag_k_corrected_IR])
plt.xlim(7.2, 9.8)
plt.ylim(7.2, 9.8)
plt.tick_params(labelsize=25)     
plt.show()

