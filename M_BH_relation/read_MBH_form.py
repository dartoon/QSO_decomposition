#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 09:34:27 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

f = open("fmos_MBH_table","r")
with f as g:
    lines = g.readlines()

porp_list = lines[0].replace('#','').split(' ')
samples = [lines[i].split(' ') for i in range(1,len(lines))]
#for i in range(len(samples)):

ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358',\
]

#==============================================================================
# ############ load the find the serial NO. for the  list##########
#==============================================================================
ID_ser_dic =  {}
#XID2202 to LID1622
#XID2138 to LID1820
#XID2396 to LID1878
#CDFS321 to ECDFS321
MB_ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'ECDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','LID1820','LID1622',\
'LID1878', 'CID206', 'ECDFS-358']
for j in range(len(ID)):
    count = 0
    for i in range(len(samples)):
        if samples[i][1] == MB_ID[j]:
            ID_ser_dic.update({ID[j]:i})
            count += 1
    if count == 0:
        ID_ser_dic.update({ID[j]: -99})

##==============================================================================
## Print on the props on one sample
##==============================================================================
#tar_in = 3
#t_name = ID[tar_in]
#ser = ID_ser_dic[t_name]
#print "information for {0}".format(t_name)
#for i in range(len(samples[0])):
#    print 'serial{0}'.format(i), porp_list[i], samples[ser][i]
#
##==============================================================================
## Comparing the Ha Hadr, Hb, Hbdr
##==============================================================================
#for target in ID:
#    t_name = target
#    if ID_ser_dic[t_name] != -99:
#        ser = ID_ser_dic[t_name]
#        print 'target, Ha Hadr, Hb, Hbdr', float(samples[ser][5])-float(samples[ser][6]), float(samples[ser][12])-float(samples[ser][13])
#
#for tar_in in range(len(ID)):   
#    #==============================================================================
#    # test M_BH by Ha
#    #6.71+0.48*(43.73806-42)+2.12*np.log10(4481.164/1000)
#    #==============================================================================
#    t_name = ID[tar_in]
#    ser = ID_ser_dic[t_name]
#    print "Cal MBH_Ha for {0}".format(t_name)
#    if samples[ser][10] != 0:
#        FWMH_a = float(samples[ser][8])
#        logLHadr = float(samples[ser][6])
#        cal_logMa = 6.71+0.48*(logLHadr-42)+2.12*np.log10(FWMH_a/1000)
#        print float(cal_logMa) - float(samples[ser][10])
#    
#    #==============================================================================
#    # test M_BH by Ha
#    #6.71+0.48*(43.73806-42)+2.12*np.log10(4481.164/1000)
#    #==============================================================================
#    t_name = ID[tar_in]
#    ser = ID_ser_dic[t_name]
#    if samples[ser][21] != 0:
#        FWMH_b = float(samples[ser][19])
#        logL5100dr = float(samples[ser][16])
#        cal_logMb = 6.91+0.5*(logL5100dr-44)+2.*np.log10(FWMH_b/1000)
#        if float(samples[ser][21]) != 0:
#            print "Cal MBH_Hb for {0}".format(t_name)
#            print float(cal_logMb) - float(samples[ser][21])

diff = []
for i in range(len(samples)):
    if float(samples[i][10]) != 0 and float(samples[i][21]) != 0:
        print i 
        FWMH_a = float(samples[i][8])
        logLHadr = float(samples[i][6])
        cal_logMa = 6.71+0.48*(logLHadr-42)+2.12*np.log10(FWMH_a/1000)
        
        FWMH_b = float(samples[i][19])
        logL5100dr = float(samples[i][16])
        cal_logMb = 6.91+0.5*(logL5100dr-44)+2.*np.log10(FWMH_b/1000)
        diff.append(cal_logMa-cal_logMb)
print diff
diff = np.asarray(diff)

