#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 17:21:53 2018

@author: Dartoon
"""
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 20:38:51 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from subprocess import call
#print call('ls', shell=True)

#IDs = ['CID1174','CID216', 'CID50','CID70','XID2138','CID3242',\
#'LID1273','XID2202','CID206','CID3570','CID543','LID1538','XID2396','CID452',\
#'LID360','CID237','CID454','CID607']
#
##CID1281, CID255, CID526, CID597 is not in our targets
##
#for key in IDs:
#    ID = key
##    print key
#    print call("cp {0}/4_analysis.py {0}/4_analysis_fix_Re_n.py".format(ID), shell=True)

#ID = ''
#runfile('/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}/1_cal_Err.py'.format(ID),
#        wdir='/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}'.format(ID))
#runfile('/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}/2_creat_msk.py'.format(ID),
#        wdir='/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}'.format(ID))
#runfile('/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}/3_esti_img_bkg.py'.format(ID),
#        wdir='/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}'.format(ID))
    

#CID1281, CID255, CID526, CID597 is not in our targets
#CID50, LID1273,LID1538,LID360

#run_list = ['CID70','XID2138','CID3242']
#run_list = ['XID2202','CID206','CID3570']
#run_list = ['XID2396','CID452','CID543']
#run_list = ['CID237','CID454','CID607']
#for i in range(len(run_list)):
#    print 'Run for: ',run_list[i]
#    runfile('/lhome/dxh/QSO_decomposition/analysis_ACS/{0}/4_analysis.py'.format(run_list[i]),
#            wdir='/lhome/dxh/QSO_decomposition/analysis_ACS/{0}'.format(run_list[i]))

run_files = ['CID1174','CID216', 'CID50','CID70','XID2138']
run_files = ['CID3242', 'LID1273','XID2202','CID206']
run_files = ['CID543','LID1538','XID2396','CID452', 'LID360']
run_files = ['CID237','CID454','CID607','CID3570']
for i in range(len(run_files)):
    print 'Run for: ',run_files[i]
    run_files('/lhome/dxh/QSO_decomposition/analysis_ACS/{0}/4_analysis_fix_Re_n.py'.format(run_files[i]),
            wdir='/lhome/dxh/QSO_decomposition/analysis_ACS/{0}'.format(run_files[i]))
#for key in run_files:
#    ID = key
##    print key
#    print "sed -n 's/fit_result_each_fix/&/p' {0}/4_analysis_fix_Re_n.py".format(ID)
