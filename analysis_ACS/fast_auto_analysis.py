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

#IDs = ['CID255','CID50','CID70','XID2138','CID1281','CID3242','CID526',\
#'LID1273','XID2202','CID206','CID3570','CID543','LID1538','XID2396','CID452',\
#'CID597','LID360','CID237','CID454','CID607']

#CID1281, CID255, CID526, CID597 is not in our targets

#for key in IDs:
#    ID = key
##    print key
#    print call("cp ../template/temp_analysis/ACS_temp/4_*.py {0}/".format(ID), shell=True)

#ID = ''
#runfile('/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}/1_cal_Err.py'.format(ID),
#        wdir='/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}'.format(ID))
#runfile('/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}/2_creat_msk.py'.format(ID),
#        wdir='/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}'.format(ID))
#runfile('/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}/3_esti_img_bkg.py'.format(ID),
#        wdir='/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}'.format(ID))
    

#CID1281, CID255, CID526, CID597 is not in our targets
IDs = ['CID50','CID70','XID2138','CID3242',\
'LID1273','XID2202','CID206','CID3570',\
'LID1538','XID2396','CID452','CID543',\
'LID360','CID237','CID454','CID607']
for i in range(len(IDs)):
    print 'Run for: ',IDs[i]
    runfile('/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}/4_analysis.py'.format(IDs[i]),
            wdir='/Users/Dartoon/Astro/QSO_decomposition/analysis_ACS/{0}'.format(IDs[i]))