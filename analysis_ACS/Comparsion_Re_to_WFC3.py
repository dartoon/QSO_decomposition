#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 17:10:00 2019

@author: Dartoon

"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import re
import matplotlib
from adjustText import adjust_text   # avoid the overlapping while ploting
import sys
sys.path.insert(0,'../py_tools')
from filter_info import redshift_info

ID = ['CID1174','CID216', 'CID50','CID70','XID2138','CID3242',\
'LID1273','XID2202','CID206','CID543','LID1538','XID2396','CID452',\
'LID360','CID237','CID454','CID607','CID3570', 'CID1281', 'CID597', 'CID255']

#CID1281, CID255, CID526, CID597 is not in our targets

import pickle

Re_results= []
flux_dict, FWHM_dict, locs_dict, filter_dict, id_stars_dict=pickle.load(open('PSFs_lib_dict','rb'))
folder = "fit_result_each" #"fit_result_each_fix" or "fit_result_each"

from load_result import load_re
WFC3_re = np.array(load_re(ID, folder = '../', flt='WFC3'))


for j in range(len(ID)):
    f = open("{0}/first_analysis/{1}/each_PSF_fit_qso.txt".format(ID[j],folder),"r")
    string = f.read()
    PSF_id = re.findall(r"by PSF(.*?):",string)
    Re = re.findall(r"R_sersic':(.*?),",string)
    Chisq = re.findall(r"redu_Chisq':(.*?),",string)
    
    Re = [float(i) for i in Re]
    Chisq = [float(i) for i in Chisq]
    
#    PSFs_dict = {}
#    for key in filt_info.keys():
#        if filt_info[key] == filt:
#            PSFs, _=pickle.load(open('{0}/{0}_PSFs_QSO'.format(key),'rb'))
#            PSFs_dict.update({'{0}'.format(key):PSFs})
    
    sort_Chisq = np.argsort(np.asarray(Chisq))
    count_n = 8
    Chisq_best = Chisq[sort_Chisq[0]]
    Chisq_last= Chisq[sort_Chisq[count_n-1]]
    inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
    weight = np.zeros(len(Chisq))
    for i in sort_Chisq[:count_n]:
        weight[i] = np.exp(-1/2. * (Chisq[i]-Chisq_best)/(Chisq_best* inf_alp))
    
    # =============================================================================
    # Weighting result
    # =============================================================================
    weighted_Re = np.sum(np.array(Re)*weight) / np.sum(weight)
    rms_Re = np.sqrt(np.sum((np.array(Re)-weighted_Re)**2*weight) / np.sum(weight))
    Re_results.append([round(weighted_Re,3), round(rms_Re,3)])
Re_results = np.asarray(Re_results)

fig, ax =  plt.subplots(figsize=(11,11))
plt.errorbar(WFC3_re[:,0], Re_results[:, 0], xerr= WFC3_re[:,1], yerr=Re_results[:,1],
            fmt='o', color='blue',ecolor='gray',zorder=1)
##plt.title('', fontsize=27)
#texts = []
#for i in range(len(ID)):
#    texts.append(ax.text(WFC3_re[i,0], Re_results[i, 0], ID[i], fontsize=12))
#adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))    

plt.xlabel("effective radius by WFC3 (arcsec)",fontsize=27)
plt.ylabel("effective radius by ACS (arcsec)",fontsize=27)
plt.tick_params(labelsize=20)
plt.xlim(0,1.1)
plt.ylim(0,1.1)
x = np.linspace(-0.5, 1.1, 21)
y = x
plt.plot(x,y,zorder=-1)
#plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.show()