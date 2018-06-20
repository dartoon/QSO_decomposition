#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 10:25:14 2018

@author: dartoon
"""

import numpy as np
import matplotlib.pylab as plt

fc_n_relax = {'CID216': 21.727, 'CID454': 21.011, 'CID607': 21.850, 
               'CID1174':21.163, 'CID3242':21.878, #'CID3570':22.190,
               'LID1273':20.407, 'LID1538':20.926, 'XID2138':21.329,
               'XID2396':20.368} 

rc_n_relax = {'CID216': 21.634, 'CID454': 20.778, 'CID607': np.nan, 
               'CID1174':21.058, 'CID3242':21.968, #'CID3570': np.nan,
               'LID1273':20.427, 'LID1538':21.134, 'XID2138': np.nan,
               'XID2396':20.305} 


fc_n_4 = {'CID216': 21.883, 'CID454': 21.169, 'CID607': 21.957, 
               'CID1174':21.239, 'CID3242':21.887, #'CID3570':22.413,
               'LID1273':20.427, 'LID1538':20.922, 'XID2138':21.496,
               'XID2396':20.568}

rc_n_4 = {'CID216': 21.885, 'CID454': 20.932, 'CID607': np.nan, 
               'CID1174':21.242, 'CID3242':22.001, #'CID3570':np.nan,
               'LID1273':20.536, 'LID1538':21.171, 'XID2138':np.nan,
               'XID2396':20.407}


fc_n_1dot5 = {'CID216': 22.271, 'CID454': 21.524, 'CID607': 22.241, 
               'CID1174':21.552, 'CID3242':22.201, #'CID3570':22.984,
               'LID1273':20.756, 'LID1538':21.301, 'XID2138':21.822,
               'XID2396':20.964} 

rc_n_1dot5 = {'CID216': 22.286, 'CID454': 21.551, 'CID607': np.nan, 
               'CID1174':21.554, 'CID3242':22.413, #'CID3570':np.nan,
               'LID1273':20.784, 'LID1538':21.460, 'XID2138':np.nan,
               'XID2396':20.964} 

diff_fc_n4 = dict()
diff_rc_n4 = dict()
diff_fc_n1dot5 = dict() 
diff_rc_n1dot5 = dict() 
for key in fc_n_relax:
    diff_fc_n4[key] = fc_n_4[key] - fc_n_relax[key]
    diff_rc_n4[key] = rc_n_4[key]- rc_n_relax[key]
    diff_fc_n1dot5[key] = fc_n_1dot5[key] - fc_n_relax[key]
    diff_rc_n1dot5[key] = rc_n_1dot5[key] - rc_n_relax[key]

fr = input('fix center of relax center?:1 fix, 2, relax:\n')

if fr == 1:
    lists_n4 = sorted(diff_fc_n4.items())
    lists_n15 = sorted(diff_fc_n1dot5.items())
elif fr == 2:
    lists_n4 = sorted(diff_rc_n4.items())
    lists_n15 = sorted(diff_rc_n1dot5.items())
target, diff_n4 = zip(*lists_n4)
_, diff_n15 = zip(*lists_n15)

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
host=plt.figure(figsize=(10,7))
ax=host.add_subplot(111)

patterns = ['.', '+', 'o', 'v', '^', '<', '>', '*', 'p','P']
colors =  ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'g', 'coral', 'darkblue']
for i in range(len(diff_n15)):
    if diff_n15[i] > -10:  ## Don't know how to define nan...
        plt.scatter(1.5, diff_n15[i], label = target[i],c=colors[i], marker=patterns[i],s=100)
        plt.scatter(4. , diff_n4[i], marker=patterns[i],c=colors[i],s=100)
ax.get_xaxis().set_visible(False)
#ax.get_yaxis().set_visible(False)
plt.axis([0,6,-0.2,1])
plt.tick_params(labelsize=25)
plt.legend(loc=1,prop={'size':15})
plt.text(1,0,'fix n=1.5', color='b', fontsize=25)
plt.text(3.5,0.4,'fix n=4', color='b', fontsize=25)
plt.ylabel('$mag_{fix\ n} - mag_{relax\ n}$', fontsize=25)
if fr == 1:
    plt.title('The offest if fix n, fix host/AGN center',  fontsize=25)
elif fr ==2:
    plt.title('The offest if fix n, relax host/AGN center',  fontsize=25)
plt.grid()
plt.show()