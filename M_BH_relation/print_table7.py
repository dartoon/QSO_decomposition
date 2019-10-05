#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 15:50:27 2018

@author: Dartoon
"""
import numpy as np
np.set_printoptions(precision=4)
import matplotlib.pyplot as plt
import matplotlib as mat
import matplotlib.lines as mlines
from matplotlib import colors
mat.rcParams['font.family'] = 'STIXGeneral'


import matplotlib as mpl
mpl.rc('image', cmap='jet')

import sys
sys.path.insert(0,'../py_tools')
########## input local data ####
#==============================================================================
# The seleting for dm and host_total and dmag are in this local
#==============================================================================
from local_MM_vz import *
#==============================================================================
#input SS13 and B11 data 
#==============================================================================
################ bulge or total relaiton? #################
####### input SS13 data ####
#######in AB system, V band#######
inp_SS13 = 1
f0 ='data/SS13_MM.txt'
ss = np.loadtxt(f0)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
if inp_SS13 ==1:
    plt.scatter(ss[:,1],ss[:,2],c='darkseagreen',marker="^",s=180,zorder=100, edgecolors='white')
  
inp_b11= 1
f1 ='data/B11_MM.txt'
b11 = np.loadtxt(f1)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
if inp_b11 ==1:
    plt.scatter(b11[:,1],b11[:,2],c='darkseagreen',marker="^",s=180,zorder=100, edgecolors='white')

inp_Cis= 1
f2 = 'data/Cisternas_data.txt'
cis11 = np.loadtxt(f2)  #0 redshift;
if inp_Cis ==1:
    plt.scatter(cis11[:,2],cis11[:,1],c='darkseagreen',marker="^",s=180,zorder=100, edgecolors='white')
plt.close()


#%%Print 
hlocname = ['M87', 'NGC1068', 'NGC3379', 'NGC4374', 'NGC4261', 'NGC6251', 'NGC7052',
            'NGC4742', 'NGC821', 'IC1459', 'M31', 'M32', 'NGC1023', 'NGC2778', 'NGC3115',
            'NGC3245', 'NGC3377', 'NGC3384', 'NGC3608', 'NGC4291', 'NGC4342', 'NGC4473',
            'NGC4564', 'NGC4594', 'NGC4649', 'NGC4697', 'NGC5845', 'NGC7332', 'NGC7457',
            'Milky Way']

for i in range(len(hloc)):
    print hlocname[i] + " & {0:.4f} & {1:.2f}$\pm${2:.2f} & {3:.2f}$\pm${4:.1f} \\\\".format(hloc[i,0], hloc[i,1], hloc[i,2], hloc[i,3], hloc[i,4])
    
#%%
blocname = ['0121$-$0102', '0206$-$0017', '0353$-$0623', '0802+3104', '0846+2522', '1042+0414',
            '1043+1105', '1049+2451', '1101+1102', '1116+4123', '1144+3653', '1210+3820',
            '1250$-$0249', '1323+2701', '1355+3834', '1405$-$0259', '1419+0754', '1434+4839',
            '1535+5754', '1545+1709', '1554+3238', '1557+0830', '1605+3305', '1606+3324', '1611+5211']
for i in range(len(bloc)):
    print blocname[i] + " & {0:.4f} & {1:.2f}$\pm${2:.2f} & {3:.2f}$\pm${4:.1f} \\\\".format(bloc[i,0], bloc[i,1], bloc[i,2], bloc[i,3], bloc[i,4])

#%%

b11name = ['J033252$-$275119', 'J033243$-$274914', 'J033239$-$274601', 'J033226$-$274035', 'J033225$-$274218', 'J033210$-$274414',
           'J033200$-$274319', 'J033229$-$274529', 'J123553$+$621037', 'J123618$+$621115', 'J123707$+$622147']
for i in range(len(b11)):
    print b11name[i] + " & {0:.3f} & {1:.2f}$\pm${2:.1f} & {3:.2f}$\pm${4:.1f} \\\\".format(b11[i,0],  b11[i,1], 0.2, b11[i,2], 0.4)

#%%
ssname = ['158', '170', '271', '273', '305', '333', '339', '348', '379', '413', '417', '465', '516', '540', '597', '712']  #ID 250 and 375 is deleted
for i in range(len(ss)):
    print  ssname[i] + " & {0:.3f} & {1:.2f}$\pm${2:.1f} & {3:.2f}$\pm${4:.1f} \\\\".format(ss[i,0], ss[i,1], 0.2, ss[i,2], 0.4)

#%%
cis11name = ['COSMOS J095817.54+021938.5','SDSS J095819.88+022903.6','COSMOS J095831.65+024901.6','COSMOS J095840.61+020426.6',
             'COSMOS J095845.80+024634.0','SDSS J095902.76+021906.5','COSMOS J095909.53+021916.5','COSMOS J095928.31+022106.9',
             'COSMOS J100002.21+021631.8','SDSS J100012.91+023522.8','COSMOS J100014.55+023852.7','COSMOS J100017.54+020012.6',
             'SDSS J100025.25+015852.2','COSMOS J100028.63+025112.7','COSMOS J100029.69+022129.7','COSMOS J100033.38+015237.2',
             'COSMOS J100033.49+013811.6','COSMOS J100037.29+024950.6','SDSS J100043.15+020637.2','COSMOS J100046.72+020404.5',
             'COSMOS J100058.71+022556.2','COSMOS J100118.52+015543.0','COSMOS J100141.09+021300.0','COSMOS J100146.49+020256.7',
             'COSMOS J100202.22+024157.8','COSMOS J100205.03+023731.5','COSMOS J100212.11+014232.4','COSMOS J100218.32+021053.1',
             'COSMOS J100230.06+014810.4','COSMOS J100230.65+024427.6','SDSS J100232.13+023537.3','COSMOS J100243.96+023428.6']
for i in range(len(cis11)):
    print  cis11name[i] + " & {0:.2f} & {1:.2f}$\pm${2:.1f} & {3:.2f}$\pm${4:.1f} \\\\".format(cis11[i,0], cis11[i,2], 0.2, cis11[i,1], 0.4)
