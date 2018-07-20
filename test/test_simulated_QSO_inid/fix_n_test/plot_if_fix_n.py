#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 11:45:16 2018

@author: Dartoon

simple plot
"""

import numpy as np
import matplotlib.pylab as plt

fix_n_list = np.linspace(0.5,7,14)
red_Chisq_list = np.array([3.237, 3.219, 3.217, 3.217, 3.218,3.219, 3.219])
host_ratio_list = np.array([25.9, 29.7, 33.0, 35.9, 38.6, 41.0, 43.1])

host=plt.figure(figsize=(10,7))
ax=host.add_subplot(111)

ax.plot(fix_n_list, fix_n_list,'bo')
#plt.title('The offest if fix n, fix host/AGN center',  fontsize=25)
plt.xlabel('fix Sersic n value', fontsize=25)
plt.ylabel('reduced Chisq', fontsize=25)
plt.grid()
plt.show()