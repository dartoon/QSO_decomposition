#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 21:46:04 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import re
import matplotlib

f = open("../analysis_Astrodrz/CID237/analysis/fit_result_each/each_PSF_fit_qso.txt","r")
#f = open("/Users/Dartoon/Astro/analysis_back_up/second_run_with_smaller_frame_size/CID3242/analysis/fit_result/each_PSF_fit_qso.txt","r")
#f = open("../analysis_Astrodrz/fit_summary_with_std**2.txt","r")

string = f.read()

# findall
labels = re.findall(r"PSF\d+", string)
#labels = re.findall(r"CID\d+", string)
S_n_list = re.findall(r"n_sersic':(.*?),",string)
Re = re.findall(r"R_sersic':(.*?),",string)
host_flux_ratio = re.findall(r"host_flux_ratio_percent':(.*?)}",string)
cmap = matplotlib.cm.get_cmap('viridis')
normalize = matplotlib.colors.Normalize(vmin=0.3, vmax=7.0)

S_n_list = [float(value) for value in S_n_list]

colors = [cmap(normalize(value)) for value in S_n_list]

Re = [float(i) for i in Re]
host_flux_ratio = [float(i) for i in host_flux_ratio]

fig, ax = plt.subplots(figsize=(10,6))
ax.scatter(Re, host_flux_ratio, color=colors)
plt.ylim(0, 100) 

cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
cbar.set_label('Sersic n',  fontsize=12)

for i in range(len(S_n_list)):
    ax.text(float(Re[i])+0.002, host_flux_ratio[i], labels[i], fontsize=12)

ax.set_xlabel('Sersic Reff')
ax.set_ylabel('Host flux ratio')
plt.show()
#host_flux_ratio = [float(value) for value in host_flux_ratio]
#Re = [float(value) for value in Re]
#plt.hist(host_flux_ratio, bins='auto')  # arguments are passed to np.histogram
#plt.title("Histogram of host flux ratio percent")
#plt.show()
#
#plt.hist(Re, bins='auto')  # arguments are passed to np.histogram
#plt.title("Histogram of Sersic Reff")
#plt.show()
#
#plt.hist(S_n_list, bins='auto')  # arguments are passed to np.histogram
#plt.title("Histogram of Sersic Index")
#plt.show()