#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 21:30:15 2018

@author: dartoon

To derive the individual image in order to Swarp
"""
import os
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import numpy as np
from subprocess import call

target=raw_input("name of the target:\n")
folder = '../Cycle25data/{0}/'.format(target)
suffix='_ima.fits'

if os.path.exists('{0}/'.format(target)) is False:
    os.mkdir("{0}".format(target))
    os.mkdir("{0}/swarp".format(target))
    

for root, dirs, files in os.walk(folder):
        files=files
        
for i in range(len(files)):
    if suffix in files[i]:
        print files[i]
        filename=folder + files[i]
        image = pyfits.open(filename)#[1].data.copy()
        input_header= image[1].header
        flux= image[1].data.copy()
#        print input_header
        plt.imshow((flux),origin='lower',vmin=-0.5,vmax=60)
        plt.colorbar()
        plt.show()
        hdu = pyfits.PrimaryHDU(flux,header=input_header)
        hdu.writeto('{0}/swarp/{1}'.format(target,files[i]),overwrite=True)
cmd_str1=("cp ../temp_swarp/* {0}/swarp/".format(target))
#cmd_str2=('./auto_swap.sh')
call(cmd_str1, shell=True)
#call("cd {0}/swarp/".format(target), shell=True)
#call("./auto_swarp.sh".format(target), shell=True)
#in the shell use "./auto_swarp.sh" to generate the final image.