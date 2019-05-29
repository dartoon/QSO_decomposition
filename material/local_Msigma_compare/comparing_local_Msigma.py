#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 17:43:42 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
plt.figure(figsize=(13.5,12))
#%%The sample from Woo+ 2009 https://ui.adsabs.harvard.edu/abs/2010ApJ...716..269W/abstract
f1 = 'Woo_local_sample.txt'
woo_loc = np.loadtxt(f1)

#woo_loc = np.concatenate((woo_loc, np.zeros(len(woo_loc))))
woo_loc[:,3] = np.log10(woo_loc[:,2]+woo_loc[:,3]) - np.log10(woo_loc[:,2])  # error
woo_loc[:,2] = np.log10(woo_loc[:,2])
#woo_loc[:,4] = np.log10(woo_loc[:,2]+woo_loc[:,3]) - np.log10(woo_loc[:,2])  # error
plt.errorbar(woo_loc[:,2],woo_loc[:,0], xerr=woo_loc[:,3] ,yerr=woo_loc[:,1],fmt='.',color='black',markersize=15, label='Woo 2011')


#%%########input local elliptical by Kormendy ############
#Distance(Mpc) logMbulge error MBH low high order flag(0/reg, 1/1, 2/2)
f2 = 'Kormendy_elliptical.txt'
Kore_loc = np.loadtxt(f2)
#Kore_loc = Kore_loc[Kore_loc[:,-1]!=2]
kor_e = np.zeros([len(Kore_loc),8])
kor_e[:,0] = Kore_loc[:,-1]
kor_e[:,1]= Kore_loc[:,1]    #LgMstar
kor_e[:,2]= Kore_loc[:,2]  # error
kor_e[:,3]= np.log10(Kore_loc[:,3]) + Kore_loc[:,6]     #LgMBH
kor_e[:,4]= np.sqrt((np.log10(Kore_loc[:,3]) - np.log10(Kore_loc[:,4]))**2 +  (np.log10(Kore_loc[:,5]) - np.log10(Kore_loc[:,3]))**2) # error
kor_e[:,4][kor_e[:,4]==np.inf] =  (np.log10(Kore_loc[:,5]) - np.log10(Kore_loc[:,3]))[kor_e[:,4]==np.inf]
kor_e[:,5] = np.log10(Kore_loc[:,7])  # velocity dispersion
kor_e[:,6] = np.log10(Kore_loc[:,7]+Kore_loc[:,8]) - np.log10(Kore_loc[:,7])  # error
kor_e[:,7] = np.log10(Kore_loc[:,7]) - np.log10(Kore_loc[:,7]-Kore_loc[:,8]) # error
plt.errorbar(kor_e[:,5],kor_e[:,3], xerr=[kor_e[:,6], kor_e[:,7]] ,yerr=kor_e[:,4],fmt='.',color='red',markersize=15, label='Kormendy 2013 elliptical')

#%%########input local elliptical by Kormendy ############
#Distance(Mpc) logMbulge error MBH low high order flag(0/reg, 1/1, 2/2)
f3 = 'Kormendy_classical_bulges.txt'
Kocb_loc = np.loadtxt(f3)
#Kocb_loc = Kocb_loc[Kocb_loc[:,-1]!=2]
kocb = np.zeros([len(Kocb_loc),8])
kocb[:,0] = Kocb_loc[:,-1]
kocb[:,1]= Kocb_loc[:,1]    #LgMstar
kocb[:,2]= Kocb_loc[:,2]  # error
kocb[:,3]= np.log10(Kocb_loc[:,3]) + Kocb_loc[:,6]     #LgMBH
kocb[:,4]= np.sqrt((np.log10(Kocb_loc[:,3]) - np.log10(Kocb_loc[:,4]))**2 +  (np.log10(Kocb_loc[:,5]) - np.log10(Kocb_loc[:,3]))**2) # error
kocb[:,5] = np.log10(Kocb_loc[:,7])  # velocity dispersion
kocb[:,6] = np.log10(Kocb_loc[:,7]+Kocb_loc[:,8]) - np.log10(Kocb_loc[:,7])  # upper error
kocb[:,7] = np.log10(Kocb_loc[:,7]) - np.log10(Kocb_loc[:,7]-Kocb_loc[:,8])  # lower error
plt.errorbar(kocb[:,5],kocb[:,3], xerr=[kocb[:,6], kocb[:,7]] ,
             yerr=kocb[:,4],fmt='.',color='blue',markersize=15, label='Kormendy 2013 classical bulge')

#############################################################
m_ml0 = 3.6 #4.48,4.04, 3.6
def lnlike(theta, x, y, yerr):
#    m, b, sint= theta
    b, sint= theta
    m = m_ml0
    model = m * x + b
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf
import scipy.optimize as op

###################fitting f1 f2 with MCMC#########################
x0=woo_loc[:,2]
y0=woo_loc[:,0]
yerr0=(woo_loc[:,1]**2+woo_loc[:,3]**2)**0.5  # 0.2 is the uncertainty level for the L_R

nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [-0.275, 0.417], args=(x0, y0, yerr0))
b_ml0, sint_ml0= result["x"]
xl = np.linspace(1.5, 2.7, 10)
plt.plot(xl, m_ml0*xl+b_ml0, color="k", linewidth=4.0,zorder=-0.5)

###################fitting f1 f3 f4 with MCMC#########################
#np.concatenate((z_ss, z_b11, z_cosmos),axis=0)
x1=np.concatenate((kor_e[:,5], kocb[:,5]),axis=0)
y1=np.concatenate((kor_e[:,3], kocb[:,3]),axis=0)
yerr1=(np.concatenate((kor_e[:,6], kocb[:,6]),axis=0)**2+np.concatenate((kor_e[:,7], kocb[:,7]),axis=0)**2)**0.5  # 0.2 is the uncertainty level for the L_R
#x = bloc[:,1]
#y = bloc[:,3]
#yerr = bloc[:,2]
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [-1.868, 0.465], args=(x1, y1, yerr1))
b_ml1,sint_ml1= result["x"]
m_ml1 = m_ml0

#%%
plt.plot(xl, m_ml1*xl+b_ml1, color="red", linewidth=4.0,zorder=-0.5)

plt.text(1.8, 5.5, "log$(M_{BH}/M_{\odot})$=%s+%slog$(\sigma_*)$"%(round(b_ml0 ,2),round(m_ml0,2)),color='black',fontsize=20)
plt.text(1.6, 9.5, "log$(M_{BH}/M_{\odot})$=%s+%slog$(\sigma_*)$"%(round(b_ml1 ,2),round(m_ml1,2)),color='red',fontsize=20)

#plt.title("Comparsion of the $M_{BH}$ - $M_*$ relation at local universe",fontsize=25)
plt.xlabel("log($\sigma_*$)",fontsize=35)
plt.ylabel("log$(M_{BH}/M_{\odot})$",fontsize=35)
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)
plt.legend(prop={'size':20},loc=2)
#plt.savefig("local_MM_comparison.pdf")
plt.show()

print round(b_ml1-b_ml0,3)
#COMMENTS:
#1. The local MBH-sigma by Kormendy are higher than the Woo2011 (fix slope value).
#2. Meaning that the logf by Kormendy would be higher than the Woo2011, that is, the 
#MBH meaused based on the recipe by Kormendy would be higher, than a value of 0.37.
#3. So, the offset is still there even by Kormendy's local.