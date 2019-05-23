#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 16:26:54 2019

@author: Dartoon

Comparing the local sample between HR+Vardha and Kormendy
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
plt.figure(figsize=(13.5,12))
########input 25 local by Bennert++2011 ############
f1 ='data/Bennert+2011_local.txt'
local25 = np.loadtxt(f1)
bloc = np.zeros([len(local25),5])
bloc[:,0]= local25[:,1]    #Redshift
bloc[:,1]= local25[:,5]    #LgMstar
bloc[:,2]= local25[:,6]    #Sigma LgMstar
bloc[:,3]= local25[:,7]    #LgMBH
bloc[:,4]= local25[:,8]    #Sigma LgMBH
plt.errorbar(bloc[:,1],bloc[:,3], xerr=bloc[:,2] ,yerr=bloc[:,4],fmt='.',color='gray',markersize=15, label='Bennert 2011')

########input 30 local by Haring 04 ############
f2 ='data/Haring04.txt'
Haring04 = np.loadtxt(f2)
hloc = np.zeros([len(Haring04),5])
hloc[:,1]= np.log10(Haring04[:,0] * 10 ** Haring04[:,1])    #LgMstar
hloc[:,2]= 0.18  # Mention in the Haring04 paper, note in fig 2
hloc[:,3]= np.log10(Haring04[:,2] * 10 ** Haring04[:,5])    #LgMBH
hloc[:,4]= (abs(np.log10(Haring04[:,2] + Haring04[:,3]) - np.log10(Haring04[:,2])) + abs(np.log10(Haring04[:,2] - Haring04[:,4]) - np.log10(Haring04[:,4])))/2
plt.errorbar(hloc[:,1],hloc[:,3], xerr=hloc[:,2] ,yerr=hloc[:,4],fmt='.',color='black',markersize=15, label='Haring 2004')

#%%########input local elliptical by Kormendy ############
#Distance(Mpc) logMbulge error MBH low high order flag(0/reg, 1/1, 2/2)
f3 = 'data/Kormendy_elliptical.txt'
Kore_loc = np.loadtxt(f3)
#Kore_loc = Kore_loc[Kore_loc[:,-1]!=2]
kor_e = np.zeros([len(Kore_loc),5])
kor_e[:,0] = Kore_loc[:,-1]
kor_e[:,1]= Kore_loc[:,1]    #LgMstar
kor_e[:,2]= Kore_loc[:,2]  # error
kor_e[:,3]= np.log10(Kore_loc[:,3]) + Kore_loc[:,6]     #LgMBH
kor_e[:,4]= np.sqrt((np.log10(Kore_loc[:,3]) - np.log10(Kore_loc[:,4]))**2 +  (np.log10(Kore_loc[:,5]) - np.log10(Kore_loc[:,3]))**2) # error
kor_e[:,4][kor_e[:,4]==np.inf] =  (np.log10(Kore_loc[:,5]) - np.log10(Kore_loc[:,3]))[kor_e[:,4]==np.inf]
plt.errorbar(kor_e[:,1],kor_e[:,3], xerr=kor_e[:,2] ,yerr=kor_e[:,4],fmt='.',color='red',markersize=15, label='Kormendy 2013 elliptical')

#%%########input local classic bulge by Kormendy ############
#Distance(Mpc) logMbulge error MBH low high order flag(0/reg, 1/1, 2/2)
f4 = 'data/Kormendy_classical_bulges.txt'
Kocb_loc = np.loadtxt(f4)
#Kocb_loc = Kocb_loc[Kocb_loc[:,-1]!=2]
kocb_e = np.zeros([len(Kocb_loc),5])
kocb_e[:,0] = Kocb_loc[:,-1]
kocb_e[:,1]= Kocb_loc[:,1]    #LgMstar
kocb_e[:,2]= Kocb_loc[:,2]  # error
kocb_e[:,3]= np.log10(Kocb_loc[:,3]) + Kocb_loc[:,6]     #LgMBH
kocb_e[:,4]= np.sqrt((np.log10(Kocb_loc[:,3]) - np.log10(Kocb_loc[:,4]))**2 +  (np.log10(Kocb_loc[:,5]) - np.log10(Kocb_loc[:,3]))**2) # error
plt.errorbar(kocb_e[:,1],kocb_e[:,3], xerr=kocb_e[:,2] ,
             yerr=kocb_e[:,4],fmt='.',color='blue',markersize=15, label='Kormendy 2013 classical bulge')


#############################################################
def lnlike(theta, x, y, yerr):
    m, b, sint= theta
    model = m * x + b
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf
import scipy.optimize as op
###################fitting f1 f2 with MCMC#########################
x0=np.append(bloc[:,1], hloc[:,1])
y0=np.append(bloc[:,3], hloc[:,3])
yerr0=(np.append(bloc[:,2], hloc[:,2])**2+np.append(bloc[:,4], hloc[:,4])**2)**0.5  # 0.2 is the uncertainty level for the L_R

nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [0.93027905, -1.95536508, 0.35], args=(x0, y0, yerr0))
m_ml0, b_ml0, sint_ml0= result["x"]
xl = np.linspace(5, 13, 100)
plt.plot(xl, m_ml0*xl+b_ml0, color="k", linewidth=4.0,zorder=-0.5)

###################fitting f1 f3 f4 with MCMC#########################
#np.concatenate((z_ss, z_b11, z_cosmos),axis=0)
x1=np.concatenate((bloc[:,1], kor_e[:,1], kocb_e[:,1]),axis=0)
y1=np.concatenate((bloc[:,3], kor_e[:,3], kocb_e[:,3]),axis=0)
yerr1=(np.concatenate((bloc[:,2], kor_e[:,2], kocb_e[:,2]),axis=0)**2+np.concatenate((bloc[:,4], kor_e[:,4], kocb_e[:,4]),axis=0)**2)**0.5  # 0.2 is the uncertainty level for the L_R
#x = bloc[:,1]
#y = bloc[:,3]
#yerr = bloc[:,2]
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [0.93027905, -1.95536508, 0.35], args=(x1, y1, yerr1))
m_ml1, b_ml1,sint_ml1= result["x"]
plt.plot(xl, m_ml1*xl+b_ml1, color="red", linewidth=4.0,zorder=-0.5)


plt.text(10.5, 6.5, "log$(M_{BH}/10^{7}M_{\odot})$=%s+%slog$(M_*/10^{10}M_{\odot})$"%(round(b_ml0+m_ml0*10-7,2),round(m_ml0,2)),color='black',fontsize=20)
plt.text(9.20, 9.5, "log$(M_{BH}/10^{7}M_{\odot})$=%s+%slog$(M_*/10^{10}M_{\odot})$"%(round(b_ml1+m_ml1*10-7,2),round(m_ml1,2)),color='red',fontsize=20)


plt.xlim(9,12.5)
plt.ylim(6.0,11)
#plt.title("Comparsion of the $M_{BH}$ - $M_*$ relation at local universe",fontsize=25)
plt.xlabel("log$(M_*/M_{\odot})$",fontsize=35)
plt.ylabel("log$(M_{BH}/M_{\odot})$",fontsize=35)
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)
plt.legend(prop={'size':20},loc=2)
#plt.savefig("local_MM_comparison.pdf")
plt.show()