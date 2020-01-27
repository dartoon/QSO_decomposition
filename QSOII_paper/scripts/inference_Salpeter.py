#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 13:59:04 2019

@author: 
    
My inference of in Salpeter version.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy
import sys
sys.path.insert(0,'../../py_tools')
from load_result import load_MBH

#Print the ID, z, M_r, L_r, M_star , MBH, in SALPETER version.
tab_list = ['CID1174', 'CID1281', 'CID206', 'CID216', 'CID237', 'CID255', 'CID3242', 'CID3570',
            'CID452', 'CID454', 'CID50', 'CID543', 'CID597', 'CID607', 'CID70', 'LID1273',
            'LID1538', 'LID360', 'XID2138', 'XID2202', 'XID2396', 'CDFS-1', 'CDFS-229',
            'CDFS-321', 'CDFS-724', 'ECDFS-358', 'SXDS-X1136', 'SXDS-X50', 'SXDS-X717', 'SXDS-X735', 'SXDS-X763', 'SXDS-X969']
ext_ID = {'XID2202':'LID1622', 'XID2138':'LID1820', 'XID2396':'LID1878', 'CDFS-321':'ECDFS-321'}
tab_sub_list = copy.deepcopy(tab_list)
for i in range(len(tab_sub_list)):
    if tab_sub_list[i] in ext_ID.keys():
        tab_sub_list[i] = ext_ID[tab_sub_list[i]]
MBH = load_MBH(tab_list,tab_sub_list, if_reportHb=0, folder = '../../')   

from load_result import load_zs, load_host_p, load_err
zs = load_zs(tab_list)

Lr, M_star, M_r = load_host_p(tab_list, folder = '../../')  # M_star by Chabrier 
#transfer M_start to Salpeter:
M_star = np.log10(10**M_star / 0.54 * 0.91)

#for i in range(len(tab_list)):
##    print tab_list[i], zs[i], M_r[i], M_star[i], MBH[i]
#    print tab_list[i], round(M_r[i],3), round(MBH[i],3),  round(M_star[i],3)
    
M_star_err = load_err(prop='Mstar', ID=tab_list)
M_star_err = (abs(M_star_err[:,0]) + abs(M_star_err[:,1]))/2
    
#%%Infer the intrinsic scatter of the data.
import linmix
x = M_star
xsig = M_star_err
y = MBH
ysig = np.ones(len(MBH))*0.4

lm = linmix.LinMix(x, y, xsig=xsig, ysig=ysig, K=2)
lm.run_mcmc(silent=True)

burnin=5000
fig,ax = plt.subplots(1,1,figsize=(10,5))

# These lines will plot 'nlines' of the instances from the linmix_err MCMC chains
xvals = np.linspace(ax.get_xlim()[0],ax.get_xlim()[1])
lmchain = lm.chain
chain = {
    'alpha': lm.chain['alpha'][burnin:],
    'beta': lm.chain['beta'][burnin:],
    'sigsqr': lm.chain['sigsqr'][burnin:],
    'corr': lm.chain['corr'][burnin:]
}
nlines = 1000
##for i in np.arange(0,len(chain['alpha']),int(len(chain['alpha'])/nlines)):
##    plt.plot(xvals, chain['alpha'][i] + chain['beta'][i]*xvals, color='g', alpha=0.01)
##plt.show()
#
#print("{}, {}".format(lm.chain['alpha'].mean(), lm.chain['alpha'].std()))
#print("{}, {}".format(lm.chain['beta'].mean(), lm.chain['beta'].std()))
#print("{}, {}".format(lm.chain['sigsqr'].mean(), lm.chain['sigsqr'].std()))
## The value of the intrinsic scatter:
sig_int = np.sqrt(lm.chain['sigsqr'].mean())
print sig_int
#Plot the result:

tot_scatter = []
plt.errorbar(M_star,MBH, xerr=M_star_err, yerr=0.4, color='blue',ecolor='orange', fmt='.',zorder=-500,markersize=1)
for i in range(0, len(lm.chain)):
    xs = np.arange(8,13)
    ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
#    plt.plot(xs, ys, color='r', alpha=0.02)
    tot_scatter.append(np.sqrt(np.mean((y - (lm.chain[i]['alpha'] + x * lm.chain[i]['beta']))**2)))
tot_scatter = np.array(tot_scatter)
   
alpha = lm.chain['alpha'].mean()
beta = lm.chain['beta'].mean()
ys = alpha + xs * beta
plt.plot(xs, ys, color='k')

plt.xlim(9,12.5)
plt.ylim(6.0,10)
plt.show()

#%%Fitting by myself.
x = M_star
xsig = M_star_err
y = MBH
ysig = np.ones(len(MBH))*0.4
m = 1

def lnlike(theta, x, y):
    m, b, sint= theta
    model = m * x + b
    sig = ((m*xsig)**2. + ysig**2.) ** 0.5 #/2
    sigma2 = (sig**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
#      print -abs((np.sum((y-model)**2/sigma2)) - 32)
#      return -abs((np.sum((y-model)**2/sigma2)) - 32)
    else:
      return -np.inf

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [0.5, 2.79, 0.1], args=(x, y))
m_ml, b_ml,sint_ml= result["x"]
print result["x"]

def lnprior(theta):
    m, b, sint	 = theta
    if -5.0 < m < 5 and -10 < b < 10.0 and 0 < sint < 10:
        return 0.0
    return -np.inf
def lnprob(theta, x, y):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y)
ndim, nwalkers = 3, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y))
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

m_mid, b_mid, sint_mid =np.percentile(samples, 50,axis=0)

plt.errorbar(x, y, xerr=M_star_err, yerr=0.4, color='blue',ecolor='orange', fmt='.',zorder=-500,markersize=1)
xs = np.arange(8,13)
ys = b_ml + xs * m_ml
plt.plot(xs, ys, color='k')
plt.xlim(9,12.5)
plt.ylim(6.0,10)
plt.show()

import corner
fig = corner.corner(samples, labels=["$m$", "$b$", "$sint$"],
                       quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})
