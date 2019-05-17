import numpy as np
import matplotlib.pyplot as plt

############### with evolution-corrected or not? #################
#select = 1
host = 1

#select= 0 #input('with evolution-corrected or not??? 0 = no;   1= yes, with dmag using dmag.py:')
#if select == 0:
#   dm=0
#if select == 1:
dm= 1
import sys
sys.path.insert(0,'../py_tools')
from dmag import pass_dmag
########input Park local data ############
f1 ='data/parklocal'
Pklc = np.loadtxt(f1)
#print np.shape(Pklc), len(Pklc)
ploc=np.zeros([len(Pklc),4])

#host= input('with bulge or total relation??? 0 = bulge;   1= total:')
if host == 0:
   ploc[:,1]= Pklc[:,6]    #LgV_sph
if host == 1:
   ploc[:,1]= Pklc[:,5]    #LgV_total
ploc[:,0]= Pklc[:,0]    #redshift
ploc[:,1]=0.4*(4.61+0.46-4.83)+ploc[:,1]-0.4*dm*pass_dmag(ploc[:,0])  #Change to LgR_sph, -0.4 is in L, be fainter
ploc[:,2]= Pklc[:,8]-0.03    #LogBHmass
ploc[:,3]= Pklc[:,9]    #delta mass
loc=np.zeros([len(ploc),3])
loc[:,0]=ploc[:,1]
loc[:,1]=ploc[:,2]
loc[:,2]=ploc[:,3]
#############################################################
#############################################################
###################fitting with MCMC#########################
x=loc[:,0]
y=loc[:,1]
yerr=loc[:,2]
def lnlike(theta, x, y, yerr):
    m, b, sint= theta
    model = m * x + b
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [1.036, -1.947, 0.3], args=(x, y, yerr))
m_ml, b_ml,sint_ml= result["x"]
#print m_ml, b_ml, sint_ml, "ka=",lnlike(theta=[m_ml, b_ml, sint_ml],x=loc[:,0], y=loc[:,1], yerr=loc[:,2])

#m_ml, b_ml = 1.11, 7.40-1.11*10  # Park's first use

xp = np.array([5, 13])
#plt.plot(xp, m_ml*xp+b_ml, 'r-')
def lnprior(theta):
    m, b, sint	 = theta
    if -5.0 < m < 5 and -10 < b < 10.0 and 0 < sint < 10:
        return 0.0
    return -np.inf
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)
ndim, nwalkers = 3, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

m_mid, b_mid, sint_mid =np.percentile(samples, 50,axis=0)
#print "lnlike=",lnlike(theta=[m_mid, b_mid, sint_mid],x=loc[:,0], y=loc[:,1], yerr=loc[:,2])
xl = np.linspace(-0.9, 13, 100)
#plt.plot(np.log10(1+xl), xl*0, color="black", linewidth=4.0,zorder=-40)


######################
#Klc=plt.errorbar(loc1[:,0]*0,loc1[:,1]-(m_mid*loc1[:,0]+b_mid),yerr=0.001,fmt='.',color='gray',markersize=10)
#Plc=plt.errorbar(loc2[:,0]*0+0.01,loc2[:,1]-(m_mid*loc2[:,0]+b_mid),yerr=0.001,fmt='.',color='black',markersize=10)
Pkc=plt.errorbar(np.log10(ploc[:,0]+1),loc[:,1]-(m_ml*loc[:,0]+b_ml),yerr=ploc[:,3],fmt='.',color='gray',markersize=10)

ty=xl*0
ty1=xl*0+np.std(loc[:,1]-(m_ml*loc[:,0]+b_ml))
ty2=xl*0-np.std(loc[:,1]-(m_ml*loc[:,0]+b_ml))
plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)

######################
