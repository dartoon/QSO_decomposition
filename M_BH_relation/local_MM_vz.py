import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
############### with evolution-corrected or not? #################
#select = 1
#host = 1

import sys
sys.path.insert(0,'../py_tools')
#from dmag import pass_dmag

########input 25 local by Bennert++2011 ############
f1 ='data/Bennert+2011_local.txt'
local25 = np.loadtxt(f1)
bloc = np.zeros([len(local25),5])
bloc[:,0]= local25[:,1]    #Redshift
bloc[:,1]= local25[:,5]    #LgMstar
bloc[:,2]= local25[:,6]    #Sigma LgMstar
bloc[:,3]= local25[:,7]    #LgMBH
bloc[:,4]= local25[:,8]    #Sigma LgMBH
#plt.errorbar(bloc[:,1],bloc[:,3], xerr=bloc[:,2] ,yerr=bloc[:,4],fmt='.',color='gray',markersize=15)
#plt.plot(bloc[:,1],bloc[:,3], '.',color='gray',markersize=15)
#Bkc=mlines.Line2D([], [], color='gray', ls='', marker='.', markersize=15)

########input 30 local by Haring 04 ############
f2 ='data/Haring04.txt'
Haring04 = np.loadtxt(f2)
hloc = np.zeros([len(Haring04),5])
hloc[:,1]= np.log10(Haring04[:,0] * 10 ** Haring04[:,1])    #LgMstar
hloc[:,2]= 0.18  # Mention in the Haring04 paper, note in fig 2
hloc[:,3]= np.log10(Haring04[:,2] * 10 ** Haring04[:,5])    #LgMBH
hloc[:,4]= (abs(np.log10(Haring04[:,2] + Haring04[:,3]) - np.log10(Haring04[:,2])) + abs(np.log10(Haring04[:,2] - Haring04[:,4]) - np.log10(Haring04[:,4])))/2

from scipy import integrate
from scipy.optimize import fsolve
def Ez(z,om):
    '''
    Actually, this is the reverse of the E(z)
    '''
    w = -1
    return   1/np.sqrt(om*(1+z)**3+(1-om)*(1+z)**(3*(1+w)))   
def r(z,om):
    return integrate.quad(Ez, 0, z, args=(om))[0]
vec_r=np.vectorize(r)
def dl(z,om=0.3,h0=70):
    """
    Calculate the luminosity distance.
    """
    c=299790.
    dl_distance = (1+z) * c/h0 * vec_r(z,om=om)
    return dl_distance
#print dl(np.array([1,2]))
def solve_z(lum_dis, om=0.3, h0=70):
    func = lambda z : (lum_dis-dl(z,om=om, h0=h0))
    zs = fsolve(func,2)
    return zs[0]
Mpc = np.array([16.1, 15.0, 10.6, 18.4, 31.6, 106.0, 58.7, 15.5, 24.1, 29.2, 0.76, 0.81, 11.4, 22.9, 9.7, 20.9, 11.2, 11.6, 23.0, 26.2, 15.3, 15.7, 15.0, 9.8, 16.8, 11.7, 25.9, 23.0, 13.2, 0.1])
solve_z = np.vectorize(solve_z)  #Calculate the redshift using the distance
zs =  solve_z(Mpc)
hloc[:,-1] = zs
# +  np.log10(Haring04[:,4] * 10 ** Haring04[:,5]))/2 #Sigma LgMBH
#plt.errorbar(hloc[:,1],hloc[:,3], xerr=hloc[:,2] ,yerr=hloc[:,4],fmt='.',color='black',markersize=15)
#plt.plot(hloc[:,1],hloc[:,3], '.',color='black',markersize=15)
#Hkc=mlines.Line2D([], [], color='black', ls='', marker='.', markersize=15)

#############################################################
###################fitting with MCMC#########################
x=np.append(bloc[:,1], hloc[:,1])
y=np.append(bloc[:,3], hloc[:,3])
yerr=(np.append(bloc[:,2], hloc[:,2])**2+np.append(bloc[:,4], hloc[:,4])**2)**0.5  # 0.2 is the uncertainty level for the L_R
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

#m_ml, b_ml = 1.12, -4.12  # Park's first use
#m_ml, b_ml = 1.02, -2.91  # Park's inference from the two sample

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
#plt.plot(np.log10(1+xl), xl*0, color="black", linewidth=4.0,zorder=-40)


######################
style = input("0 as SS13, 1 as delta(logMBH):\n")
if style == 0:
    plt.plot(np.log10(hloc[:,0]+1),
             10**(hloc[:,3]-hloc[:,1]), '.',color='black',markersize=10)
    
    plt.plot(np.log10(bloc[:,0]+1),
                 10**(bloc[:,3]-bloc[:,1]),'.',color='gray',markersize=10)
elif style ==1:
    xl = np.linspace(-0.9, 13, 100)
    plt.errorbar(np.log10(hloc[:,0]+1),
                 hloc[:,3]-(m_ml*hloc[:,1]+b_ml),yerr=(hloc[:,2]**2 + hloc[:,4]**2)**0.5 ,fmt='.',color='black',markersize=10)
    plt.errorbar(np.log10(bloc[:,0]+1),
                 bloc[:,3]-(m_ml*bloc[:,1]+b_ml),yerr=(bloc[:,2]**2 + bloc[:,4]**2)**0.5 ,fmt='.',color='gray',markersize=10)
    ty=xl*0
    ty1=xl*0+np.std(y-(m_ml*x+b_ml))
    ty2=xl*0-np.std(y-(m_ml*x+b_ml))
    plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)

Bkc=mlines.Line2D([], [], color='gray', ls='', marker='.', markersize=15)
Hkc=mlines.Line2D([], [], color='black', ls='', marker='.', markersize=15)
######################
##%%
##calcualte the mean offset for the local:
#y_local=np.append(hloc[:,3]-(m_ml*hloc[:,1]+b_ml),bloc[:,3]-(m_ml*bloc[:,1]+b_ml))
#y_local_err = np.append((hloc[:,2]**2 + hloc[:,4]**2)**0.5, (bloc[:,2]**2 + bloc[:,4]**2)**0.5)
#weighted_offset = np.sum(np.asarray(y_local)*y_local_err) / np.sum(y_local_err)                              
#rms_offset = np.sqrt(np.sum((np.asarray(y_local)-weighted_offset)**2*y_local_err) / np.sum(y_local_err))
#print weighted_offset, rms_offset