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
f1 ='data/Bentz_MM.txt'
Ben_loc = np.loadtxt(f1)
ben = np.zeros([len(Ben_loc),5])
ben[:,0] = Ben_loc[:,0]
ben[:,1]= Ben_loc[:,1]    #LgMstar
ben[:,2]= Ben_loc[:,2]  # error
ben[:,3]= Ben_loc[:,3] + 0.08     #LgMBH, the f in Bentz is log(f) = 0.63, which is lower then ours 0.71 by 0.08
ben[:,4]= (abs(Ben_loc[:,4])+ abs(Ben_loc[:,5]))/2 #err
#plt.errorbar(ben[:,1],ben[:,3], xerr=ben[:,2] ,
#             yerr=ben[:,4],fmt='.',color='blue',markersize=15, label='Kormendy 2013 classical bulge')

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
solve_z = np.vectorize(solve_z)  #Calculate the redshift using the distance
ben[:,0] =  solve_z(ben[:,0])
# +  np.log10(Haring04[:,4] * 10 ** Haring04[:,5]))/2 #Sigma LgMBH
#plt.errorbar(ben[:,1],ben[:,3], xerr=ben[:,2] ,yerr=ben[:,4],fmt='.',color='black',markersize=15)
#plt.plot(ben[:,1],ben[:,3], '.',color='black',markersize=15)
#Hkc=mlines.Line2D([], [], color='black', ls='', marker='.', markersize=15)

#############################################################
###################fitting with MCMC#########################
x=ben[:,1]
y=ben[:,3]
yerr=(ben[:,2]**2+ ben[:,4]**2)**0.5  # 0.2 is the uncertainty level for the L_R

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
    plt.plot(np.log10(ben[:,0]+1),
             10**(ben[:,3]-ben[:,1]), '.',color='black',markersize=10)
    
elif style ==1:
    xl = np.linspace(-0.9, 13, 100)
    plt.errorbar(np.log10(ben[:,0]+1),
                 ben[:,3]-(m_ml*ben[:,1]+b_ml),yerr=(ben[:,2]**2 + ben[:,4]**2)**0.5 ,fmt='.',color='black',markersize=10)
    ty=xl*0
    ty1=xl*0+np.std(y-(m_ml*x+b_ml))
    ty2=xl*0-np.std(y-(m_ml*x+b_ml))
    plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)

ben=mlines.Line2D([], [], color='black', ls='', marker='.', markersize=15)
######################
##%%
##calcualte the mean offset for the local:
#y_local=np.append(ben[:,3]-(m_ml*ben[:,1]+b_ml),kor_e[:,3]-(m_ml*kor_e[:,1]+b_ml))
#y_local_err = np.append((ben[:,2]**2 + ben[:,4]**2)**0.5, (kor_e[:,2]**2 + kor_e[:,4]**2)**0.5)
#weighted_offset = np.sum(np.asarray(y_local)*y_local_err) / np.sum(y_local_err)                              
#rms_offset = np.sqrt(np.sum((np.asarray(y_local)-weighted_offset)**2*y_local_err) / np.sum(y_local_err))
#print weighted_offset, rms_offset