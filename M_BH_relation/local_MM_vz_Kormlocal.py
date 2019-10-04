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
f1 ='data/Kormendy_classical_bulges.txt'
Kocb_loc = np.loadtxt(f1)
kocb = np.zeros([len(Kocb_loc),5])
kocb[:,0] = Kocb_loc[:,0]
kocb[:,1]= Kocb_loc[:,1]    #LgMstar
kocb[:,2]= Kocb_loc[:,2]  # error
kocb[:,3]= np.log10(Kocb_loc[:,3]) + Kocb_loc[:,6]     #LgMBH
kocb[:,4]= np.sqrt((np.log10(Kocb_loc[:,3]) - np.log10(Kocb_loc[:,4]))**2 +  (np.log10(Kocb_loc[:,5]) - np.log10(Kocb_loc[:,3]))**2) # error
#plt.errorbar(kocb[:,1],kocb[:,3], xerr=kocb[:,2] ,
#             yerr=kocb[:,4],fmt='.',color='blue',markersize=15, label='Kormendy 2013 classical bulge')

########input 30 local by Haring 04 ############
f2 ='data/Kormendy_elliptical.txt'
Kore_loc = np.loadtxt(f2)
#Kore_loc = Kore_loc[Kore_loc[:,-1]!=2]
kor_e = np.zeros([len(Kore_loc),5])
kor_e[:,0] = Kore_loc[:,0]
kor_e[:,1]= Kore_loc[:,1]    #LgMstar
kor_e[:,2]= Kore_loc[:,2]  # error
kor_e[:,3]= np.log10(Kore_loc[:,3]) + Kore_loc[:,6]     #LgMBH
kor_e[:,4]= np.sqrt((np.log10(Kore_loc[:,3]) - np.log10(Kore_loc[:,4]))**2 +  (np.log10(Kore_loc[:,5]) - np.log10(Kore_loc[:,3]))**2) # error
kor_e[:,4][kor_e[:,4]==np.inf] =  (np.log10(Kore_loc[:,5]) - np.log10(Kore_loc[:,3]))[kor_e[:,4]==np.inf]
#plt.errorbar(kor_e[:,1],kor_e[:,3], xerr=kor_e[:,2] ,yerr=kor_e[:,4],fmt='.',color='red',markersize=15, label='Kormendy 2013 elliptical')

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
Mpc = np.concatenate((kor_e[:,0], kocb[:,0]),axis=0)
solve_z = np.vectorize(solve_z)  #Calculate the redshift using the distance
kocb[:,0] =  solve_z(kocb[:,0])
kor_e[:,0] =  solve_z(kor_e[:,0])
# +  np.log10(Haring04[:,4] * 10 ** Haring04[:,5]))/2 #Sigma LgMBH
#plt.errorbar(kocb[:,1],kocb[:,3], xerr=kocb[:,2] ,yerr=kocb[:,4],fmt='.',color='black',markersize=15)
#plt.plot(kocb[:,1],kocb[:,3], '.',color='black',markersize=15)
#Hkc=mlines.Line2D([], [], color='black', ls='', marker='.', markersize=15)

#############################################################
###################fitting with MCMC#########################
x=np.concatenate((kor_e[:,1], kocb[:,1]),axis=0)
y=np.concatenate((kor_e[:,3], kocb[:,3]),axis=0)
yerr=(np.concatenate((kor_e[:,2], kocb[:,2]),axis=0)**2+np.concatenate((kor_e[:,4], kocb[:,4]),axis=0)**2)**0.5  # 0.2 is the uncertainty level for the L_R

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
    plt.plot(np.log10(kocb[:,0]+1),
             10**(kocb[:,3]-kocb[:,1]), '.',color='black',markersize=10)
    
    plt.plot(np.log10(kor_e[:,0]+1),
                 10**(kor_e[:,3]-kor_e[:,1]),'.',color='gray',markersize=10)
elif style ==1:
    xl = np.linspace(-0.9, 13, 100)
    plt.errorbar(np.log10(kocb[:,0]+1),
                 kocb[:,3]-(m_ml*kocb[:,1]+b_ml),yerr=(kocb[:,2]**2 + kocb[:,4]**2)**0.5 ,fmt='.',color='black',markersize=10)
    plt.errorbar(np.log10(kor_e[:,0]+1),
                 kor_e[:,3]-(m_ml*kor_e[:,1]+b_ml),yerr=(kor_e[:,2]**2 + kor_e[:,4]**2)**0.5 ,fmt='.',color='gray',markersize=10)
    ty=xl*0
    ty1=xl*0+np.std(y-(m_ml*x+b_ml))
    ty2=xl*0-np.std(y-(m_ml*x+b_ml))
    plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)

Kor_ell=mlines.Line2D([], [], color='gray', ls='', marker='.', markersize=15)
Kocb=mlines.Line2D([], [], color='black', ls='', marker='.', markersize=15)
######################
##%%
##calcualte the mean offset for the local:
#y_local=np.append(kocb[:,3]-(m_ml*kocb[:,1]+b_ml),kor_e[:,3]-(m_ml*kor_e[:,1]+b_ml))
#y_local_err = np.append((kocb[:,2]**2 + kocb[:,4]**2)**0.5, (kor_e[:,2]**2 + kor_e[:,4]**2)**0.5)
#weighted_offset = np.sum(np.asarray(y_local)*y_local_err) / np.sum(y_local_err)                              
#rms_offset = np.sqrt(np.sum((np.asarray(y_local)-weighted_offset)**2*y_local_err) / np.sum(y_local_err))
#print weighted_offset, rms_offset