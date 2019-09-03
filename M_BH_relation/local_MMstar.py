import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mat
#mat.rcParams['font.family'] = 'STIXGeneral'
#plt.figure(figsize=(14.5,12))
############### with evolution-corrected or not? #################
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
plt.errorbar(bloc[:,1],bloc[:,3], xerr=bloc[:,2] ,yerr=bloc[:,4],fmt='.',color='gray',markersize=15)
#plt.plot(bloc[:,1],bloc[:,3], '.',color='gray',markersize=15)
Bkc=mlines.Line2D([], [], color='gray', ls='', marker='.', markersize=15)

########input 30 local by Haring 04 ############
f2 ='data/Haring04.txt'
Haring04 = np.loadtxt(f2)
hloc = np.zeros([len(Haring04),5])
hloc[:,1]= np.log10(Haring04[:,0] * 10 ** Haring04[:,1])    #LgMstar
hloc[:,2]= 0.18  # Mention in the Haring04 paper, note in fig 2
hloc[:,3]= np.log10(Haring04[:,2] * 10 ** Haring04[:,5])    #LgMBH
hloc[:,4]= (abs(np.log10(Haring04[:,2] + Haring04[:,3]) - np.log10(Haring04[:,2])) + abs(np.log10(Haring04[:,2] - Haring04[:,4]) - np.log10(Haring04[:,4])))/2
# +  np.log10(Haring04[:,4] * 10 ** Haring04[:,5]))/2 #Sigma LgMBH
plt.errorbar(hloc[:,1],hloc[:,3], xerr=hloc[:,2] ,yerr=hloc[:,4],fmt='.',color='black',markersize=15)
#plt.plot(hloc[:,1],hloc[:,3], '.',color='black',markersize=15)
Hkc=mlines.Line2D([], [], color='black', ls='', marker='.', markersize=15)

#############################################################
###################fitting all togethera with MCMC#########################
x=np.append(bloc[:,1], hloc[:,1])
y=np.append(bloc[:,3], hloc[:,3])
yerr=(np.append(bloc[:,2], hloc[:,2])**2+np.append(bloc[:,4], hloc[:,4])**2)**0.5  # 0.2 is the uncertainty level for the L_R
def lnlike(theta, x, y, yerr):
    m, b, sint= theta
    model = m * x + b
#    yerr=(m*(np.append(bloc[:,2], hloc[:,2]))**2+np.append(bloc[:,4], hloc[:,4])**2)**0.5  # right error level
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [0.93027905, -1.95536508, 0.35], args=(x, y, yerr))
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
sampler.run_mcmc(pos, 1000)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
xl = np.linspace(5, 13, 100)
m, b, sint =np.percentile(samples, 50,axis=0)
plt.plot(xl, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)

def find_n(array,value):           #get the corresponding b for a given m 
    idx= (np.abs(array-value)).argmin()
    return array[idx]
m=np.percentile(samples,50,axis=0)[0]
#print samples[:,1][samples[:,0]==find_n(samples[:,0],m)]
for i in range(100):
    posi=np.random.uniform(16,84)
    m=np.percentile(samples,posi,axis=0)[0]
    b=samples[:,1][samples[:,0]==find_n(samples[:,0],m)][0]   #may find out many numbers
    plt.plot(xl, m*xl+b, color="lightgray", alpha=0.2,linewidth=7.0,zorder=-1000)

plt.text(9.3, 6.24, "log$(M_{BH}/10^{7}M_{\odot})$=%s+%slog$(M_*/10^{10}M_{\odot})$"%(round(b_ml+m_ml*10-7,2),round(m_ml,2)),color='blue',fontsize=25)

#import corner
#fig = corner.corner(samples, labels=["$m$", "$b$", "$sint$"],
#                       quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})