import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mat
#mat.rcParams['font.family'] = 'STIXGeneral'
#plt.figure(figsize=(14.5,12))
############### with evolution-corrected or not? #################
host = 1

select= 0 #input('with evolution-corrected or not??? 0 = no;   1= yes, with dmag using dmag.py:')
if select == 0:
   dm=0
if select == 1:
   dm=1

import sys
sys.path.insert(0,'../py_tools')
from dmag import pass_dmag
########input Park local data ############
f1 ='../M_BH_relation/data/parklocal'
Pklc = np.loadtxt(f1)
ploc=np.zeros([len(Pklc),4])

#host= input('with bulge or total relation??? 0 = bulge;   1= total:')
if host == 0:
   ploc[:,1]= Pklc[:,6]    #LgV_sph
if host == 1:
   ploc[:,1]= Pklc[:,5]    #LgV_total
ploc[:,0]= Pklc[:,0]    #redshift
ploc[:,1]=0.4*(4.61+0.46-4.83)+ploc[:,1]-0.4*dm*pass_dmag(ploc[:,0])  #Change to LgR_sph, -0.4 times fainter
ploc[:,2]= Pklc[:,8]-0.03    #LogBHmassc, 0.03 is the updated recipe according to Daeseong's Email.
ploc[:,3]= Pklc[:,9]    #delta mass
#host_LR = np.log10(10 ** (0.4*(4.61-host_Mags)))
ploc[:,1] = 4.61 - ploc[:,1]/0.4
#Pkc=plt.errorbar(ploc[:,1],ploc[:,2],yerr=ploc[:,3],fmt='*',color='green',markersize=18)
plt.errorbar(ploc[:,1],ploc[:,2], xerr=0.2/0.4, yerr=ploc[:,3],fmt='.',color='gray',ecolor='gray',markersize=18)
Pkc=mlines.Line2D([], [], color='gray', ls='', marker='.', markersize=18)
loc=np.zeros([len(ploc),3])
loc[:,0]=ploc[:,1]
loc[:,1]=ploc[:,2]
loc[:,2]=ploc[:,3]
#print np.around(loc[:,0],2)

#############################################################
###################fitting with MCMC#########################
x=loc[:,0]
y=loc[:,1]
yerr=(loc[:,2]**2+(0.2/0.4)**2)**0.5  # 0.2 is the uncertainty level for the L_R
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
result = op.minimize(nll, [0.899, -1.509, 0.2], args=(x, y, yerr))
m_ml, b_ml,sint_ml= result["x"]
#print m_ml, b_ml, sint_ml, "ka=",lnlike(theta=[m_ml, b_ml, sint_ml],x=loc[:,0], y=loc[:,1], yerr=loc[:,2])

#m_ml, b_ml = 1.11, -3.7  # Park's first use

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
xl = np.linspace(-26, -19, 100)
m, b, sint =np.percentile(samples, 50,axis=0)
plt.plot(xl, m*xl+b, color="k", linewidth=4.0,zorder=-0.5)

def find_n(array,value):           #get the corresponding b for a given m 
    idx= (np.abs(array-value)).argmin()
    return array[idx]
m=np.percentile(samples,50,axis=0)[0]
#print samples[:,1][samples[:,0]==find_n(samples[:,0],m)]
rec=np.zeros([100,2])
for i in range(100):
    posi=np.random.uniform(16,84)
    m=np.percentile(samples,posi,axis=0)[0]
    b=samples[:,1][samples[:,0]==find_n(samples[:,0],m)][0]   #may find out many numbers
    rec[i,0],rec[i,1]=m,b
    plt.plot(xl, m*xl+b, color="lightgray", alpha=0.2,linewidth=7.0,zorder=-1000)
sm,sb=(np.amax(rec[:,0])-np.amin(rec[:,0]))/2,(np.amax(rec[:,1])-np.amin(rec[:,1]))/2
#plt.text(-26, 6.5, r"log(M$_{\rm BH}/$10$^{7}$M$_{\odot})$="+"{0:.2f}+{1:.2f}".format(round(b_ml+m_ml*10-7,2),round(m_ml,2))+r"log(L$_{\rm R}/$10$^{10}$L$_{\odot}$)"
#         ,color='blue',fontsize=25)
'''
if host == 0:
  plt.xlabel("$log(L_{R,bulge}/L_{\odot})$",fontsize=35)
if host == 1:
  plt.xlabel("$log(L_{R,total}/L_{\odot})$",fontsize=35)
plt.ylabel("$log(M_{BH}/M_{\odot})$",fontsize=35)
plt.xticks(np.arange(5,16.5,0.5))
plt.yticks(np.arange(6,16.5,0.5))
plt.axis([8.5,12.5,6.25,11])
plt.grid()
plt.tick_params(labelsize=25)
plt.show()
#for m, b, sint in samples[np.random.randint(len(samples), size=100)]:
#    plt.plot(xl, m*xl+b, color="k", alpha=0.2)
'''
'''
import corner
fig = corner.corner(samples, labels=["$m$", "$b$", "$sint$"],
                       quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})
'''
######################
