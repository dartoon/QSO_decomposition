# to derive the absolute magnitue of RXJ1131
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
zs=0.65
H=70             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
dl=(1+zs)*c*EE(zs)/H   #in Mpc
dl=dl*10**6
mr_bulge=21.8
Mr_bulge=mr_bulge-0.7-5*(np.log10(dl)-1)
mr_disk=20.067
Mr_disk=mr_disk-0.3-5*(np.log10(dl)-1)

#print "M_bulge=",Mr_bulge,"M_disk=",Mr_disk
