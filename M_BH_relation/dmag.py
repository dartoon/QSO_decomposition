import numpy as np
import matplotlib.pyplot as plt
###### use Taka data ####

def dmag(z): 
     return 3.70*np.log10(1+z)
#test=np.array([1,2,3])
#print dmag(test)
#print dmag(np.linspace(0,5,50))/np.linspace(0,5,50)
#z1=1.69
#z2=2.69
#print "z=1.69",dmag(z1),y0+0.69*y1,"z=2.69",dmag(z2),y0+y1+0.69*y2
#plt.plot(np.linspace(0,5,50),dmag(np.linspace(0,5,50)),'k')
#plt.show()
