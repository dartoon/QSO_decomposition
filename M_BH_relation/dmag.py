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

def d_kcorr_R(z, filt = 'F140w'):
    from scipy import interpolate
    if filt == 'F140w':
        k_grid = np.array([[1, 0.849], [1.05, 0.842], [1.1, 0.834], [1.15, 0.814], [1.2, 0.800], [1.25, 0.790], [1.3, 0.770], [1.35, 0.744], [1.4, 0.720], [1.45, 0.697], [1.5, 0.673], [1.55, 0.649], [1.6, 0.625], [1.65, 0.594], [1.7, 0.560], [1.75, 0.520], [1.8, 0.472], [1.85, 0.425], [1.9, 0.379], [1.95, 0.335], [2, 0.288]])
        k_grid = k_grid.T
    if filt == 'F125w':
        k_grid = np.array([[1, 0.849], [1.05, 0.842], [1.1, 0.834], [1.15, 0.814], [1.2, 0.800], [1.25, 0.790], [1.3, 0.770], [1.35, 0.744], [1.4, 0.720], [1.45, 0.697], [1.5, 0.673], [1.55, 0.649], [1.6, 0.625], [1.65, 0.594], [1.7, 0.560], [1.75, 0.520], [1.8, 0.472], [1.85, 0.425], [1.9, 0.379], [1.95, 0.335], [2, 0.288]])
        k_grid = k_grid.T
    x = k_grid[0]
    y = k_grid[1]
    f = interpolate.interp1d(x, y)
    dm = f(z)
    return dm