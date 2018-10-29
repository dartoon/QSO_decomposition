import numpy as np
####For Mg########
#####Recipi######
#########Peng's ##############
#a=np.log10(6.1)+6 # =6.7853   
#b=0.47
#############use McGill APJ 2008 modifiy with Park's work#########
#log(Mbh/Mo)=(6.767-0.144) +/- 0.055 + 0.47log(L3000/10^44 erg/s) + 2 log (FWHM/1000 km/s) +/- 0.233
#####Calibrate function for Mg##################
f0 ='data/Mg_Peng'
Mg = np.loadtxt(f0)
Mg[:,4]=Mg[:,4]+0.21   #change Mg from Vega to AB
z_Mg=Mg[:,1]             #redshift
L_Mg=Mg[:,7]             #continuum luminiosities
F_A_Mg=Mg[:,6]           # FWHM in A (wavelength unit)
real_Mg=9+np.log10(Mg[:,9])         # data are in 10**9 unit
FWHM_Mg=F_A_Mg*2.99792*10**5/(1+z_Mg)/2798   # 2798 is MgII wavelength in restframe
def Mgcal(L,FWHM):
    a=6.767-0.144
    b=0.47
    logM=a + b*(L-44.) +2 * np.log10(FWHM/1000.)
    return logM
vec_Mgcal=np.vectorize(Mgcal)
err_Mg= (0.055**2 + 0.233 **2) **0.5
#print Mgcal(43.95,3600)

#########For Hb########
#########Recipi######
#########Peng's ##############
#a=np.log10(5.9)+6         ## =6.7853   
#b=0.69
#####Calibrate tool##################
#######for Hb#######
f1 ='data/Hb_Peng'
H = np.loadtxt(f1)
H[:,4]=H[:,4]+0.21   #change Mg from Vega to AB
z_hb=H[:,1]             #redshift
L_hb=H[:,7]             #continuum luminiosities
F_A_hb=H[:,6]           # FWHM in A (wavelength unit)
real_hb=9+np.log10(H[:,9])         # data are in 10**9 unit
FWHM_hb=F_A_hb*2.99792*10**5/(1+z_hb)/4861   # 4861 is HII wavelength in restframe
def Hcal(L,FWHM):
    a=7.026-0.144
    b=0.518
    logM=a + b*(L-44.) +2 * np.log10(FWHM/1000.)
    return logM
vec_Hcal=np.vectorize(Hcal)
err_hb= (0.056**2 + 0.235 **2) **0.5

#########For C#########
########Recipi#########
#########Peng's ##############
#####Calibrate tool##################
#######for C#######
f2 ='data/C_Peng'
C= np.loadtxt(f2)
C[:,4]=C[:,4]+0.21   #change Mg from Vega to AB
z_c=C[:,1]             #redshift
L_c=C[:,7]             #continuum luminiosities
F_A_c=C[:,6]           # FWHM in A (wavelength unit)
real_c=9+np.log10(C[:,9])         # data are in 10**9 unit
FWHM_c=F_A_c*2.99792*10**5/(1+z_c)/1549   # 1549 is C wavelength in restframe
def Ccal(L,FWHM):
    a=np.log10(4.5)+6      -0.331475893525   ## =6.3217   ###including the calibrate...
    b=0.53
    logM=a + b*(L-44.) +2 * np.log10(FWHM/1000.)
    return logM
vec_Ccal=np.vectorize(Ccal)


