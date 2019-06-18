#paper: https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.5410D/graphics
#file from: https://lerma.obspm.fr/huertas/form_CANDELS?

col1 : ID from the CANDELS catalog
col2 : RA
col3 : DEC
col4 : FIELD
col5 : B_T_m, stellar mass bulge-to-total ratio, corrected as explained in the paper; fig.5. It is set to 0 (1) when the best profile is identified to be a pure disk (pure bulge)
col6-col9  : Profile probability 
P_B : probability for single sersic profile with n free 
P_D : probability for Pure Exponential profile
P_BD :probability for Bulge+disk profile
P_PD : Pseudo-bulge + disk profile

col10-col13 : Morphological probability (Huertas-Company+ 2015)
F_SPH   : probability that the galaxy is a spheroid
F_DISK  : probability that the galaxy is a disk
F_IRR    : probability that the galaxy is irregular
F_PS     : probability that the galaxy is a point spread function/star

col14 : setup, used for the correction of ‘ambiguous’ cases [see fig 5 in the paper]
col15 : n_comp, number of components 
col16 : profile -  1: Pure Sersic ; 2 Pure Exponential ; 3 Bulge+disk ; 4 Pseudo-bulge + disk ; 5 Unclassified ; 6 Star ;
col17 : n_profile, if n_profile=1 more profiles are possible according to the probability threshold. 
col18 : q_fit, quality of the fit done comparing the total magnitude in the H band with the one from CANDELS.

# Total
col19-col58 : Galfit result for the single sersic fit
col59 : lmass: log[mass/Msol]
col60 : l68_lmass : lower estimation
col61 : u68_lmass : upper estimation
col62 : uv rest frame color
col63 : vj rest frame color
col64 : chi2 from the SED fitting

# Bulge
col65-col104 : result from the Galfit fit
col105 : lmass: log[mass/Msol]
col106 : l68_lmass : lower estimation
col107 : u68_lmass : upper estimation
col108 : uv rest frame color
col109 : error uv
col110 : vj rest frame color
col111 : error vj
col112 : chi2 from the SED fitting

#Disk
col113-col152 : result from the Galfit fit
col153 : lmass: log[mass/Msol]
col154 : l68_lmass : lower estimation
col155 : u68_lmass : upper estimation
col156 : uv rest frame color
col157 : error uv
col158 : vj rest frame color
col159 : error uv
col160 : chi2 from the SED fitting

#Error
col161-col256 

# Bias 
col 257-col334


