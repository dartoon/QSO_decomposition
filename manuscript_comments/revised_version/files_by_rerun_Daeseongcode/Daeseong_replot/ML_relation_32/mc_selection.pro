
;Thu May  1 20:39:22 PDT 2014
;	
;	re-implementation of Monte Carlo simulation code for taking into account selection effects in M-L relation evolution signal
;	based on the SM code from Vardha Bennert (originally written by Tommaso Treu)

;Mon May 19 20:51:58 PDT 2014
;	starting from BHMF (not GLF) to bypass the active fraction effect

;Fri Jan 13 14:22:14 KST 2017
;	running simulation for the M-L project by Xuheng


function fct_BHMF, logMBH
	;Schulze & Wisotzki (2010) active local BHMF with broad-line AGNs
	phi_char = 2.75D-6 ;Mpc^-3
	M_char = 10D^8.11 ;Msun
	alpha = -2.11D
	beta = 0.50D

	MBH = 10D^logMBH

	;BHMF modified Schechter function fit in d log M
	Phi_M = MBH * alog(10D) * (phi_char/M_char) * (MBH/M_char)^alpha * exp( -(MBH/M_char)^beta)

	return,Phi_M
end


function gen_sample, py, N,sig_int,dellogM,lowlogM,upplogM
common baseline, local

	;measurement errors on X and Y
	elogM = 0.4D
	elogL = 0.2D

	;baseline local relation with redshift evolution effect (dellogM = gamma*log(1+z))
	;log M = alpha + beta*(log L - 10) + gamma*log(1+z)
	;py = local.alpha + local.beta * px + dellogM
	;[log M - alpha - gamma*log(1+z)] / beta = (log L - 10)
	px = (py - local.alpha - dellogM) / local.beta

	;add Gaussian intrinsic scatter on y
	gnoise = randomN(seed,N,/double)
	py += (gnoise - mean(gnoise,/DOUBLE))/stddev(gnoise,/DOUBLE) * sig_int

	;add Gaussian measurement error on x and y (uncorrelated)
	gnoise = randomN(seed,N,/double)
	px += (gnoise - mean(gnoise,/DOUBLE))/stddev(gnoise,/DOUBLE) * elogL 
	gnoise = randomN(seed,N,/double)
	py += (gnoise - mean(gnoise,/DOUBLE))/stddev(gnoise,/DOUBLE) * elogM 

	;add selection effects by hard thresholds in logMBH
	idx_sel = where(py ge lowlogM and py le upplogM, c_sel)
		if (c_sel lt 1) then stop,'No sample at the criteria'

	return, {X:px, Y:py, NS:c_sel, SX:px[idx_sel], SY:py[idx_sel]}
end


function make_sample, py, N,sig_int
common baseline, local

	;measurement errors on X and Y
	elogM = 0.4D
	elogL = 0.2D

	;baseline local relation with redshift evolution effect (dellogM = gamma*log(1+z))
	;[log M - alpha - gamma*log(1+z)] / beta = (log L - 10)
	px = (py - local.alpha) / local.beta


	;add Gaussian intrinsic scatter on y
	gnoise = randomN(seed,N,/double)
	py += (gnoise - mean(gnoise,/DOUBLE))/stddev(gnoise,/DOUBLE) * sig_int

	;add Gaussian measurement error on x and y (uncorrelated)
	gnoise = randomN(seed,N,/double)
	px += (gnoise - mean(gnoise,/DOUBLE))/stddev(gnoise,/DOUBLE) * elogL 
	gnoise = randomN(seed,N,/double)
	py += (gnoise - mean(gnoise,/DOUBLE))/stddev(gnoise,/DOUBLE) * elogM 

	return, {X:px, Y:py}
end


PRO mc_selection
;--------------------  extra jobs  -------------------------------{{{
journal,'log_mc_sel.pro' 
;profiler,/reset & profiler,/system & profiler                    ;
start_mem = memory(/current)   ;to know how much memory required ;
itime = systime(/sec)          ;to compute elapsed time          ;
spawn,'date'
;-----------------------------------------------------------------}}}

CPU,TPOOL_NTHREADS = 4

common baseline, local
;Bennert+10 local RM AGN M-L relation (logM = alpha + beta log L/10^10)
local = { $
	alpha   :7.40D, $
	beta    :1.11D, $
	int_scat:0.27D}

;---------------------------------------------------------
;*make appropriate random sample using nonparametric method (transformation method)
;---------------------------------------------------------
;logMBH = arrgen(6.7D,9.5D,0.00001D) ; => logMBH = [5x10^6 - 3x10^9]
logMBH = arrgen(6.5D,11D,0.00001D) ; => logMBH = [5x10^6 - 3x10^9]
help,logMBH 
;stop

Mprior = fct_BHMF(logMBH)


;---------------------------------------------------------
;-checkprior
;---------------------------------------------------------
;plot,logMBH, alog10(Mprior)
;null = get_kbrd()


;---------------------------------------------------------
;-calculate integral function of the prior
;---------------------------------------------------------
;normalized cumulative Mprior
cMprior = total(Mprior,/cumulative)/total(Mprior)   ;cumulative probability distribution

;plot, cMprior, psym=1
;null = get_kbrd()

;---------------------------------------------------------
;-draw a set of random points from the Mprior
;---------------------------------------------------------
N_real = 2D+8
Xran = randomu(seed, N_real,/double) 
Xran = Xran[sort(Xran)]
;plothist, Xran, bin=0.01
;null = get_kbrd()

;-linear interpolation 
pMBH = INTERPOL(logMBH,cMprior,Xran,/SPLINE) ; rank-based mapping using inverse cumulative distribution
;plot,pMBH,psym=1
;stop
;plot, cMprior,pMBH ,/xsty,/ysty
;plot, cMprior,Xran ,/xsty,/ysty
;plot, cMprior,logMBH 
;null = get_kbrd()
;plot, logMBH, pMBH,/xsty,/ysty
;plot, Xran, pMBH
;null = get_kbrd()
;stop
undefine, Xran
undefine, Mprior
undefine, cMprior


help,pMBH
PLOTHIST, pMBH,/AUTOBIN,tit='!6',ytit='!6N',xtit=textoidl('!6log !8M!6_{BH}'),charsize=2,/ylog,xsty=3,ysty=3
;null = get_kbrd()
;stop


;---------------------------------------------------------
;-produce Monte Calro sample
;---------------------------------------------------------
py = temporary(pMBH)
N = N_real

;sig_int = 0.4D
;dellogM = 0D
;lowlogM = 7.5D
;upplogM = 9D
;
;help,py,N,sig_int,dellogM,lowlogM,upplogM
;
;sample = gen_sample(py, N,sig_int,dellogM,lowlogM,upplogM)
;help,sample,/str
;snpx = sample.SX
;snpy = sample.SY
;npx = sample.X
;npy = sample.Y
;plot,npx+10D,npy,psym=1,xsty=1,ysty=1,xr=[6,13],yr=[3,11]
;oplot,snpx+10D,snpy,psym=1,color=cgcolor('red')
;oplot,(py - local.alpha - dellogM) / local.beta +10D,py,color=cgcolor('green')
;
;stop


;---------------------------------------------------------
;-evaluate likelihood function in 2D grids (gamma, sig_int)
;---------------------------------------------------------
;-read data
READCOL,'../../data_MBH-LRhost_NewSampleOnly.txt',name,z,logMBH,elogMBH,logL,elogL,flag,format='A,D,D,D,D,D,L',count=Ndata,comment='#'

;-selection threshold in log M for each sample
lowlogM = dblarr(Ndata)
upplogM = dblarr(Ndata)

i_RM = where(flag eq 0, c_RM)
if c_RM gt 0 then begin
	lowlogM[i_RM] = 7D 
	upplogM[i_RM] = 9D
endif

i_SS = where(flag eq 2, c_SS)
if c_SS gt 0 then begin
	lowlogM[i_SS] = 7.3D 
	upplogM[i_SS] = 8.2D
endif

i_SW = where(flag eq 1 or flag eq 3, c_SW)
if c_SW gt 0 then begin
	lowlogM[i_SW] = 7.7D 
	upplogM[i_SW] = 9.1D
endif


i_B11 = where(flag eq 5, c_B11)
if c_B11 gt 0 then begin
	lowlogM[i_B11] = 7.8D 
	upplogM[i_B11] = 9.3D
endif

i_SS13 = where(flag eq 6, c_SS13)
if c_SS13 gt 0 then begin
	lowlogM[i_SS13] = 7.1D 
	upplogM[i_SS13] = 9.3D
endif


;new HST QSOs from Xuheng (32)
i_newHSTQSO_Xuheng = where(flag eq 7, c_Xu)
if c_Xu gt 0 then begin
	lowlogM[i_newHSTQSO_Xuheng] = 7.7D 
	upplogM[i_newHSTQSO_Xuheng] = 8.8D
endif




XX = logL - 10D
YY = logMBH
EX = elogL
EY = elogMBH


;---------------------------------------------------------
Nsample = N_real                             ;how many?
hbinsize = 0.1D;0.04D;0.2D;                  ;how small?
hbinsize = 0.05D

;-slope range (0, 5) & stepsize 0.01
;input_gam = arrgen(-2.0D, 3.0D, 0.1D)         ;how fine?
;input_gam = arrgen(-1.0D, 1.0D, 0.05D)         ;how fine?
input_gam = arrgen(-1.0D, 6.0D, 0.1D)
;-sig_int range (0, 1) & stepsize 0.01       ;how wide?
input_sig = arrgen(0.02D, 0.98D, 0.02D)
;---------------------------------------------------------

;-2D array for saving 2lnP=-chi2 in (gam,sig) grids
grid_2lnP = dblarr(n_elements(input_gam),n_elements(input_sig))
help, grid_2lnP 


;;-read AF and interpolate
;READCOL,'../dutycycle/S09_BHdutycycle.txt',AF_z,AF_logMBH,AF_P0,format='D,D,D',/sil,comment='#'
;P0GridInterp = MIN_CURVE_SURF(AF_P0, AF_z, AF_logMBH, GS=[0.01, 0.01],xpout=xout,ypout=yout, /DOUBLE)


it = systime(/sec)
for i=0L,n_elements(input_gam)-1L do begin
   if i eq 1 then print,'-elapsed time for one setup = ',strtrim((systime(/sec)-it)/60D,2), ' min'
   print, string(i,format='(I0)')+' of '+string(n_elements(input_gam)-1L,format='(I0)')+': current input gamma = '+string(input_gam[i],format='(D5.2)')

	for j=0L,n_elements(input_sig)-1L do begin

		sample = make_sample(py, Nsample,input_sig[j])
		npx = temporary(sample.X)
		npy = temporary(sample.Y)

		for k=1L,n_elements(XX)-1L do begin

			dellogM = input_gam[i]*alog10(1D + z[k])
			znpy = npy + dellogM

			;add selection effects by hard thresholds in logMBH
			idx_sel = where(znpy ge lowlogM[k] and znpy le upplogM[k], c_sel)
				if (c_sel lt 1) then stop,'No sample at the criteria'
			snpx = npx[idx_sel]
			snpy = znpy[idx_sel]

			;select the sample in [logL-elogL, logL+elogL] 
			i_slogL = where(snpx ge XX[k]-EX[k] and snpx le XX[k]+EX[k], c_slogL)
				if c_slogL lt 1 then begin
					grid_2lnP[i,j] += 2D*alog(1D-20)
					continue
				endif
			hist = histogram(snpy[i_slogL], BINSIZE=hbinsize, LOCATIONS=yh, MIN=lowlogM[k], MAX=upplogM[k])
			nhist = hist/(total(hist)*hbinsize)
			;plot,yh,nhist,psym=10

			prob = INTERPOL(nhist,yh,YY[k],/SPLINE) ; could be imporved? cspline? hermite?

;         nhistA = nhist*0D
;         norm = 0D
;         for iii=0L,n_elements(yh)-1L do begin
;            dist = sqrt((xout - z[k])^2 + (yout - yh[iii])^2)
;            min = min(dist, location)
;            ind = ARRAY_INDICES(dist, location)
;            nhistA[iii] = nhist[iii] * P0GridInterp[ind[0],ind[1]]
;            norm += nhist[iii] * P0GridInterp[ind[0],ind[1]]
;         endfor
;
;         nhistA /= (norm*hbinsize)
;        ;oplot,yh,nhistA,psym=10,color=cgCOLOR('red')
;
;         prob = INTERPOL(nhistA,yh,YY[k],/SPLINE)


			if prob le 1D-20 $
				then grid_2lnP[i,j] += 2D*alog(1D-20) $
				else grid_2lnP[i,j] += 2D*alog(prob) ; 2linP = - chi2 => should be maximized

		endfor
	endfor
endfor


SAVE,grid_2lnP,input_gam,input_sig,filename='mcrun_result.sav'
SPAWN,'date'


print
print,'*********************************************************'
print,'using flat prior for int. scat.'
print,'*********************************************************'


;-find the maximum likelihood point in 2D grid
max_2lnP = max(grid_2lnP, location,/NaN)
i_max = array_indices(grid_2lnP, location)
best_gam = input_gam[i_max[0]]
best_sig = input_sig[i_max[1]]
print, 'max grid index   = ',i_max
print, 'max 2lnP         = ',max_2lnP
print, 'gamma at max P   = ',best_gam 
print, 'sig_int at max P = ',best_sig 


;plot
thick = 4
charsize = 2.5
charthick = 1
symsize = 2


window,2,xs=1000,ys=700
CONTOUR, max_2lnP - grid_2lnP, input_gam,input_sig, levels=[0,2.30,6.17,11.80],$
	thick=thick,charsize=charsize,charthick=charthick,$
	/fill,C_COLORS=[cgcolor('orange'),cgcolor('lime green'),cgcolor('royal blue'),cgcolor('black')] ,$
	xsty=3,ysty=3,$
	tit='!6',$
	xtit=textoidl('\gamma'),$
	ytit=textoidl('\sigma_{int}')

CONTOUR, max_2lnP - grid_2lnP, input_gam,input_sig, levels=[0,2.30,6.17,11.80],/over,$
	C_LABELS=[1,1,1,1],C_ANNOTATION=['best','68.3%','95.4%','99.7%'],c_charsize=charsize-1,c_charthick=charthick
oPLOT,[best_gam],[best_sig],psym=7,symsize=symsize+2,thick=thick+1,color=cgcolor('magenta')

WRITE_PNG,'sel_eff_test.png',tvrd(true=1)


;---------------------------------------------------------
;-marginalizes likelihood function to get 1D prob on each paremeter and determine confidence interval
;---------------------------------------------------------
;-1) gamma
marginal_Pgam = dblarr(n_elements(input_gam))
for i=0L,n_elements(input_gam)-1L do marginal_Pgam[i] = total(exp(0.5D*(grid_2lnP[i,*] - max_2lnP)))  ; = sum(p/pmax)

	maxPgam = max(marginal_Pgam,i_maxPgam)
	i_sorted = sort(marginal_Pgam)
	marPgam_sorted = marginal_Pgam[i_sorted]
	input_gam_sorted = input_gam[i_sorted]

	int_i = arrgen(maxPgam,0D,maxPgam/100D)
	int_i[n_elements(int_i)-1] = 1D-10
	;help,int_i

	cint_p = dblarr(n_elements(int_i))

	;for i=0L,n_elements(int_i)-1L do begin
	;	i_cum = where(marPgam_sorted ge int_i[i], c_cum)
	;		if c_cum lt 1 then stop,'no elements in cum. sum.'
	;	cint_p[i] = total(marPgam_sorted[i_cum])
	;endfor
	for i=0L,n_elements(int_i)-1L do cint_p[i] = total(marPgam_sorted[where(marPgam_sorted ge int_i[i], c_cum)])
	cint_p /= total(marPgam_sorted)

	CL_per = [0.683, 0.954, 0.997] ;1sigma 2sigma 3sigma

	;pCL_per = INTERPOL(int_i,cint_p,CL_per)
	LINTERP, cint_p, int_i, CL_per, pCL_per
	print,pCL_per

	;for 1sigma error
	gam_1sig = input_gam_sorted[where(marPgam_sorted ge pCL_per[0])]
	print,'best gamma = ',best_gam
	print,'gam error  = ',(max(gam_1sig) - min(gam_1sig))/2D
	print,'gamma-1sig = ',min(gam_1sig)
	print,'gamma+1sig = ',max(gam_1sig)


;-2) intrinsic scatter
marginal_Psig = dblarr(n_elements(input_sig))
for j=0L,n_elements(input_sig)-1L do marginal_Psig[j] = total(exp(0.5D*(grid_2lnP[*,j] - max_2lnP)))  ; = sum(p/pmax)

	maxPsig = max(marginal_Psig,i_maxPsig)
	i_sorted = sort(marginal_Psig)
	marPsig_sorted = marginal_Psig[i_sorted]
	input_sig_sorted = input_sig[i_sorted]

	int_i = arrgen(maxPsig,0D,maxPsig/100D)
	int_i[n_elements(int_i)-1] = 1D-10
	;help,int_i

	cint_p = dblarr(n_elements(int_i))

	;for i=0L,n_elements(int_i)-1L do begin
	;	i_cum = where(marPsig_sorted ge int_i[i], c_cum)
	;		if c_cum lt 1 then stop,'no elements in cum. sum.'
	;	cint_p[i] = total(marPsig_sorted[i_cum])
	;endfor
	for i=0L,n_elements(int_i)-1L do cint_p[i] = total(marPsig_sorted[where(marPsig_sorted ge int_i[i], c_cum)])
	cint_p /= total(marPsig_sorted)

	CL_per = [0.683, 0.954, 0.997] ;1sigma 2sigma 3sigma

	;pCL_per = INTERPOL(int_i,cint_p,CL_per)
	LINTERP, cint_p, int_i, CL_per, pCL_per
	print,pCL_per

	;for 1sigma error
	sig_1sig = input_sig_sorted[where(marPsig_sorted ge pCL_per[0])]
	print,'best sig = ',best_sig
	print,'sig error  = ',(max(sig_1sig) - min(sig_1sig))/2D
	print,'sig-1sig = ',min(sig_1sig)
	print,'sig+1sig = ',max(sig_1sig)



;null = get_kbrd()
print
print,'*********************************************************'
print,'using informative prior for int. scat.'
print,'*********************************************************'


;---------------------------------------------------------
;-compute posterior using informative prior for int. scat. (above priors for gamma & int. scat. are assumed to be uniform)
;---------------------------------------------------------
;1) Gultekin+09 M-L relation sig_int = 0.38+/-0.09
;sig_int = 0.38D
;esig_int = 0.09D
;2) McConnell+13 M-L relation sig_int = 0.52+/-0.06
;sig_int = 0.52D
;esig_int = 0.06D
;3) Bennert+10 local RM AGNs M-L relation sig_int = 0.21+/-0.08                                                                                        
;-rederived for R-band wiht linmix_err
sig_int = 0.27D                                                                                                                                        
esig_int = 0.10D                                                                                                                                       


grid_2lnP_IP = grid_2lnP
for j=0L,n_elements(input_sig)-1L do grid_2lnP_IP[*,j] -= ((input_sig[j]-sig_int)/esig_int)^2  ;-chi2 = 2lnP

;-find the maximum likelihood point in 2D grid
max_2lnP_IP = max(grid_2lnP_IP, location,/NaN)
i_max = array_indices(grid_2lnP_IP, location)
best_gam = input_gam[i_max[0]]
best_sig = input_sig[i_max[1]]
print, 'max grid index   = ',i_max
print, 'max 2lnP_IP         = ',max_2lnP_IP
print, 'gamma at max P   = ',best_gam 
print, 'sig_int at max P = ',best_sig 


window,3,xs=1000,ys=700
CONTOUR, max_2lnP_IP - grid_2lnP_IP, input_gam,input_sig, levels=[0,2.30,6.17,11.80],$
	thick=thick,charsize=charsize,charthick=charthick,$
	/fill,C_COLORS=[cgcolor('orange'),cgcolor('lime green'),cgcolor('royal blue'),cgcolor('black')] ,$
	xsty=3,ysty=3,$
	tit='!6',$
	xtit=textoidl('\gamma'),$
	ytit=textoidl('\sigma_{int}')

CONTOUR, max_2lnP_IP - grid_2lnP_IP, input_gam,input_sig, levels=[0,2.30,6.17,11.80],/over,$
	C_LABELS=[1,1,1,1],C_ANNOTATION=['best','68.3%','95.4%','99.7%'],c_charsize=charsize-1,c_charthick=charthick
oPLOT,[best_gam],[best_sig],psym=7,symsize=symsize+2,thick=thick+1,color=cgcolor('magenta')

WRITE_PNG,'sel_eff_test_IP.png',tvrd(true=1)


;---------------------------------------------------------
;-marginalizes likelihood function to get 1D prob on each paremeter and determine confidence interval
;---------------------------------------------------------
;-1) gamma
marginal_Pgam = dblarr(n_elements(input_gam))
for i=0L,n_elements(input_gam)-1L do marginal_Pgam[i] = total(exp(0.5D*(grid_2lnP_IP[i,*] - max_2lnP_IP)))  ; = sum(p/pmax)

	maxPgam = max(marginal_Pgam,i_maxPgam)
	i_sorted = sort(marginal_Pgam)
	marPgam_sorted = marginal_Pgam[i_sorted]
	input_gam_sorted = input_gam[i_sorted]

	int_i = arrgen(maxPgam,0D,maxPgam/100D)
	int_i[n_elements(int_i)-1] = 1D-10
	;help,int_i

	cint_p = dblarr(n_elements(int_i))

	;for i=0L,n_elements(int_i)-1L do begin
	;	i_cum = where(marPgam_sorted ge int_i[i], c_cum)
	;		if c_cum lt 1 then stop,'no elements in cum. sum.'
	;	cint_p[i] = total(marPgam_sorted[i_cum])
	;endfor
	for i=0L,n_elements(int_i)-1L do cint_p[i] = total(marPgam_sorted[where(marPgam_sorted ge int_i[i], c_cum)])
	cint_p /= total(marPgam_sorted)

	CL_per = [0.683, 0.954, 0.997] ;1sigma 2sigma 3sigma

	;pCL_per = INTERPOL(int_i,cint_p,CL_per)
	LINTERP, cint_p, int_i, CL_per, pCL_per
	print,pCL_per

	;for 1sigma error
	gam_1sig = input_gam_sorted[where(marPgam_sorted ge pCL_per[0])]
	print,'best gamma = ',best_gam
	print,'gam error  = ',(max(gam_1sig) - min(gam_1sig))/2D
	print,'gamma-1sig = ',min(gam_1sig)
	print,'gamma+1sig = ',max(gam_1sig)


;-2) intrinsic scatter
marginal_Psig = dblarr(n_elements(input_sig))
for j=0L,n_elements(input_sig)-1L do marginal_Psig[j] = total(exp(0.5D*(grid_2lnP_IP[*,j] - max_2lnP_IP)))  ; = sum(p/pmax)

	maxPsig = max(marginal_Psig,i_maxPsig)
	i_sorted = sort(marginal_Psig)
	marPsig_sorted = marginal_Psig[i_sorted]
	input_sig_sorted = input_sig[i_sorted]

	int_i = arrgen(maxPsig,0D,maxPsig/100D)
	int_i[n_elements(int_i)-1] = 1D-10
	;help,int_i

	cint_p = dblarr(n_elements(int_i))

	;for i=0L,n_elements(int_i)-1L do begin
	;	i_cum = where(marPsig_sorted ge int_i[i], c_cum)
	;		if c_cum lt 1 then stop,'no elements in cum. sum.'
	;	cint_p[i] = total(marPsig_sorted[i_cum])
	;endfor
	for i=0L,n_elements(int_i)-1L do cint_p[i] = total(marPsig_sorted[where(marPsig_sorted ge int_i[i], c_cum)])
	cint_p /= total(marPsig_sorted)

	CL_per = [0.683, 0.954, 0.997] ;1sigma 2sigma 3sigma

	;pCL_per = INTERPOL(int_i,cint_p,CL_per)
	LINTERP, cint_p, int_i, CL_per, pCL_per
	print,pCL_per

	;for 1sigma error
	sig_1sig = input_sig_sorted[where(marPsig_sorted ge pCL_per[0])]
	print,'best sig = ',best_sig
	print,'sig error  = ',(max(sig_1sig) - min(sig_1sig))/2D
	print,'sig-1sig = ',min(sig_1sig)
	print,'sig+1sig = ',max(sig_1sig)



;--------------------  extra jobs  -------------------------------{{{
print
etime = (systime(/sec)-itime)/60D
if etime ge 1D then print,'*Elapsed time= ',strtrim(etime,2), ' min (',strtrim(etime/60D,2),' hrs)' $
   else print,'*Elapsed time= ',strtrim(etime*60D,2), ' sec'
print,'*Memory used= ',strtrim((memory(/highwater)-start_mem)/2D^20D,2),' MB'
;profiler,/report,filename='profiler_report.txt' & profiler,/reset
;print
spawn,'date'
journal
;-----------------------------------------------------------------}}}

;STOP,' *** END of CODE *** '

END
