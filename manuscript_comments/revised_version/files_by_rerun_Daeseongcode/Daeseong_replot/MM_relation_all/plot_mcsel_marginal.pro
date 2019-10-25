
;Sun Jul 21 23:55:06 KST 2013
;	plot 2D distribution of parameters for M-L relation from Monte Carlo simulation results for selection effects in M-L relation evolution signal


PRO plot_mcsel_marginal

RESTORE,'../mcrun_result.sav'

;--------------------  extra jobs  -------------------------------{{{
journal,'log_mc_sel.pro' 
;profiler,/reset & profiler,/system & profiler                    ;
start_mem = memory(/current)   ;to know how much memory required ;
itime = systime(/sec)          ;to compute elapsed time          ;
;-----------------------------------------------------------------}}}


print
print,'*********************************************************'
print,'using flat prior for int. scat.'
print,'*********************************************************'
help,grid_2lnP


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



	;---------------------------------------------------------
	;eps version

	thick = 4
	charsize = 1.5
	charthick = 3
	symsize = 2

	VSYM,24,thick=thick;,/FILL

	SET_PLOT,'PS'

	!p.multi = 0
	 !x.margin = [6,2.5]
	 !y.margin = [3.0,1.5]
	 ;!P.FONT = -1
	 !x.thick = 3
	 !y.thick = 3
	 !x.ticklen = 0.04
	 !y.ticklen = 0.02

	 !x.margin = [12,2]
	 !y.margin = [8,2]

	DEVICE,/PORTRAIT,/COLOR,BITS=8,/ISOLatin1,$
		xsize=25,ysize=15,/ENCAPSUL ;/Times/Helvetica/Palantino
	;DEVICE,/PORTRAIT,/COLOR,BITS=8,/ISOLatin1,$
	;		xsize=18,xoffset=1.5,ysize=20,yoffset=7.5 ;/Times/Helvetica/Palantino

	DEVICE,FileName='MC_seleff_flatprior.eps'

;	MULTIPLOT,[2,2],$
;		;ygap=0.2,$
;		mtitle='!6',$
;		mytitsize=2.5,mytitthick=4,mytitoffset=2.5,mytitle=textoidl('!6\sigma_{int}'),$
;		mxtitsize=2.5,mxtitthick=4,mxtitoffset=2.5,mxtitle=textoidl('!6\gamma')


	xr = [-1,6]
	yr = [0,1]

	xr = [0,4]
	yr = [0.10,0.90]


	;up to 2sigma
	CONTOUR, max_2lnP - grid_2lnP, input_gam,input_sig, levels=[0,2.30,6.17],$
		position=[0.1,0.12,0.83,0.75],/normal,/noerase,$
		thick=thick,charsize=charsize,charthick=charthick,$
		/fill,C_COLORS=[fsc_color('orange'),fsc_color('lime green'),fsc_color('white')] ,$
		xsty=1,ysty=1,$
		;xr=xr,yr=yr,$
		tit='!6',$
		xtit=textoidl('\gamma'),$
		ytit=textoidl('\sigma_{int}')

	CONTOUR, max_2lnP - grid_2lnP, input_gam,input_sig, levels=[0,2.30,6.17],/over;,$
		;thick=1,charsize=charsize,charthick=charthick,$
		;C_LABELS=[0, 1, 1],C_ANNOTATION=['zero',textoidl('1\sigma'),textoidl('2\sigma')]
	;CONTOUR, smooth(max_2lnP - grid_2lnP,5), input_gam,input_sig, levels=[0,0.4],/over;,$
		;C_LABELS=[1,1,1],C_ANNOTATION=['best','68.3%','95.4%'],c_charsize=charsize,c_charthick=charthick
		;C_LABELS=[0,0,0,0],C_ANNOTATION=['best','68.3%','95.4%','99.7%'],c_charsize=charsize,c_charthick=charthick
	;oPLOT,[best_gam],[best_sig],psym=7,symsize=symsize,thick=thick+2,color=fsc_color('black')

;	PLOTs,[best_gam],[best_sig],psym=7,symsize=symsize,thick=thick+2,color=fsc_color('black')
;

;	AL_LEGEND,[$
;		textoidl('uniform prior for \sigma_{int}')],$
;		box=0,/TOP,/RIGHT,charsize=charsize,charthick=charthick,thick=thick




;---------------------------------------------------------
;-marginalizes likelihood function to get 1D prob on each paremeter and determine confidence interval
;---------------------------------------------------------
;-1) gamma
marginal_Pgam = dblarr(n_elements(input_gam))
for i=0L,n_elements(input_gam)-1L do marginal_Pgam[i] = total(exp(0.5D*(grid_2lnP[i,*] - max_2lnP)))  ; = sum(p/pmax)

	maxPgam = max(marginal_Pgam,i_maxPgam)
	print,'==> max marginal Pgamma = ',maxPgam
	print,'==> gam at max marginal Pgamma = ',input_gam[i_maxPgam]
	best_gam_marginal = input_gam[i_maxPgam]
	PLOT,input_gam,marginal_Pgam,psym=10,thick=thick-1,$
		ytit=textoidl('!8p!6(\gamma)'),$
		charsize=charsize-0.5,charthick=charthick,$
		xsty=1,ysty=3,$
		;xr=xr,yr=yr,$
		xtickformat='(A1)',ytickformat='(A1)',xticklen=0.08,yticklen=10D-10,$
		position=[0.1,0.75,0.83,0.98],/normal,/noerase

	PLOTs,[best_gam_marginal,best_gam_marginal],[min(marginal_Pgam),max(marginal_Pgam)],linesty=5,thick=thick,color=fsc_color('red')
	;stop

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

	pCL_per = INTERPOL(int_i,cint_p,CL_per)
	;LINTERP, cint_p, int_i, CL_per, pCL_per
	print,pCL_per

	;for 1sigma error
	gam_1sig = input_gam_sorted[where(marPgam_sorted ge pCL_per[0])]
	print,'best gamma = ',best_gam
	print,'gam error  = ',(max(gam_1sig) - min(gam_1sig))/2D
	print,'gamma-1sig = ',min(gam_1sig)
	print,'gamma+1sig = ',max(gam_1sig)

	XYOUTs,best_gam_marginal+0.02,(max(marginal_Pgam)+min(marginal_Pgam))/2D*0.9D,textoidl('\gamma=')+string(best_gam_marginal,format='(D4.1)')+textoidl('\pm')+string((max(gam_1sig) - min(gam_1sig))/2D,format='(D3.1)'),$
		charsize=charsize-0.5,charthick=charthick+4.5,alignment=0,color=cgCOLOR('white')

	XYOUTs,best_gam_marginal+0.02,(max(marginal_Pgam)+min(marginal_Pgam))/2D*0.9D,textoidl('\gamma=')+string(best_gam_marginal,format='(D4.1)')+textoidl('\pm')+string((max(gam_1sig) - min(gam_1sig))/2D,format='(D3.1)'),$
		charsize=charsize-0.5,charthick=charthick-0.5,alignment=0


;-2) intrinsic scatter
marginal_Psig = dblarr(n_elements(input_sig))
for j=0L,n_elements(input_sig)-1L do marginal_Psig[j] = total(exp(0.5D*(grid_2lnP[*,j] - max_2lnP)))  ; = sum(p/pmax)

	maxPsig = max(marginal_Psig,i_maxPsig)
	print,'==> max marginal Psig = ',maxPsig
	print,'==> sig at max marginal Psig = ',input_sig[i_maxPsig]
	best_sig_marginal = input_sig[i_maxPsig]
	PLOT,marginal_Psig,input_sig,psym=10,thick=thick-1,$
		xtit=textoidl('!8p!6(\sigma_{int})'),$
		charsize=charsize-0.5,charthick=charthick,$
		xsty=3,ysty=1,$
		;xr=xr,yr=yr,$
		xtickformat='(A1)',ytickformat='(A1)',xticklen=10D-10,yticklen=0.06,$
		position=[0.83,0.12,0.98,0.75],/normal,/noerase

	PLOTs,[min(marginal_Psig),max(marginal_Psig)],[best_sig_marginal,best_sig_marginal],linesty=5,thick=thick,color=fsc_color('red')

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

	pCL_per = INTERPOL(int_i,cint_p,CL_per)
	;LINTERP, cint_p, int_i, CL_per, pCL_per
	print,pCL_per

	;for 1sigma error
	sig_1sig = input_sig_sorted[where(marPsig_sorted ge pCL_per[0])]
	print,'best sig = ',best_sig
	print,'sig error  = ',(max(sig_1sig) - min(sig_1sig))/2D
	print,'sig-1sig = ',min(sig_1sig)
	print,'sig+1sig = ',max(sig_1sig)

	XYOUTs,(max(marginal_Psig)+min(marginal_Psig))/2D*0.9D,best_sig_marginal-0.035,textoidl('\sigma_{int}=')+string(best_sig_marginal,format='(D4.2)')+textoidl('\pm')+string((max(sig_1sig) - min(sig_1sig))/2D,format='(D4.2)'),$
		charsize=charsize-0.5,charthick=charthick-0.5,alignment=0.5



;	PLOTs,[best_gam_marginal],[best_sig_marginal],psym=7,symsize=symsize,thick=thick+2,color=fsc_color('red')



	DEVICE,/Close_File
	
	CLEANPLOT,/Silent
		
	SET_PLOT,'X'

;=========================================================


spawn,'epstopdf MC_seleff_flatprior.eps'








	;---------------------------------------------------------
	;eps version

	thick = 4
	charsize = 1.5
	charthick = 3
	symsize = 2

	VSYM,24,thick=thick;,/FILL

	SET_PLOT,'PS'

	!p.multi = 0
	 !x.margin = [6,2.5]
	 !y.margin = [3.0,1.5]
	 ;!P.FONT = -1
	 !x.thick = 3
	 !y.thick = 3
	 !x.ticklen = 0.04
	 !y.ticklen = 0.02

	 !x.margin = [12,2]
	 !y.margin = [8,2]

	DEVICE,/PORTRAIT,/COLOR,BITS=8,/ISOLatin1,$
		xsize=25,ysize=15,/ENCAPSUL ;/Times/Helvetica/Palantino
	;DEVICE,/PORTRAIT,/COLOR,BITS=8,/ISOLatin1,$
	;		xsize=18,xoffset=1.5,ysize=20,yoffset=7.5 ;/Times/Helvetica/Palantino

	DEVICE,FileName='MC_seleff_lognormprior.eps'








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
;sig_int = 0.27D                                                                                                                                        
;esig_int = 0.10D                                                                                                                                       
;4) Bennert+11 local inactive  M-M relation sig_int = 0.38+/-0.1
;sig_int = 0.38D                                                                                                                                        
;esig_int = 0.10D                                                                                                                                       
;5) combined sample of B+11 and HR04 local inactive  M-Mbrl relation sig_int = 0.36+/-0.07
sig_int = 0.36D
esig_int = 0.07D




grid_2lnP_IP = grid_2lnP
for j=0L,n_elements(input_sig)-1L do grid_2lnP_IP[*,j] -= ((input_sig[j]-sig_int)/esig_int)^2  ;-chi2 = 2lnP


help,grid_2lnP_IP
;-find the maximum likelihood point in 2D grid
max_2lnP_IP = max(grid_2lnP_IP, location,/NaN)
i_max = array_indices(grid_2lnP_IP, location)
best_gam = input_gam[i_max[0]]
best_sig = input_sig[i_max[1]]
print, 'max grid index   = ',i_max
print, 'max 2lnP_IP         = ',max_2lnP_IP
print, 'gamma at max P   = ',best_gam 
print, 'sig_int at max P = ',best_sig 





	;up to 2sigma
	;CONTOUR, max_2lnP_IP - grid_2lnP_IP, input_gam,input_sig, levels=[0,2.30,6.17,11.80],$
	CONTOUR, max_2lnP_IP - grid_2lnP_IP, input_gam,input_sig, levels=[0,2.30,6.17],$
		position=[0.1,0.12,0.83,0.75],/normal,/noerase,$
		thick=thick,charsize=charsize,charthick=charthick,$
		;/fill,C_COLORS=[fsc_color('orange'),fsc_color('lime green'),fsc_color('royal blue'),fsc_color('white')] ,$
		/fill,C_COLORS=[fsc_color('orange'),fsc_color('lime green'),fsc_color('white')] ,$
		xsty=1,ysty=1,$
		;xr=xr,yr=yr,$
		tit='!6',$
		xtit=textoidl('\gamma'),$
		ytit=textoidl('\sigma_{int}')

	CONTOUR, max_2lnP_IP - grid_2lnP_IP, input_gam,input_sig, levels=[0,2.30,6.17],/over;,$
		;thick=1,charsize=charsize,charthick=charthick,$
		;C_LABELS=[0, 1, 1],C_ANNOTATION=['zero',textoidl('1\sigma'),textoidl('2\sigma')]
	;CONTOUR, smooth(max_2lnP_IP - grid_2lnP_IP,5), input_gam,input_sig, levels=[0,1.9],/over;,$
		;C_LABELS=[1,1,1],C_ANNOTATION=['best','68.3%','95.4%'],c_charsize=charsize,c_charthick=charthick
		;C_LABELS=[0,0,0,0],C_ANNOTATION=['best','68.3%','95.4%','99.7%'],c_charsize=charsize,c_charthick=charthick
	;oPLOT,[best_gam],[best_sig],psym=7,symsize=symsize,thick=thick+2,color=fsc_color('black')
;	PLOTs,[best_gam],[best_sig],psym=7,symsize=symsize,thick=thick+2,color=fsc_color('black')


	AL_LEGEND,[$
		textoidl('log-normal prior for \sigma_{int}')],$
		box=0,/TOP,/RIGHT,charsize=charsize,charthick=charthick,thick=thick



;---------------------------------------------------------
;-marginalizes likelihood function to get 1D prob on each paremeter and determine confidence interval
;---------------------------------------------------------
;-1) gamma
marginal_Pgam = dblarr(n_elements(input_gam))
for i=0L,n_elements(input_gam)-1L do marginal_Pgam[i] = total(exp(0.5D*(grid_2lnP_IP[i,*] - max_2lnP_IP)))  ; = sum(p/pmax)

	maxPgam = max(marginal_Pgam,i_maxPgam)
	print,'==> max marginal Pgamma = ',maxPgam
	print,'==> gam at max marginal Pgamma = ',input_gam[i_maxPgam]
	best_gam_marginal = input_gam[i_maxPgam]
	PLOT,input_gam,marginal_Pgam,psym=10,thick=thick-1,$
		ytit=textoidl('!8p!6(\gamma)'),$
		charsize=charsize-0.5,charthick=charthick,$
		xsty=1,ysty=3,$
		;xr=xr,yr=yr,$
		xtickformat='(A1)',ytickformat='(A1)',xticklen=0.08,yticklen=10D-10,$
		position=[0.1,0.75,0.83,0.98],/normal,/noerase

	PLOTs,[best_gam_marginal,best_gam_marginal],[min(marginal_Pgam),max(marginal_Pgam)],linesty=5,thick=thick,color=fsc_color('red')

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

	pCL_per = INTERPOL(int_i,cint_p,CL_per)
	;LINTERP, cint_p, int_i, CL_per, pCL_per
	print,pCL_per

	;for 1sigma error
	gam_1sig = input_gam_sorted[where(marPgam_sorted ge pCL_per[0])]
	print,'best gamma = ',best_gam
	print,'gam error  = ',(max(gam_1sig) - min(gam_1sig))/2D
	print,'gamma-1sig = ',min(gam_1sig)
	print,'gamma+1sig = ',max(gam_1sig)

	XYOUTs,best_gam_marginal+0.02,(max(marginal_Pgam)+min(marginal_Pgam))/2D*0.9D,textoidl('\gamma=')+string(best_gam_marginal,format='(D4.1)')+textoidl('\pm')+string((max(gam_1sig) - min(gam_1sig))/2D,format='(D3.1)'),$
		charsize=charsize-0.5,charthick=charthick+4.5,alignment=0,color=cgCOLOR('white')

	XYOUTs,best_gam_marginal+0.02,(max(marginal_Pgam)+min(marginal_Pgam))/2D*0.9D,textoidl('\gamma=')+string(best_gam_marginal,format='(D4.1)')+textoidl('\pm')+string((max(gam_1sig) - min(gam_1sig))/2D,format='(D3.1)'),$
		charsize=charsize-0.5,charthick=charthick-0.5,alignment=0



;-2) intrinsic scatter
marginal_Psig = dblarr(n_elements(input_sig))
for j=0L,n_elements(input_sig)-1L do marginal_Psig[j] = total(exp(0.5D*(grid_2lnP_IP[*,j] - max_2lnP_IP)))  ; = sum(p/pmax)

	maxPsig = max(marginal_Psig,i_maxPsig)
	print,'==> max marginal Psig = ',maxPsig
	print,'==> sig at max marginal Psig = ',input_sig[i_maxPsig]
	best_sig_marginal = input_sig[i_maxPsig]
	PLOT,marginal_Psig,input_sig,psym=10,thick=thick-1,$
		xtit=textoidl('!8p!6(\sigma_{int})'),$
		charsize=charsize-0.5,charthick=charthick,$
		xsty=3,ysty=1,$
		;xr=xr,yr=yr,$
		xtickformat='(A1)',ytickformat='(A1)',xticklen=10D-10,yticklen=0.06,$
		position=[0.83,0.12,0.98,0.75],/normal,/noerase

	PLOTs,[min(marginal_Psig),max(marginal_Psig)],[best_sig_marginal,best_sig_marginal],linesty=5,thick=thick,color=fsc_color('red')


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

	pCL_per = INTERPOL(int_i,cint_p,CL_per)
	;LINTERP, cint_p, int_i, CL_per, pCL_per
	print,pCL_per

	;for 1sigma error
	sig_1sig = input_sig_sorted[where(marPsig_sorted ge pCL_per[0])]
	print,'best sig = ',best_sig
	print,'sig error  = ',(max(sig_1sig) - min(sig_1sig))/2D
	print,'sig-1sig = ',min(sig_1sig)
	print,'sig+1sig = ',max(sig_1sig)


	XYOUTs,(max(marginal_Psig)+min(marginal_Psig))/2D*0.9D,best_sig_marginal-0.035,textoidl('\sigma_{int}=')+string(best_sig_marginal,format='(D4.2)')+textoidl('\pm')+string((max(sig_1sig) - min(sig_1sig))/2D,format='(D4.2)'),$
		charsize=charsize-0.5,charthick=charthick-0.5,alignment=0.5


;	PLOTs,[best_gam_marginal],[best_sig_marginal],psym=7,symsize=symsize,thick=thick+2,color=fsc_color('red')


	DEVICE,/Close_File
	
	CLEANPLOT,/Silent
		
	SET_PLOT,'X'

;=========================================================

spawn,'epstopdf MC_seleff_lognormprior.eps'





;--------------------  extra jobs  -------------------------------{{{
print
etime = (systime(/sec)-itime)/60D
if etime ge 1D then print,'*Elapsed time= ',strtrim(etime,2), ' min (',strtrim(etime/60D,2),' hrs)' $
   else print,'*Elapsed time= ',strtrim(etime*60D,2), ' sec'
print,'*Memory used= ',strtrim((memory(/highwater)-start_mem)/2D^20D,2),' MB'
;profiler,/report,filename='profiler_report.txt' & profiler,/reset
;print
journal
;-----------------------------------------------------------------}}}



STOP,' *** END of CODE *** '

END
