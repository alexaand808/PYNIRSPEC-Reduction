PRO convsamp,rfmspec,xout=xout,yout=yout,specres=specres,sampling=sampling,$
             psc=psc,xrange=xrange,sh=sh,ksh=ksh,doplot=doplot
@natconst
  
  rfmrd,rfmspec,inx,iny

  IF NOT KEYWORD_SET(psc) THEN psc = 1.
  IF NOT KEYWORD_SET(sh) THEN sh  = 0
  IF NOT KEYWORD_SET(ksh) THEN ksh  = 0
  
  fwhm   = specres*inx[0]/(cc*1d-5*(inx[1]-inx[0]))
  nk     = fwhm*5.
  sp_KERNEL = exp(-(findgen(nk)-(nk-1)/2.-ksh/2.)^2.*2.77/fwhm^2)+exp(-(findgen(nk)-(nk-1)/2.+ksh/2.)^2.*2.77/fwhm^2)
  sp_KERNEL = sp_KERNEL/TOTAL(sp_KERNEL)
  iny_conv  = CONVOL(iny,sp_KERNEL,TOTAL(sp_KERNEL),/edge_truncate)        
  
;  obs=MRDFITS('/Users/pontoppi/WORK/OBSPROGS/CRIRES/REDUCTIONS/2008-08-04/WAVECAL/BS6084_200178927_4730.0_0_0_wav.fits',3)
  obs=MRDFITS('/Users/pontoppi/WORK/OBSPROGS/CRIRES/REDUCTIONS/2007-08-29/WAVECAL/BS6084_200168536_4730.0_0_0_wav.fits',3)
  obs.flux = obs.flux/max(obs.flux)

  IF KEYWORD_SET(doplot) THEN BEGIN
     tvlct, [[0,200,0,0],[0,0,200,0],[0,0,0,200]]
     @plot_setup
     set_plot, 'ps'
     device, filename='trans_model.eps',/encapsulated, /color,xsize=17,ysize=10
;  !p.charsize=2
;  window, 1, ysize=800,xsize=1200
     !p.multi=[0,1,2]
     obs.flux = obs.flux*psc
     plot, obs.wavelength, obs.flux,psym=10,xrange=xrange,/xs,ymargin=[1,1],ytitle='!6transmission'
     oplot, 1d7/inx-sh, iny_conv,color=1,linestyle=0, thick=4
     
     plot, obs.wavelength, obs.flux/INTERPOL(iny_conv,1d7/inx-sh,obs.wavelength),yrange=[0.8,1.2],xrange=xrange,/xs,psym=10,$
           ytitle='!6standard/model',xtitle='!6nm'
;  plot, obs.wavelength, obs.flux-INTERPOL(iny_conv,1d7/inx-sh,obs.wavelength),yrange=[-0.2,0.2],xrange=xrange,/xs,psym=10
     device, /close
     set_plot, 'x'
     @plot_clear
  ENDIF
  xout = 1d4/inx
  yout = iny_conv
END
