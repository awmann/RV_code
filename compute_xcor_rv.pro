PRO compute_xcor_rv,wave, $              ; vavelength array
                    flux, $              ; flux array 
                    spt, $               ; spectral type for auto-find template. This requires other programs and data that are not on github... so ignore this
                    rv, $                ; output RV
                    wavelim=wavelim, $   ; 2 element array specifying the wavelength region to search around
                    norm_reg=norm_reg, $ ; 2 element array specifying the normalization region
                    twave=twave, $       ; Template wavelength array. Supplying this will cause the program to ignore the spt variable
                    tflux=tflux, $       ; template flux
                    maxshift=maxshift, $ ; maximum allowed shift (IN PIXELS)
                    startshift=startshift, $ ; initial shift, sometimes we use this for applying barycentric and other corrections to the data before doing the RV, but I suggest letting it sit at 0 (the default)
                    showplot=showplot, $     ; plot the chi^2 distribution as a function of the offset above the spectrum and shifted spectrum
                    redo=redo ; set to 1 to turn off the redo module
;Eric J. Hilton
;modified 20120307
;enable user to manually adjust fit a bit, so see if it locks on to a
;better answer
;also allows for more plotting

;20111015
;compute the radial velocity using cross correlation to templates. 
;this pro is mostly general, although is meant for UH22 SNIFS or MDM spectra.
;(can be easily adapted for other sets)


;inputs: wave (in linear dispersion), flux, spectral type
;keywords: 
;wavelim = [x1,x2] is the wavelength region used for the correlation (in angstroms)
;default is 5300, 8500]
;norm_reg = [x1,x2] is wavelength region used to normalize the
;spectrum and the template spectrum (in angs) default is [6450,6650]
;twave and tflux - if present, then don't read in the template
;this allows for faster computation. the wrapper can read them all in
;once, and just pass the relevant one when necessary.
;note that in this case, the spt variable is not necessary

  oldwave = wave
  oldflux = flux

  c = 299792.458D
  if n_elements(showplot) eq 0 then showplot=0
  if n_elements(startshift) eq 0 then startshift=0

  wave = wave/(startshift/c + 1.)

  if not(keyword_set(wavelim)) then wavelim = [5300, 8500]
  if not(keyword_set(norm_reg)) then norm_reg = [7150,7350]

  if n_elements(maxshift) eq 0 then maxshift=300

  ;;read in Bochanski's template
  if not(keyword_set(tflux)) then $
     readtemplatespec, spt, twave,tflux,/silent ;; the code readtemplatespec is missing from the


  ;;normalize both template and spectrum
  tflux = tflux/ median(tflux[where(twave gt norm_reg[0] and twave lt norm_reg[1])])
  flux = flux/ median(flux[where(wave gt norm_reg[0] and wave lt norm_reg[1])])

  ind=where(wave gt wavelim[0] and wave lt wavelim[1])
  ind2=where(twave gt wavelim[0] and twave lt wavelim[1])

  ; CREATE WAVELENGTH ARRAY THAT IS EVENLY SPACED IN LOG SPACE
  n_bins = double((size(wave[ind]))[1])
  ln_wave = double(alog(wave[ind]))
  ln_range = double(max(ln_wave)-min(ln_wave))
  ln_grid = findgen(n_bins)*ln_range/(n_bins-1) + min(ln_wave)
  newwave = exp(1)^ln_grid

  tnewflux=interpol(tflux,twave,newwave)
  newflux=interpol(flux,wave,newwave)

  ;;Find the velocity per pixel
  ind_med = where(newwave eq median(newwave))
  velpix = 299792.458D*(newwave[ind_med] - newwave[ind_med-1])/newwave[ind_med-1]  
  
  minvel=min([velpix,velpix])

  ;;create a new wavelength array in log-space
  ;;so the pixel shift is constant in velocity space
  ;;ca=alog10(double(min(wave[ind])))
  ;;cb=alog10(double(min(wave[ind]))+minvel/2.99792458d5*min(wave[ind]))-ca
  ;;newwave=10^(ca+cb*dindgen(100000))
  ;;newwave = newwave[where(newwave gt wavelim[0] and newwave lt wavelim[1])]


  ;;do the actual cross correlation  
  if showplot eq 0 then begin
     am_xcorl, newflux,tnewflux,maxshift,sft,xp,chi
     rv= sft*minvel
     if n_elements(redo) eq 0 then begin
        if abs(rv) gt 1000 then begin
           nwavelim = [5300, 7000]
           nind=where(wave gt nwavelim[0] and wave lt nwavelim[1])
           nind2=where(twave gt nwavelim[0] and twave lt nwavelim[1])
           am_xcorl, newflux,tnewflux,maxshift,sft,xp,chi
           rv= sft*minvel
           if abs(rv) gt 1000 then begin
              nwavelim = [7700, 9000]
              nind=where(wave gt nwavelim[0] and wave lt nwavelim[1])
              nind2=where(twave gt nwavelim[0] and twave lt nwavelim[1])
              am_xcorl, newflux,tnewflux,maxshift,sft,xp,chi
              rv= sft*minvel
           endif
        endif
        if abs(rv) gt 1000 then rv = 0.0 ;; do no harm algorithm
     endif
  endif else begin 
     !p.multi=[0,1,2]
     loadct,39,/silent
     am_xcorl, newflux,tnewflux,maxshift,sft,xp,chi,/pl
     rv= sft*minvel
     plot, newwave, tnewflux,xrange=norm_reg,/xsty
     oplot, newwave, newflux,color=cgcolor('red')
     oplot, newwave/(rv/c+1.),newflux,color=cgcolor('blue')  
     !p.multi=0
  endelse
  wave = oldwave
  flux = oldflux

end
