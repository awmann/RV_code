;+
; NAME:  
;       AM_xcorl
; PURPOSE:
;       Measure the offset in pixel units between two arrays. Good to
;       sub-pixel (maybe 0.1 pixels, roughly) shifts. Very useful for
;       measuring the Radial Velocity shift between two spectra, such
;       as an observed and template.
;
;       If you want an RV the array needs in velocity space. Eric
;       Hilton wrote compute_xcor_rv to help with this.
;
;       If you get an error "did not converge", this means there
;       was no obvious absolute minimum (e.g., a chi^2 that is
;       constantly decreasing with shift). This can happen for lots of
;       reasons, but the most common is that there are not enough
;       unambiguous features in the spectral region. Another common
;       problem is that maxshift is too small.
;
;
; CALLING SEQUENCE: 
;       am_xcorl,star,template,maxshift,shft,xpout,chiout
;
; INPUT PARAMETERS: 
;       star       spectrum in arbitrary units. 
;       template   template spectrum in the same units as star
;       maxshift   the maximum allowed shift (in pixels) allowed
;
;
; OUTPUT PARAMETERS: 
;       shft       final computed shift in pixels 
;       xpout      x-array (shift numbers) corresponding to chiout
;       chiout     chi^2 values for the xpout
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       plot     Plot the original and slightly shifted spectrum 
;   
; EXAMPLES:
;       
;
; RESTRICTIONS:
;
;
; DEPENDENCIES:
;       
; HISTORY: 
;      I wrote this based on another routine, xcorl. It's not
;      entirely clear who wrote the original... I've
;      significantly edited it for our use.
;      Rewritten for github A. Mann 04/07/2014
;-

Pro  am_xcorl,star, $           ; spectrum of the star, must have even (wavelength) separations in log-space (so the shift is constant in RV space)
              template, $       ; template spectrum, must be a pixel-for-pixel match (wavelength) with star
              maxshift, $       ; maximum pixel shift allowed
              shft, $           ; the shift, in pixels
              xpout, $          ;
              chiout, $         ;
              plot=pl           ;

  if n_params() lt 3 then begin   
     print,'am_xcorl,star,template,maximum shift'
     return 
  endif
  if n_elements(star) ne n_elements(template) then begin
     print,'The star and the template should have the same number of pixels'
     return
  endif
  
  len = n_elements(template)
  if maxshift gt (len-1)/2 then begin
     print,'Specificed maximum shift is impossible, switching to maxshift = '+strcompress((len-1)/2)
     maxshift = (len-1)/2
  endif

  newln = len - 2*maxshift        
  te = template/(total(template)/len) ; normalization
  st = star/(total(star)/len)         ; normalization
  newend = maxshift + newln - 1
  x =findgen(2 * maxshift+1) - maxshift
  chi = fltarr(2 * maxshift+1)
  
  for j = -maxshift,maxshift do begin 
     dif = te(maxshift:newend) - st(maxshift+j:newend+j)
     chi(j+maxshift) = total(dif*dif)
  endfor
  xcr = chi
  if keyword_set(pl) then begin
     npt=n_elements(chi)
     xp=indgen(npt)-npt/2
     if pl eq 2 then !p.multi=[0,1,2]
     plot,xp,chi,/ynozero,ps=-4
     xpout = xp
     chiout = chi
     if pl eq 2 then !p.multi=0
  endif
  len = n_elements(x) * 100
  xl = findgen(len)
  xl = xl/100. - maxshift
  xp = xl(0:len-100)
  cp = spline(x,chi,xp)
  minchi=min(cp,mm)
  shft = xp(mm(0))

  if abs(shft) ge maxshift then begin
     print,'fit did not converge, setting shift to absurdly high number.'
     shift = 1d10
     return 
  endif


  nf=21.
  rf=10.
  nc=-1.
  fchi=fltarr(nf)
  xl=fltarr(nf)
  for j=-rf,rf do begin
     xl(nc+1)=shft+j/10.
     am_shfour,st,-xl(nc+1),nst 
     nc=nc+1
     if keyword_set(mult) then begin
        dif=nst(maxshift:newend)*te(maxshift:newend)
        fchi(nc)=total(abs(dif))
     endif else begin
        dif=nst(maxshift:newend)-te(maxshift:newend)
        fchi(nc)=total(dif*dif)
     endelse
  endfor
  xp=findgen((nf-1)*100+1)/(1000.)+shft-1
  cp=spline(xl,fchi,xp)
  if keyword_set(mult) then minchi = max(cp,mm) else minchi=min(cp,mm)
  fshft=xp(mm(0))
  if keyword_set(pl) then begin
     if max(cp) gt max(chi) and keyword_set(mult) then cp=cp*max(chi)/max(cp)
     oplot,xp,cp,ps=0,co=2
  endif
  shft=fshft
  
end

