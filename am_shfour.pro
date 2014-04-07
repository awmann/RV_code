;+
; NAME:  
;       AM_shfour
; PURPOSE:
;        Shifts a spectrum by a small user-specified (non-integer) amount.
;
;
; CALLING SEQUENCE: 
;       am_shfour,sp,shift,newsp,plot=plot
;
; INPUT PARAMETERS: 
;       sp       spectrum in arbitrary units
;       shift    shift in pixels (array index units)
;
; OUTPUT PARAMETERS: 
;       newsp    shifted spectrum
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
;      I wrote this based on another routine, shfour. It's not
;      entirely clear who wrote the original...
;      Rewritten for github A. Mann 04/07/2014
;-

pro am_shfour, sp,shift,newsp,plot=plot

  if n_params() lt 3 then begin
     print, 'please provide spectrum, shift, and a variable to store the new spectrum'  
     return
  endif
  ln = n_elements(sp)
  nsp = sp
  fourtr = fft(nsp,1)
  sig = findgen(ln)/ln-.5
  sig = shift(sig,(ln/2))
  sh = sig*2. *!pi*shift
  shfourtr = complex(cos(sh),sin(sh)) * fourtr
  newsp = fft(shfourtr,-1)
  newsp = float(newsp(0:ln-1))
  if keyword_set(plot) then begin
     plot,sp,yrange=[min(sp)-0.5,max(sp)]
     oplot,newsp-.5
  endif
  return

end
