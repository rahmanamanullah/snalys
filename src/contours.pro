PRO contours,fname,omlow,omhigh,nom,oxlow,oxhigh,nox,hard=hard
readcol,fname,tmp
omvec = findgen(nom)*(omhigh-omlow)/(nom-1) + omlow
oxvec = findgen(nox)*(oxhigh-oxlow)/(nox-1) + oxlow
chi2  =fltarr(nom,nox)

n=long(0)
FOR I=0,nom-1 DO BEGIN
 FOR J=0,nox-1 DO BEGIN
;  chi2(j,i) = tmp(n)
   chi2(i,j) = tmp(n)
  n = n + 1
 ENDFOR
ENDFOR

chi2 = chi2 - min(chi2)

level=[-2.3/2.]
thick=[5]
loadct,4
colors=[220]

!x.margin=13
!y.margin=5
IF keyword_set(hard) THEN BEGIN
set_plot,'ps'
device,bits_per_pixel=8, /color, filename='multi.ps',$
SET_FONT = 'Palatino-Roman'
ENDIF
xtitle='!7X!D!6M'
;ytitle='!7X!DK'
ytitle='!8 w!D0'
contour,-chi2,omvec,oxvec,levels=level,c_colors=colors,c_thick=thick,/fill,$
    xcharsize=1.5,ycharsize=1.5,xtitle=xtitle,ytitle=ytitle
colors=[0]
IF keyword_set(hard) THEN BEGIN
device,/close
set_plot,'x'
ENDIF
RETURN
END