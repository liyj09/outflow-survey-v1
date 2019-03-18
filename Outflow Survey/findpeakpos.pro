function findpeakpos,dat,sm=sm,limit=limit
if ~keyword_set(sm) then sm=3
if ~keyword_set(limit) then limit=1
sz=size(dat)
if sz[0] eq 1 && sz[1] ge sm then begin
  ;print,sz
  peakpos=[0]
  smdat=gauss_smooth(dat,sm,/edge_mirror)>limit
  d_dat=gauss_smooth(deriv(smdat),sm*2+1,/edge_mirror)
  ;dd_dat=smooth(deriv(d_dat),5,/edge_mirror)
  count=0
  for i=1,sz[1]-2 do begin
    if d_dat[i] gt 0 && d_dat[i+1] lt 0 then begin 
    ;if smdat[i] gt smdat[i-1] && smdat[i] gt smdat[i+1] then begin
      peakpos=[peakpos,i]
      count++
    ;endif
    endif
    ;pos=where(d_dat gt 0)
  endfor
  if count ge 1 then begin
    peakpos=peakpos[1:-1]
    if count gt 1 then begin
    for i=0, n_elements(peakpos)-2 do begin
      if peakpos[i+1]-peakpos[i] le round(sm/2.)+1 then begin
        newpos=round(mean(peakpos[i:i+1]))
        peakpos[i]=newpos
        peakpos[i+1]=newpos
      endif
    endfor
    peakpos=peakpos[uniq(peakpos)]
    endif
  endif
endif

return,peakpos

end