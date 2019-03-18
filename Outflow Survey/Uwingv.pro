function auto_range,spec,peakpos,flat
range=where(spec ge flat,/null)
peak=where(range eq peakpos)
;print,n_elements(range),peak
peak=peak[0]
for i=peak, n_elements(range)-2 do begin
  if range[i+1]-range[i] le 2 then begin
    range_r=range[i]+2
  endif else break
endfor
for i=peak, 1, -1 do begin
  if range[i]-range[i-1] le 2 then begin
    range_l=range[i]-2
  endif else break
endfor
if (size(range_l))[-1] eq 0 then range_l=peakpos-1
if (size(range_r))[-1] eq 0 then range_r=peakpos+1
return,[range_l,range_r]
end

pro Uwingv,region=region,Lvrange=Lvrange

;if ~keyword_set(region) then region='BFS52'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;if ~keyword_set(region) then region='GGMC1'
;if ~keyword_set(Lvrange) then Lvrange=[2,9]
if ~keyword_set(region) then region='region_C_III'
if ~keyword_set(Lvrange) then Lvrange=[-40,-20]
;if ~keyword_set(region) then region='GGMC3'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;if ~keyword_set(region) then region='GGMC4'
;if ~keyword_set(Lvrange) then Lvrange=[-3,5]
;if ~keyword_set(region) then region='lynds'
;if ~keyword_set(Lvrange) then Lvrange=[-3,3]
;if ~keyword_set(region) then region='west'
;if ~keyword_set(Lvrange) then Lvrange=[-1,4]
;if ~keyword_set(region) then region='swallow'
;if ~keyword_set(Lvrange) then Lvrange=[12,18]
;if ~keyword_set(region) then region='horn'
;if ~keyword_set(Lvrange) then Lvrange=[12,18]
;if ~keyword_set(region) then region='remote'
;if ~keyword_set(Lvrange) then Lvrange=[18,28]
rmsU0=0.5
rmsU1=0.9
rmsU2=1.2
cd,'/home/lee/W3/'+region
fits_read,'U_C.fits',datUa,hdrUa
fits_read,'L_C.fits',datLa,hdrLa
fits_read,'Lpeakv0.fits',peakvmap,vhdr

pixU=double(sxpar(hdrUa,'CRPIX3'))
delU=double(sxpar(hdrUa,'CDELT3'))
crvU=double(sxpar(hdrUa,'CRVAL3'))
pixL=double(sxpar(hdrLa,'CRPIX3'))
delL=double(sxpar(hdrLa,'CDELT3'))
crvL=double(sxpar(hdrLa,'CRVAL3'))

specU=reform(datUa[0,0,*])
specL=reform(datLa[0,0,*])
vU=((indgen(n_elements(specU))+1-pixU)*delU+crvU)/1000.0
vL=((indgen(n_elements(specL))+1-pixL)*delL+crvL)/1000.0

szU=size(datUa)
szL=size(datLa)
sz=size(peakvmap)

Urvmap=make_array(sz[1],sz[2],value=0.)
Ubvmap=make_array(sz[1],sz[2],value=0.)

for i=0,sz[1]-1 do begin
  for j=0,sz[2]-1 do begin
    if peakvmap[i,j] ne !values.f_nan then begin
      specU=gauss_smooth(reform(datUa[i,j,*]),1,/edge_mirror)
      specL=gauss_smooth(reform(datLa[i,j,*]),1,/edge_mirror)     
      peakposU=where(vU ge peakvmap[i,j]-0.085 and vU le peakvmap[i,j]+0.085)
      peakposL=where(vL ge peakvmap[i,j]-0.085 and vL le peakvmap[i,j]+0.085);,/null)
      if peakposU[0] ne -1 && peakposL[0] ne -1 then begin
        widthU0=vU[auto_range(specU,peakposU[0],rmsU0)]
        widthU1=vU[auto_range(specU,peakposU[0],rmsU1)]
        widthU2=vU[auto_range(specU,peakposU[0],rmsU2)]
        Urvmap[i,j]=0.3*widthU0[1]+0.3*widthU1[1]+0.4*widthU2[1]-Lvrange[0]
        Ubvmap[i,j]=Lvrange[1]-(0.3*widthU0[0]+0.3*widthU1[0]+0.4*widthU2[0])
        ;Urvmap[i,j]=widthU0[1]-Lvrange[0]
        ;Ubvmap[i,j]=Lvrange[1]-widthU0[0]
      endif
    endif
  endfor
endfor

Ubvmap=congrid(Ubvmap,3*n_elements(Ubvmap[*,0]),3*n_elements(Ubvmap[0,*]),/MINUS_ONE);,cubic=-0.5)
Ubvmap=smooth(Ubvmap,[3,3],/edge_mirror)
Ubvmap=round(Ubvmap/0.2)*0.2
Urvmap=congrid(Urvmap,3*n_elements(Urvmap[*,0]),3*n_elements(Urvmap[0,*]),/MINUS_ONE);,cubic=-0.5)
Urvmap=smooth(Urvmap,[3,3],/edge_mirror)
Urvmap=round(Urvmap/0.2)*0.2
mkhdr,rhdr,Urvmap
sxaddpar,rhdr,'CTYPE1',sxpar(vhdr,'CTYPE1')
sxaddpar,rhdr,'CRVAL1',sxpar(vhdr,'CRVAL1')
sxaddpar,rhdr,'CDELT1',sxpar(vhdr,'CDELT1')/3.
sxaddpar,rhdr,'CRPIX1',sxpar(vhdr,'CRPIX1')*3.
sxaddpar,rhdr,'CTYPE2',sxpar(vhdr,'CTYPE2')
sxaddpar,rhdr,'CRVAL2',sxpar(vhdr,'CRVAL2')
sxaddpar,rhdr,'CDELT2',sxpar(vhdr,'CDELT2')/3.
sxaddpar,rhdr,'CRPIX2',sxpar(vhdr,'CRPIX2')*3.

fits_write,'Urvmap.fits',Urvmap,rhdr
fits_write,'Ubvmap.fits',Ubvmap,rhdr

end