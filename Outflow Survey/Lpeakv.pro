pro lpeakv,region=region,Lvrange=Lvrange

if ~keyword_set(region) then region='region_C_III'
if ~keyword_set(Lvrange) then Lvrange=[-34,-28.5]
cd,'/home/lee/W3/'+region
;fits_read,region+'_Ua_C.fits',datUa,hdrUa 
fits_read,'L_C_III.fits',datLa,hdrLa
cpx=double(sxpar(hdrLa,'CRPIX3'))
del=double(sxpar(hdrLa,'CDELT3'))
crv=double(sxpar(hdrLa,'CRVAL3'))
szL=size(datLa)
peakvmap=make_array(szL[1],szL[2],3,value=!values.F_NAN)
spec=datLa[0,0,*]
v=((indgen(n_elements(spec))+1-cpx)*del+crv)/1000.0
peakvpos=[0]
for i=0,szL[1]-1 do begin
  for j=0,szl[2]-1 do begin
    spec=mean(datLa[(i-1)>0:(i+1)<(szL[1]-1),(j-1)>0:(j+1)<(szL[2]-1),*],dim=1)
    spec=mean(spec,dim=1)
    spec=reform(spec)
    peakvpos=findpeakpos(spec,sm=3,limit=0.3)
    if peakvpos[0] ne 0 then begin
    peakv=v[peakvpos]
    ;for test
    pos0=where(peakv gt Lvrange[0] and peakv lt Lvrange[1],/null)
    if n_elements(pos0) ge 1 then begin
      ;print,n_elements(pos0)
      peakv0=(peakv[pos0])[0]
      peakvmap[i,j,0]=peakv0
    endif
;    pos1=where(peakv ge 5 and peakv lt 11,/null)
;    if n_elements(pos1) ge 1 then begin
;      peakv1=(peakv[pos1])[0]
;      peakvmap[i,j,1]=peakv1
;    endif
;    pos2=where(peakv ge 11 and peakv lt 30,/null)
;    if n_elements(pos2) ge 1 then begin
;      peakv2=(peakv[pos2])[0]
;      peakvmap[i,j,2]=peakv2
;    endif
;    peakv=peakv[0]
;    peakvmap[i,j,0]=peakv
    endif
  endfor
endfor

mkhdr,hdr,reform(peakvmap[*,*,0])
sxaddpar,hdr,'NAXIS1',sxpar(hdrLa,'NAXIS1')
sxaddpar,hdr,'CTYPE1',sxpar(hdrLa,'CTYPE1')
sxaddpar,hdr,'CRVAL1',sxpar(hdrLa,'CRVAL1')
sxaddpar,hdr,'CDELT1',sxpar(hdrLa,'CDELT1')
sxaddpar,hdr,'CRPIX1',sxpar(hdrLa,'CRPIX1')
sxaddpar,hdr,'NAXIS2',sxpar(hdrLa,'NAXIS2')
sxaddpar,hdr,'CTYPE2',sxpar(hdrLa,'CTYPE2')
sxaddpar,hdr,'CRVAL2',sxpar(hdrLa,'CRVAL2')
sxaddpar,hdr,'CDELT2',sxpar(hdrLa,'CDELT2')
sxaddpar,hdr,'CRPIX2',sxpar(hdrLa,'CRPIX2')
sxaddpar,hdr,'UNIT','VELOCITY(km/s)'
fits_write,'Lpeakv0.fits',reform(peakvmap[*,*,0]),hdr
;fits_write,'Lpeakv1.fits',reform(peakvmap[*,*,1]),hdr
;fits_write,'Lpeakv2.fits',reform(peakvmap[*,*,2]),hdr

end
    