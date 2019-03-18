pro Uwingvpeak,region=region,Lvrange=Lvrange

;if ~keyword_set(region) then region='BFS52'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;distance=2
;if ~keyword_set(region) then region='GGMC1'
;if ~keyword_set(Lvrange) then Lvrange=[2,9]
;distance=2
if ~keyword_set(region) then region='region_C_III'
if ~keyword_set(Lvrange) then Lvrange=[-34,-28.5]
distance=2.0
;if ~keyword_set(region) then region='GGMC3'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;distance=2
;if ~keyword_set(region) then region='GGMC4'
;if ~keyword_set(Lvrange) then Lvrange=[-3,5]
;distance=2
;if ~keyword_set(region) then region='lynds'
;if ~keyword_set(Lvrange) then Lvrange=[-3,3]
;distance=0.4
;if ~keyword_set(region) then region='west'
;if ~keyword_set(Lvrange) then Lvrange=[-1,4]
;distance=0.6
;if ~keyword_set(region) then region='swallow'
;if ~keyword_set(Lvrange) then Lvrange=[12,18]
;distance=3.8
;if ~keyword_set(region) then region='horn'
;if ~keyword_set(Lvrange) then Lvrange=[12,18]
;distance=3.8
;if ~keyword_set(region) then region='remote'
;if ~keyword_set(Lvrange) then Lvrange=[18,28]
;distance=8.5
cd,'/home/lee/W3/'+region
fits_read,'Ubvmap.fits',Ubvmap,bhdr
fits_read,'Urvmap.fits',Urvmap,rhdr
;fits_read,region+'_Ua_C.fits',datUa,hdrUa
;fits_read,region+'_La_mask.fits',datLa,hdrLa
;fits_read,region+'_L2a_mask.fits',datL2a,hdrL2a
;dathdrL=list(datLa,hdrLa)
;dathdrL2=list(datL2a,hdrL2a)
;Lvrange=[-0.3,10]
;L2vrange=[2,9]

conf13='fellwalker2d.conf'
openw,c13,/get_lun,conf13
vjump=round(18./distance)>3
vjump=vjump<20.
;print,vjump
;vjump=round(12./distance)>4
;vjump=vjump<10
printf,c13,'FellWalker.AllowEdge=0'
printf,c13,'FellWalker.CleanIter=1'
printf,c13,'FellWalker.FlatSlope=2*RMS';2*RMS'
printf,c13,'FellWalker.FwhmBeam=9'
printf,c13,'FellWalker.MaxBad=0.1'
printf,c13,'FellWalker.MaxJump='+num2str(vjump,format='(i0)')
printf,c13,'FellWalker.MinHeight=3';10*RMS'
printf,c13,'FellWalker.MinPix=13'
printf,c13,'FellWalker.MinDip=0.5*RMS';.5*RMS'
printf,c13,'FellWalker.Noise=2'
;printf,c13,'FellWalker.VeloRes=0'
printf,c13,'FellWalker.RMS=0.17'
free_lun,c13

;convert fits to ndf; starlink-2015A
bfile_fits='Ubvmap.fits'
bfile_sdf='Ubvmap.sdf'
bfile_mask='Ubvmap_mask.sdf'
bpeaks_cat='Ubvmap_peaks.fit'
bcmdfile='Ubvmap_cmd'
openw,bcmd,/get_lun,bcmdfile
;printf,cmd,'export STARLINK_DIR=/usr/local/astrosoft/star-2015A'
printf,bcmd,'export STARLINK_DIR=/home/lee/star-2015B'
printf,bcmd,'source $STARLINK_DIR/etc/profile'
printf,bcmd,'convert'
printf,bcmd,'cupid'
bconvert_str='fits2ndf '+bfile_fits+' '+bfile_sdf
bfindclumps_str="findclumps in="+bfile_fits+" deconv=no method=fellwalker out="+bfile_mask+" outcat="+bpeaks_cat+" config=""'^"+conf13+"'"" wcspar=yes"+" SHAPE=""Ellipse"""
printf,bcmd,bfindclumps_str
free_lun,bcmd
spawn,'chmod +x '+bcmdfile
spawn,'./'+bcmdfile

;read catalog in wcs units pix*cdelt
;peak[1|2] cen[1|2]: galactic degree
;size[1|2]: arcsec
ftab_ext,bpeaks_cat,[1,2,3,4,5,6,7,9,11],$
  bidx,bpeak1,bpeak2,bcen1,bcen2,bsize1,bsize2,bpeak,bshape
openw,bcat,/get_lun,'bluepeaks.txt'
openw,belli,/get_lun,'bshape.txt'
printf,bcat,'#num','peak1','peak2','cen1','cen2','size1','size2','Bwing',$
  format='(a4,2(2a9),3a6)' 
printf,bcat,'#   ','[deg]','[deg]','[deg]','[deg]','['']','['']','[km/s]',$
  format='(a4,2(2a9),3a6)'

bcount=0
for i=0, n_elements(bidx)-1 do begin
  ;if min([bsize1[i],bsize2[i]]) gt 0.4*60. then begin
    bcount++
    printf,bcat,bcount,bpeak1[i],bpeak2[i],bcen1[i],bcen2[i],bsize1[i]/60.,bsize2[i]/60,bpeak[i],$
    format='(i4,2(2f9.3),3f6.1)'
    printf,belli,bcount,bshape[i],format='(i4,x,a0)'
  ;endif
endfor
delvar,bcat,bcount,bpeak1,bpeak2,bcen1,bcen2,bsize1,bsize2,bpeak

free_lun,bcat,belli

readcol,'bluepeaks.txt',bnum,bpeak1,bpeak2,bcen1,bcen2,bsize1,bsize2,bpeak,$
    format='I,F,F,F,F,F,F,F',stringskip='#',/silent
readcol,'bshape.txt',bnum,bstr1,bstr2,bstr3,bepos1,bepos2,bmajor,bminor,bposangle,format='I,A,A,A,F,F,F,F,F',stringskip='#',/silent
bepos1[where(bepos1 lt 0)]+=360.
openw,bpcat,/get_lun,'bluepeaks.cat'
printf,bpcat,'#num','peak1','peak2','cen1','cen2','major','minor','posang','size1','size2','Bwing',$
  format='(a4,2(a9,a7),2a6,a7,3a6)' 
bcount=0
for i=0, n_elements(bnum)-1 do begin
  if min([bsize1[i],bsize2[i]]) ge 0.1 then begin
    bcount++
    printf,bpcat,bcount,bpeak1[i],bpeak2[i],bcen1[i],bcen2[i],bmajor[i]<5./60,bminor[i],bposangle[i],bsize1[i],bsize2[i],bpeak[i],$
      format='(i4,2(f9.3,f7.3),2f6.3,i7,3f6.1)'
  endif
endfor
free_lun,bpcat;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;redlobe;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
rfile_fits='Urvmap.fits'
rfile_sdf='Urvmap.sdf'
rfile_mask='Urvmap_mask.sdf'
rpeaks_cat='Urvmap_peaks.fit'
rcmdfile='Urvmap_cmd'
openw,rcmd,/get_lun,rcmdfile
;printf,cmd,'export STARLINK_DIR=/usr/local/astrosoft/star-2015A'
printf,rcmd,'export STARLINK_DIR=/home/lee/star-2015B'
printf,rcmd,'source $STARLINK_DIR/etc/profile'
printf,rcmd,'convert'
printf,rcmd,'cupid'
rconvert_str='fits2ndf '+rfile_fits+' '+rfile_sdf
rfindclumps_str="findclumps in="+rfile_fits+" deconv=no method=fellwalker out="+rfile_mask+" outcat="+rpeaks_cat+" config=""'^"+conf13+"'"" wcspar=yes"+" SHAPE=""Ellipse"""
printf,rcmd,rfindclumps_str
free_lun,rcmd
spawn,'chmod +x '+rcmdfile
spawn,'./'+rcmdfile

ftab_ext,rpeaks_cat,[1,2,3,4,5,6,7,9,11],$
  ridx,rpeak1,rpeak2,rcen1,rcen2,rsize1,rsize2,rpeak,rshape
openw,rcat,/get_lun,'redpeaks.txt'
openw,relli,/get_lun,'rshape.txt'
printf,rcat,'#num','peak1','peak2','cen1','cen2','size1','size2','Rwing',$
  format='(a4,2(2a9),3a6)' 
printf,rcat,'#   ','[deg]','[deg]','[deg]','[deg]','['']','['']','[km/s]',$
  format='(a4,2(2a9),3a6)'

rcount=0
for i=0, n_elements(ridx)-1 do begin
  ;if min([rsize1[i],rsize2[i]]) gt 0.4*60. then begin
    rcount++
    printf,rcat,rcount,rpeak1[i],rpeak2[i],rcen1[i],rcen2[i],rsize1[i]/60.,rsize2[i]/60,rpeak[i],$
    format='(i4,2(2f9.3),3f6.1)'
    printf,relli,rcount,rshape[i],format='(i4,x,a0)'
  ;endif
endfor
delvar,rcat,rcount,rpeak1,rpeak2,rcen1,rcen2,rsize1,rsize2,rpeak

free_lun,rcat,relli

readcol,'redpeaks.txt',rnum,rpeak1,rpeak2,rcen1,rcen2,rsize1,rsize2,rpeak,$
    format='I,F,F,F,F,F,F,F',stringskip='#',/silent
readcol,'rshape.txt',rnum,rstr1,rstr2,rstr3,repos1,repos2,rmajor,rminor,rposangle,format='I,A,A,A,F,F,F,F,F',/silent
repos1[where(repos1 lt 0)]+=360.
openw,rpcat,/get_lun,'redpeaks.cat'
printf,rpcat,'#num','peak1','peak2','cen1','cen2','major','minor','posang','size1','size2','Rwing',$
  format='(a4,2(a9,a7),2a6,a7,3a6)' 
rcount=0
for i=0, n_elements(rnum)-1 do begin
  if min([rsize1[i],rsize2[i]]) ge 0.1 then begin
    rcount++
    printf,rpcat,rcount,rpeak1[i],rpeak2[i],rcen1[i],rcen2[i],rmajor[i]<5./60,rminor[i],rposangle[i],rsize1[i],rsize2[i],rpeak[i],$
      format='(i4,2(f9.3,f7.3),2f6.3,i7,3f6.1)'
  endif
endfor
free_lun,rpcat;,rscat


end