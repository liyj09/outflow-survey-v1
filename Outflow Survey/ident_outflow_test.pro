function auto_wing,spec,peakpos,flat
range=where(spec ge flat,/null)
peak=where(range eq peakpos,/null)
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
return,[range_l,range_r]
end

function multi_check,peak3,specU,vU,FWHMU,wingU,wingL,rmsU,rmsL,log,pixsz
pixsz=(2.*pixsz+1)/3.0
B_bad=0
R_bad=0
B_wing=[wingU[0],wingL[0]]
R_wing=[wingL[1],wingU[1]]
if wingL[0]-wingU[0] le 1.0 then begin 
  B_bad++
  printf,log,' blue wing too narrow!'
endif
if wingU[1]-wingL[1] le 1.0 then begin
  R_bad++
  printf,log,' red wing too narrow!'
endif
;where(vL ge peak3[i]-0.1 and vL le peak3[i]+0.1,/null)
Fpos0=where(vU ge FWHMU[0]-0.085 and vU le FWHMU[0]+0.085,/null) & Fpos0=Fpos0[0]
Fpos1=where(vU ge FWHMU[1]-0.085 and vU le FWHMU[1]+0.085,/null) & Fpos1=Fpos1[0]
Bpos0=where(vU ge B_wing[0]-0.085 and vU le B_wing[0]+0.085,/null) & Bpos0=Bpos0[0]
Bpos1=where(vU ge B_wing[1]-0.085 and vU le B_wing[1]+0.085,/null) & Bpos1=Bpos1[0]
Rpos0=where(vU ge R_wing[0]-0.085 and vU le R_wing[0]+0.085,/null) & Rpos0=Rpos0[0]
Rpos1=where(vU ge R_wing[1]-0.085 and vU le R_wing[1]+0.085,/null) & Rpos1=Rpos1[0]
iicore=total(specU[Fpos0:Fpos1])
iib=total(specU[Bpos0:Bpos1])
iir=total(specU[Rpos0:Rpos1])
Fspec=specU[Fpos0:Fpos1]
Bspec=specU[Bpos0:Bpos1]
Bspec=Bspec[where(Bspec gt rmsL)]
if n_elements(Bspec) ge 5 then Bspec=smooth(Bspec[0:-2],3,/edge_mirror)
Rspec=specU[Rpos0:Rpos1]
Rspec=Rspec[where(Rspec gt rmsL)]
if n_elements(Rspec) ge 5 then Rspec=smooth(Rspec[1:-1],3,/edge_mirror)
Binter=interpol([Bspec[0],Bspec[-1]],n_elements(Bspec))
Rinter=interpol([Rspec[0],Rspec[-1]],n_elements(Rspec))
dB=Binter-Bspec
dR=Rinter-Rspec
if peak3-B_wing[1] gt 1.75*(FWHMU[1]-FWHMU[0]) then begin
  B_bad++
  printf,log,' blue wing too far!'
endif
if R_wing[0]-peak3 gt 1.75*(FWHMU[1]-FWHMU[0]) then begin
  R_bad++
  printf,log,' red wing too far!'
endif
if (iib/iicore gt 0.6) || (n_elements(where(dB lt -2.0*pixsz*rmsU,/null)) gt 1) then begin
  B_bad++
  printf,log,' blue wing contaminated!'
endif
if (iir/iicore gt 0.6) || (n_elements(where(dR lt -2.0*pixsz*rmsU,/null)) gt 1) then begin
  R_bad++
  printf,log,' red wing contaminated!'
endif
if min(Fspec) lt max(Bspec) then begin
  B_bad++
  printf,log,' blue wing too high!'
endif
if min(Fspec) lt max(Rspec) then begin
  R_bad++
  printf,log,' red wing too high!'
endif
  
return,[B_bad,R_bad]
end

function coord2pix,hdr,value,dim,fitspix=fitspix
value=float(value)
dim=strcompress(string(dim),/remove_all)
if dim eq '1' && long(sxpar(hdr,'CRPIX1')) le 0l then pix=(value-360.-sxpar(hdr,'CRVAL'+dim))/sxpar(hdr,'CDELT'+dim)+sxpar(hdr,'CRPIX'+dim) $
 else pix=round((value-sxpar(hdr,'CRVAL'+dim))/sxpar(hdr,'CDELT'+dim)+sxpar(hdr,'CRPIX'+dim),/l64)
if keyword_set(fitspix) then return,round(pix,/l64) else return,round(pix,/l64)-1
end

function pix2coord,hdr,pixel,dim,fitspix=fitspix
if keyword_set(fitspix) then pixel=round(pixel,/l64) else pixel=round(pixel,/l64)+1l
;if keyword_set(fitspix) then pixel=long(pixel) else pixel=long(pixel)+1l
dim=strcompress(string(dim),/remove_all)
return,(pixel-sxpar(hdr,'CRPIX'+dim))*sxpar(hdr,'CDELT'+dim)+sxpar(hdr,'CRVAL'+dim)
end

function iirms,hdr,sigma
return,sigma*sqrt(sxpar(hdr,'NAXIS3')*(sxpar(hdr,'CDELT3')/1000.)^2)
end

function co_scl,dat,min=min,max=max,sm=sm,flat=flat
if keyword_set(sm) then dat=smooth(dat,sm,/EDGE_MIRROR,/NAN)
if keyword_set(flat) then dat[where(dat le flat,/null)]=0
if ~keyword_set(min) || ~keyword_set(max) then begin
  min=0
  max=mean(dat)
endif
;min=min
;max=max
;dat=bytscl(dat,min=mini,max=maxi)
redat=reform(dat,1,n_elements(dat[*,0]),n_elements(dat[0,*]))
return,bytscl(redat,min=min,max=max)
end

function array_band, array, band=band
if ~keyword_set(band) then band=0 else band=0.5*(array[1]-array[0])
sz=size(array)
band_arr=make_array(2*sz[1])
for i=0,sz[1]-1 do begin
  band_arr[2*i:2*i+1]=[array[i]-band,array[i]+band]
endfor
return,band_arr
end

function mkpvdat, dat, hdr, a, d, gal=gal, step=step
;fits_read,fitsfile,dat,hdr
if keyword_set(gal) then gal = 1b else gal = 0b
fitsgal = strcmp(sxpar(hdr,'CTYPE1'), 'GL', 2, /fold_case)
if fitsgal eq gal then at=a & dt=d
if fitsgal and ~gal then glactc,a,d,2000,at,dt,1,/deg
if ~fitsgal and gal then glactc,at,dt,2000,a,d,2,/deg
sxaddpar,hdr,'CTYPE1',repstr(sxpar(hdr,'CTYPE1'),'GLS','SFL')
sxaddpar,hdr,'CTYPE2',repstr(sxpar(hdr,'CTYPE2'),'GLS','SFL')
adxy,hdr,at,dt,x,y
defstep=0.5 ;unit pixel
if ~keyword_set(step) then step = defstep
if step le 0 then step = defstep
length = sqrt( ((x-shift(x,1))[1:*])^2+((y-shift(y,1))[1:*])^2 )
nstep = ceil(total(length)/step)
step = total(length)/nstep
;print,'Use step '+string(step)+' * pixelsize'
xs = x[0]+(x[1]-x[0])*findgen(nstep+1)/nstep
ys = y[0]+(y[1]-y[0])*findgen(nstep+1)/nstep
intlen = length
for i=n_elements(length)-1,0,-1 do intlen[i]=total(length[0:i])
for i=0, nstep do begin
  way = i*step
  node = (where(way le intlen))[0]
  xs[i] = x[node+1]-(intlen[node]-way)/length[node]*(x[node+1]-x[node])
  ys[i] = y[node+1]-(intlen[node]-way)/length[node]*(y[node+1]-y[node])
endfor
nx1 = sxpar(hdr,'NAXIS3')
nx2 = n_elements(xs)
slice = make_array(nx1, nx2, type=size(dat,/type))
for i=0,nx1-1 do slice[i,*] = interpolate(dat[*,*,i],xs,ys,missing=0)
mkhdr,pvhdr,slice
sxaddpar,pvhdr,'CTYPE1','VELOCITY'
sxaddpar,pvhdr,'CRPIX1',sxpar(hdr,'CRPIX3')
sxaddpar,pvhdr,'CRVAL1',sxpar(hdr,'CRVAL3')/1000d
sxaddpar,pvhdr,'CDELT1',sxpar(hdr,'CDELT3')/1000d
sxaddpar,pvhdr,'CTYPE2','POSITION'
sxaddpar,pvhdr,'CRPIX2',1
sxaddpar,pvhdr,'CRVAL2',0d
sxaddpar,pvhdr,'CDELT2',step*abs(sxpar(hdr,'CDELT1'))
;sxaddhist,'PV file: '+fitsfile,pvhdr
sxaddhist,'PV path:',pvhdr
for i=0,n_elements(a)-1 do sxaddhist,string(a[i])+' '+string(d[i]),pvhdr
sxaddhist,'Position in Degree',pvhdr
sxaddhist,'Velocity in km/s',pvhdr
outdat={dat:slice,hdr:pvhdr}
;fits_write,'pvslice.fits',slice,pvhdr
return, outdat
end

PRO IDENT_OUTFLOW_TEST

;cd,'/home/data/Galaxy/-18_30/TEST/'
cd,'/home/Alpha/Astrodata/galaxysurvey/-18_30/'

;RESOLVE_ROUTINE,'auto_wing',/is_function;,/COMPILE_FULL_FILE
;RESOLVE_ROUTINE,'multi_check',/is_function;,/COMPILE_FULL_FILE

;find 13CO peaks.
region='GMCC1'
lrange=[192.5,195.25]
brange=[-2,0.6]
Uvrange=[-3,10]
Lvrange=[-0.3,10]
L2vrange=[2,9]

if ~file_test('./'+region) then spawn,'mkdir ./'+region
cd,'./'+region

outnameU=region+'_U'
outnameL=region+'_L'
outnameL2=region+'_L2'

if ~file_test(region+'_U_C.fits') then cropfits,'../mosaic_U.fits',lrange,brange,outfile=region+'_U',Uvrange
if ~file_test(region+'_L_C.fits') then cropfits,'../mosaic_L.fits',lrange,brange,outfile=region+'_L',Lvrange
if ~file_test(region+'_L2_C.fits') then cropfits,'../mosaic_L2.fits',lrange,brange,outfile=region+'_L2',L2vrange

if ~file_test(region+'_Ua_C.fits') then cropfits,'../mosaic_U.fits',lrange,brange,outfile=region+'_Ua'
if ~file_test(region+'_La_C.fits') then cropfits,'../mosaic_L.fits',lrange,brange,outfile=region+'_La'
if ~file_test(region+'_L2a_C.fits') then cropfits,'../mosaic_L2.fits',lrange,brange,outfile=region+'_L2a'

conf13='fellwalker.conf'
openw,c13,/get_lun,conf13
printf,c13,'FellWalker.AllowEdge=1'
printf,c13,'FellWalker.CleanIter=1'
printf,c13,'FellWalker.FlatSlope=2*RMS'
printf,c13,'FellWalker.FwhmBeam=2'
printf,c13,'FellWalker.MaxBad=0.2'
printf,c13,'FellWalker.MaxJump=7'
printf,c13,'FellWalker.MinHeight=3';10*RMS'
printf,c13,'FellWalker.MinPix=16'
printf,c13,'FellWalker.MinDip=0*RMS';.5*RMS'
printf,c13,'FellWalker.Noise=1*RMS'
printf,c13,'FellWalker.VeloRes=0'
printf,c13,'FellWalker.RMS=0.3'
free_lun,c13

;convert fits to ndf; starlink-2015A
L_file_fits=region+'_La_C.fits'
L_file_sdf=region+'_L.sdf'
L_file_mask=region+'_L_mask.sdf'
L_clumps_cat=region+'_L_clumps.fit'
cmdfile=region+'_cmd'
openw,cmd,/get_lun,cmdfile
;printf,cmd,'export STARLINK_DIR=/usr/local/astrosoft/star-2015A'
printf,cmd,'export STARLINK_DIR=/usr/local/Astrosoft/star-2015A'
printf,cmd,'source $STARLINK_DIR/etc/profile'
printf,cmd,'convert'
printf,cmd,'cupid'
convert_str='fits2ndf '+L_file_fits+' '+L_file_sdf
findclumps_str="findclumps in="+L_file_fits+" deconv=no method=fellwalker out="+L_file_mask+" outcat="+L_clumps_cat+" config=""'^"+conf13+"'"" wcspar=yes"+" SHAPE=""Ellipse"""
printf,cmd,findclumps_str
free_lun,cmd
spawn,'chmod +x '+cmdfile
spawn,'./'+cmdfile

;read catalog in wcs units pix*cdelt
;peak[1|2] cen[1|2]: galactic degree
;peak3 cen3: velocity m/s
;size[1|2]: arcsec
;size3: m/s
ftab_ext,L_clumps_cat,[1,2,3,4,5,6,7,8,9,10,12,14],$
  idx,peak1,peak2,peak3,cen1,cen2,cen3,size1,size2,size3,peak,shape
openw,cat,/get_lun,'L_clumps.cat'
openw,elli,/get_lun,'L_clumps_shape.txt'
printf,cat,'#num','peak1','peak2','peak3','cen1','cen2','cen3','size1','size2','size3','Tpeak',$
  format='(a4,2(2a9,a8),2a6,2a7)'
printf,cat,'#   ','[deg]','[deg]','[km/s]','[deg]','[deg]','[km/s]','['']','['']','[km/s]','[K]',$
  format='(a4,2(2a9,a8),2a6,2a7)'

count=0
for i=0, n_elements(idx)-1 do begin
  if peak3[i]/1000. le Lvrange[1] && size3[i] gt 300 then begin
    count++
    printf,cat,count,peak1[i],peak2[i],peak3[i]/1000.,$
    cen1[i],cen2[i],cen3[i]/1000.,$
    size1[i]/60.,size2[i]/60,size3[i]/1000.,peak[i],$
    format='(i4,2(2f9.3,f8.3),2f6.1,f7.1,f7.2)'
    printf,elli,count,shape[i],format='(i4,x,a0)'
  endif
endfor
free_lun,cat,elli

delvar,idx,peak1,peak2,peak3,cen1,cen2,cen3,size1,size2,size3,peak,shape

readcol,'L_clumps.cat',idx,peak1,peak2,peak3,cen1,cen2,cen3,size1,size2,size3,peak,format='I,F,F,F,F,F,F,F,F,F,F',stringskip='#',/silent
readcol,'L_clumps_shape.txt',num,str1,str2,str3,epos1,epos2,major,minor,posangle,format='I,A,A,A,F,F,F,F,F',/silent
epos1[where(epos1 lt 0)]+=360.
openw,scat,/get_lun,'L_clumps_shape.cat'
printf,scat,'#num','cen1','cen2','major','minor','posang',format='(a4,a9,a7,2a6,a7)'
;printf,scat,'#   ','[deg]'
for i=0, n_elements(num)-1 do $
  printf,scat,num[i],epos1[i],epos2[i],major[i],minor[i],posangle[i],$
  format='(i4,f9.3,f7.3,2f6.3,f7.1)'
free_lun,scat

;draw 13CO clumps on false RGB color image
fits_read,outnameU+'_C.fits',datU,hdrU
fits_read,outnameL+'_C.fits',datL,hdrL
fits_read,outnameL2+'_C.fits',datL2,hdrL2
;print,sxpar(hdrU,'NAXIS3')
iirmsU=iirms(hdrU,0.5)
iirmsL=iirms(hdrL,0.3)
iirmsL2=iirms(hdrL2,0.3)
;print,iirmsU,iirmsL,iirmsL2
iidatU=total(datU,3)*sxpar(hdrU,'CDELT3')/1000.
iidatL=total(datL,3)*sxpar(hdrL,'CDELT3')/1000.
iidatL2=total(datL2,3)*sxpar(hdrL2,'CDELT3')/1000.
;print,max(iidatU),max(iidatL),max(iidatL2)
;print,mean(iidatU),mean(iidatL),mean(iidatL2)
imR=co_scl(iidatL2,min=0.3*2,max=3,sm=2)
imG=co_scl(iidatL,min=0.4*2,max=20,sm=1)
imB=co_scl(iidatU,min=0.7*2,max=70)
imRGB=[imR,imG,imB]

crp1=sxpar(hdrU,'CRPIX1')
crv1=sxpar(hdrU,'CRVAL1')
del1=sxpar(hdrU,'CDELT1')
l_l=(360.+(0.5-crp1)*del1+crv1) mod 360
l_r=(360.+(sxpar(hdrU,'NAXIS1')+0.5-crp1)*del1+crv1) mod 360
   
crp2=sxpar(hdrU,'CRPIX2')
crv2=sxpar(hdrU,'CRVAL2')
del2=sxpar(hdrU,'CDELT2')
b_d=(0.5-crp2)*del2+crv2
b_u=(sxpar(hdrU,'NAXIS2')+0.5-crp2)*del2+crv2

x_range=[l_l,l_r]
y_range=[b_d,b_u]
ratio=abs((b_u-b_d)/(l_l-l_r))
psname=region+'.eps'
cgps_open,psname,font=!p.font,/quiet,default_thickness=1.0,charsize=1.0;,/encapsulated;,/portrait
xsize=1000. & ysize=1000.*ratio
cgDisplay, xsize=round(xsize), ysize=round(ysize)
pos0=[100./xsize,100.*ratio/ysize,1.-50./xsize,1-50./ysize]
cgplot,[0],[0],xrange=x_range,yrange=y_range,$
  xtickinterval=0.5,ytickinterval=0.5,aspect=ratio,AxisColor='black',xthick=5,ythick=5,$
  ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})'),position=pos0
cgimage,imRGB,/overplot;,/noerase,position=pos0

cgplot,[0],[0],xrange=x_range,yrange=y_range,position=pos0,$
  xtickinterval=0.5,ytickinterval=0.5,aspect=ratio,$
  AxisColor='white',/noerase,xtickformat='(a1)',ytickformat='(a1)';,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')
;locate core peak

;cgsymcat

for i=0, n_elements(num)-1 do begin
  cgPlotS,draw_ellipse(epos1[i],epos2[i],major[i],minor[i],pa=-90-posangle[i]), Color='green',NOCLIP=0,thick=2
  cgplotS,draw_ellipse(epos1[i],epos2[i],major[i],minor[i],pa=-90-posangle[i]), Color='white',NOCLIP=0,thick=1
  cgplot,/overplot,peak1[i],peak2[i],psym=7,symsize=0.5,symcolor='green'
  cgtext,epos1[i],epos2[i],color='red',/data,num2str(num[i]),charsize=0.5,alignment=0.5,charthick=2
  cgtext,epos1[i],epos2[i],color='white',/data,num2str(num[i]),charsize=0.5,alignment=0.5,charthick=1
endfor

cgps_close
  
;extract 12CO and 13CO spectrum at 13CO peak

fits_read,region+'_Ua_C.fits',datUa,hdrUa 
fits_read,region+'_La_C.fits',datLa,hdrLa
fits_read,region+'_L2a_C.fits',datL2a,hdrL2a

;timestr=STRJOIN(STRSPLIT(systime(), /EXTRACT), '_')
openw,log,/get_lun,'log.txt'
openw,out,/get_lun,'out_12CO.cat'

printf,out,"# AY[AN] = Auto Yes[NO]
printf,out,"# B[R,D]u = blue[red,double] uncertain"
printf,out,"# B[R,D]m = blue[red,double] manucheck"
printf,out,'# AY/AN',' out_c','B/R/D','peak_n','peak_l','peak_b','peak_v','bluewing','redwing',$
  format='(4a7,a8,2a7,2a12)'
printf,out,'#      ','','','','','','','L','R','L','R',format='(4a7,a8,2a7,4a6)'
printf,out,'#      ','','','','[deg]','[deg]','[km/s]','[km/s]','[km/s]','[km/s]','[km/s]',$
  format='(4a7,a8,2a7,4a6)'
  
ocyuc=0; Outflow Candidates Yes Uncertain Count
ocnuc=0; Outflow Candidates No Uncertain Count

for i=0, n_elements(num)-1 do begin
  pixsz1=round(size1[i]/120.)>1
  pixsz1=pixsz1<2
  pixsz2=round(size2[i]/120.)>1
  pixsz2=pixsz2<2
  ;;;;;;;;;;;;;;;;;
  ;pixsz1=1 & pixsz2=1
  ;;;;;;;;;;;;;;;;;
  szMAXU=size(datUa)
  szMAXL=size(datLa)
  szMAXL2=size(datL2a)
  specU=datUa[(coord2pix(hdrUa,peak1[i],1)-pixsz1)>0:(coord2pix(hdrUa,peak1[i],1)+pixsz1)<(szMAXU[1]-1),$
    (coord2pix(hdrUa,peak2[i],2)-pixsz2)>0:(coord2pix(hdrUa,peak2[i],2)+pixsz2)<(szMAXU[2]-1),*]
  specL=datLa[(coord2pix(hdrLa,peak1[i],1)-pixsz1)>0:(coord2pix(hdrLa,peak1[i],1)+pixsz1)<(szMAXL[1]-1),$
    (coord2pix(hdrLa,peak2[i],2)-pixsz2)>0:(coord2pix(hdrLa,peak2[i],2)+pixsz2)<(szMAXL[2]-1),*]
  specL2=datL2a[(coord2pix(hdrL2a,peak1[i],1)-pixsz1)>0:(coord2pix(hdrL2a,peak1[i],1)+pixsz1)<(szMAXL2[1]-1),$
    (coord2pix(hdrL2a,peak2[i],2)-pixsz2)>0:(coord2pix(hdrL2a,peak2[i],2)+pixsz2)<(szMAXL2[2]-1),*]
  specU=mean(specU,dim=1)
  specU=mean(specU,dim=1)
  specL=mean(specL,dim=1)
  specL=mean(specL,dim=1)
  specL2=mean(specL2,dim=1)
  specL2=mean(specL2,dim=1)
  sw=3
  specU=smooth(specU,sw,/EDGE_MIRROR)
  specL=smooth(specL,sw,/EDGE_MIRROR)
  specL2=smooth(specL2,sw,/EDGE_MIRROR)  
  vU=((indgen(n_elements(specU))+1-double(sxpar(hdrUa,'CRPIX3')))*double(sxpar(hdrUa,'CDELT3'))+double(sxpar(hdrUa,'CRVAL3')))/1000.0
  vL=((indgen(n_elements(specL))+1-double(sxpar(hdrLa,'CRPIX3')))*double(sxpar(hdrLa,'CDELT3'))+double(sxpar(hdrLa,'CRVAL3')))/1000.0
  vL2=((indgen(n_elements(specL2))+1-double(sxpar(hdrL2a,'CRPIX3')))*double(sxpar(hdrL2a,'CDELT3'))+double(sxpar(hdrL2a,'CRVAL3')))/1000.0
  
  peakposU=where(vU ge peak3[i]-0.085 and vU le peak3[i]+0.085,/null)
  peakposL=where(vL ge peak3[i]-0.085 and vL le peak3[i]+0.085,/null) 
  printf,log, 'Candidate'+num2str(i+1)
  FWHMU=vU[auto_wing(specU,peakposU[0],specU[peakposU[0]]/1.7)]
  width=FWHMU[1]-FWHMU[0]
  vrange=[peak3[i]-3.5*width,peak3[i]+3.5*width] 
  rmsU=0.5/sqrt(sw*(2.*pixsz1+1.)*(2.*pixsz2+1.))
  rmsL=(0.06*specL[peakposL[0]])>(0.3/sqrt((2.*pixsz1+1.)*(2.*pixsz2+1.))) 
  wingU=vU[auto_wing(specU,peakposU[0],rmsU)]
  wingL=vL[auto_wing(specL,peakposL[0],rmsL)]
  baseline=make_array(n_elements(specU),value=0)
  
  check=multi_check(peak3[i],specU,vU,FWHMU,wingU,wingL,rmsU,rmsL,log,pixsz1)
  
  specpsname=region+'_'+num2str(idx[i])+'.eps'
  cgps_open,specpsname,font=!p.font,/quiet,default_thickness=1.0,charsize=1.0;,/encapsulated;,/portrait
  cgplot,vU,specU,psym=10,/nodata,xrange=vrange,color='blue',yrange=[0-max(specU)*0.15,max(specU)*1.15],$
    xtickinterval=5,xtitle='LSR Velocity (km s!U-1!N)',ytitle='T!DMB!N (K)',position=p0,$
    title=region+' CLUMP '+num2str(idx[i])+': '+num2str(peak1[i],format='(f7.3)')+num2str(peak2[i],format='(f+6.3)')
  
  if ~check[0] then begin
    bpolyrg=where(vU ge wingU[0] and vU le wingL[0])
    bx1_arr=array_band(vU[bpolyrg],/band)
    bx2_arr=array_band(Reverse(vU[bpolyrg]),/band)
    by1_arr=array_band(specU[bpolyrg])
    by2_arr=array_band(baseline[bpolyrg])
    cgColorFill,[bx1_arr, bx2_arr, bx1_arr[0]],[by1_arr,by2_arr,by1_arr[0]],Color='Sky Blue'
  endif

  if ~check[1] then begin
    rpolyrg=where(vU ge wingL[1] and vU le wingU[1])
    rx1_arr=array_band(vU[rpolyrg],/band)
    rx2_arr=array_band(Reverse(vU[rpolyrg]),/band)
    ry1_arr=array_band(specU[rpolyrg])
    ry2_arr=array_band(baseline[rpolyrg])
    cgColorFill,[rx1_arr, rx2_arr, rx1_arr[0]],[ry1_arr,ry2_arr,ry1_arr[0]],Color='pink'
  endif
  
  if check[0] && check[1] then begin
    ocnuc++
    printf,out,'#    AN',string(ocnuc,format="(i0,'*')"),'---',$
      i+1,peak1[i],peak2[i],peak3[i],wingU[0],wingL[0],wingL[1],wingU[1],$
      format='(3a7,i7,f8.3,f7.3,f7.2,4f6.1)'
  endif 
  if ~check[0] && check[1] then begin
    ocyuc++
    printf,out,'AY',ocyuc,'Bu',i+1,peak1[i],peak2[i],peak3[i],wingU[0],wingL[0],wingL[1],wingU[1],$
      format='(2(a7,i7),f8.3,f7.3,f7.2,4f6.1)'
  endif
  if check[0] && ~check[1] then begin
    ocyuc++
    printf,out,'AY',ocyuc,'Ru',i+1,peak1[i],peak2[i],peak3[i],wingU[0],wingL[0],wingL[1],wingU[1],$
      format='(2(a7,i7),f8.3,f7.3,f7.2,4f6.1)'
  endif
  if ~check[0] && ~check[1] then begin
    ocyuc++
    printf,out,'AY',ocyuc,'Du',i+1,peak1[i],peak2[i],peak3[i],wingU[0],wingL[0],wingL[1],wingU[1],$
      format='(2(a7,i7),f8.3,f7.3,f7.2,4f6.1)'
  endif
  
  cgplot,vU,specU,psym=10,color='blue',/overplot
  cgplot,vL,specL,psym=10,color='green',/overplot
  cgplot,vL2,specL2,psym=10,color='red',/overplot
  ;AL_Legend, items, LineStyle=lines, Color=colors, box=0, position=[0.575*(!x.CRange[1]-!x.Crange[0])+!x.Crange[0],0.95*(!y.CRange[1]-!y.Crange[0])+!y.Crange[0]]
  cgPlot, !x.CRange,[0.,0.], LineStyle=0, Color='black',/overplot
  cgPlot, wingU,[rmsU,rmsU], LineStyle=1, Color='blue',/overplot,thick=5
  cgPlot, wingL,[rmsL,rmsL], LineStyle=1, Color='red',/overplot,thick=5
  cgPlot, [peak3[i],peak3[i]],[0.1,!y.CRange[1]], LineStyle=1,Color='black',/overplot
  cgtext, peak3[i],0-max(specU)*0.075,num2str(peak3[i],format='(f7.1)'),alignment=0.5
  cgps_close
endfor

free_lun,log,out

;Draw outflow contours
;'out_12CO.catm'
readcol,'out_12CO.cat',autosign,out_c,lobesign,peak_n,peak_l,peak_b,peak_v,blue_l,blue_r,red_l,red_r,$
  format='A,I,A,I,F,F,F,F,F,F,F',stringskip='#',/silent

psname=region+'_out_12CO.eps'
cgps_open,psname,font=!p.font,/quiet,default_thickness=1.0,charsize=1.0;,/encapsulated;,/portrait
xsize=1000. & ysize=1000.*ratio
cgDisplay, xsize=round(xsize), ysize=round(ysize)
pos0=[100./xsize,100.*ratio/ysize,1.-50./xsize,1-50./ysize]
cgplot,[0],[0],xrange=x_range,yrange=y_range,$
  xtickinterval=0.5,ytickinterval=0.5,aspect=ratio,AxisColor='black',xthick=5,ythick=5,$
  ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})'),position=pos0
cgimage,imRGB,/overplot;,/noerase,position=pos0
cgplot,[0],[0],xrange=x_range,yrange=y_range,position=pos0,$
  xtickinterval=0.5,ytickinterval=0.5,aspect=ratio,$
  AxisColor='white',/noerase,xtickformat='(a1)',ytickformat='(a1)'
for i=0, n_elements(out_c)-1 do begin
  ;cgPlotS,draw_ellipse(epos1[i],epos2[i],major[i],minor[i],pa=posangle[i]), Color='green',NOCLIP=0,thick=2
  ;cgplotS,draw_ellipse(epos1[i],epos2[i],major[i],minor[i],pa=posangle[i]), Color='white',NOCLIP=0,thick=1
  cgplot,/over,peak_l[i],peak_b[i],psym=6,symsize=0.5,symcolor='dark green',thick=3
  cgplot,/over,peak_l[i],peak_b[i],psym=6,symsize=0.5,symcolor='white',thick=1
  ;cgplot,/over,peak_l[i],peak_b[i],psym=1,symsize=0.5,symcolor='green',thick=3
  ;cgplot,/over,peak_l[i],peak_b[i],psym=1,symsize=0.5,symcolor='white',thick=1
  cgtext,peak_l[i]-0.005,peak_b[i]+0.01,color='red',/data,'c'+num2str(num[i]),charsize=0.5,alignment=0,charthick=3
  cgtext,peak_l[i]-0.005,peak_b[i]+0.01,color='white',/data,'c'+num2str(num[i]),charsize=0.5,alignment=0,charthick=1
endfor
cgps_close

;cgsymcat
;cgcolor
;adxy
;cropfits

  ;fits_read,region+'_Ua_C.fits',datUa,hdrUa 
  ;fits_read,region+'_La_C.fits',datLa,hdrLa
  ;fits_read,region+'_L2a_C.fits',datL2a,hdrL2a
;readcol,'out_12CO.cat',autosign,out_c,lobesign,peak_n,peak_l,peak_b,peak_v,blue_l,blue_r,red_l,red_r,$
  ;format='A,I,A,I,F,F,F,F,F,F,F',stringskip='#',/silent
;    biicon0=total(blueimage[*,*,*],3)*abs(sxpar(hdrU,'CDELT3'))/1000.0
;    biisz=size(biicon0)
;    ;biicon=temporary(congrid(biicon0,5*biisz[1],5*biisz[2],cubic=-0.5,/MINUS_ONE))
;    biicon=temporary(congrid(biicon0,4*biisz[1],4*biisz[2],cubic=-0.3,/minus_one))
;    bmax=max(biicon0[fix(biisz[1]/3):fix(2*biisz[1]/3),fix(biisz[2]/3):fix(2*biisz[2]/3)])
;    levelsign=where(sign3 eq sign1[i])
;    levelsign=long(mean(levelsign))
;    if levelsign ne -1 then begin
;      dbcon=bmax*dblev[levelsign]/10.0
;      bcon0=bmax*blev0[levelsign]
;    endif else begin
;      dbcon=bmax*1.0/10.0
;      bcon0=bmax*0.4
;    endelse
;    levels1=indgen(10)*dbcon+bcon0

;psname=region+'_out_12CO.eps'
;cgps_open,psname,font=!p.font,/quiet,default_thickness=1.0,charsize=1.0;,/encapsulated;,/portrait
;xsize=1000. & ysize=1000.*ratio
;cgDisplay, xsize=round(xsize), ysize=round(ysize)
;pos0=[100./xsize,100.*ratio/ysize,1.-50./xsize,1-50./ysize]
;cgplot,[0],[0],xrange=x_range,yrange=y_range,$
;  xtickinterval=0.5,ytickinterval=0.5,aspect=ratio,AxisColor='black',xthick=5,ythick=5,$
;  ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})'),position=pos0
;cgimage,imRGB,/noerase,position=pos0
;cgplot,[0],[0],xrange=x_range,yrange=y_range,position=pos0,$
;  xtickinterval=0.5,ytickinterval=0.5,aspect=ratio,$
;  AxisColor='white',/noerase,xtickformat='(a1)',ytickformat='(a1)'
;for i=0, n_elements(out_c)-1 do begin
;  ;cgPlotS,draw_ellipse(epos1[i],epos2[i],major[i],minor[i],pa=posangle[i]), Color='green',NOCLIP=0,thick=2
;  ;cgplotS,draw_ellipse(epos1[i],epos2[i],major[i],minor[i],pa=posangle[i]), Color='white',NOCLIP=0,thick=1
;  cgplot,/over,peak_l[i],peak_b[i],psym=6,symsize=0.6,symcolor='black',thick=1
;  cgplot,/over,peak_l[i],peak_b[i],psym=1,symsize=0.5,symcolor='green',thick=3
;  cgplot,/over,peak_l[i],peak_b[i],psym=1,symsize=0.5,symcolor='white',thick=1
;  cgtext,peak_l[i]-0.005,peak_b[i]+0.01,color='red',/data,num2str(num[i]),charsize=0.5,alignment=0,charthick=3
;  cgtext,peak_l[i]-0.005,peak_b[i]+0.01,color='white',/data,num2str(num[i]),charsize=0.5,alignment=0,charthick=1
;endfor
;cgps_close


if ~file_test('figures') then spawn,'mkdir figures'
spawn,'rm ./figures/*.eps'
spawn,'mv *.eps figures'


END