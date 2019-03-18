function coord2pix,hdr,value,dim,fitspix=fitspix
value=float(value)
dim=strcompress(string(dim),/remove_all)
if dim eq '1' && long(sxpar(hdr,'CRPIX1')) lt 0l then pix=(value-360.-sxpar(hdr,'CRVAL'+dim))/sxpar(hdr,'CDELT'+dim)+sxpar(hdr,'CRPIX'+dim) $
 else pix=round((value-sxpar(hdr,'CRVAL'+dim))/sxpar(hdr,'CDELT'+dim)+sxpar(hdr,'CRPIX'+dim),/l64)
if keyword_set(fitspix) then return,round(pix,/l64) else return,round(pix,/l64)-1
end

function pix2coord,hdr,pixel,dim,fitspix=fitspix
if keyword_set(fitspix) then pixel=round(pixel,/l64) else pixel=round(pixel,/l64)+1l
;if keyword_set(fitspix) then pixel=long(pixel) else pixel=long(pixel)+1l
dim=strcompress(string(dim),/remove_all)
return,(pixel-sxpar(hdr,'CRPIX'+dim))*sxpar(hdr,'CDELT'+dim)+sxpar(hdr,'CRVAL'+dim)
end

function getfitsname,gl,gb
gl = gl mod 360
return,string(round(gl*2)*5,format='(I04)')+string(round(gb*2)*5,format='(I+004)')
end

function boxsmooth,array,width,edge_mirror=edge_mirror
sz=size(array)
width=float(round(width>1))
sz_mod=sz[1] mod round(width)
if keyword_set(edge_mirror) then n_sm=ceil(sz[1]/width) else n_sm=ceil((sz[1]-sz_mod)/width)
sm_arr=make_array(n_sm,/double)
for i=0,sz[1]-sz_mod-1,round(width) do begin
  sm_arr[round(i/width)]=total(array[i:i+round(width)-1])/width
endfor
if keyword_set(edge_mirror) && (sz_mod gt 0) then begin
  sm_arr[n_elements(sm_arr)-1]=(total(array[sz[1]-sz_mod-1:sz[1]-1])+total(array[sz[1]-1-(round(width)-sz_mod):sz[1]-1]))/width
endif
return,sm_arr
end

function cubefilter,dat,threshold
iimap=total(dat[*,*,*],3)
pos=where(iimap lt threshold)
if max(pos) ne -1 then begin
  iimap[pos]=0
  cov=reform(iimap) & cov[*]=1
  cov[where(iimap eq 0)]=0
  cov=total(cov)
  for k=0,n_elements(dat[0,0,*])-1 do begin
    temp=dat[*,*,k]
    temp[where(iimap eq 0)]=0
    dat[*,*,k]=temp
  endfor
endif else begin
  cov=reform(iimap) & cov[*]=1
  cov=total(cov)
endelse
strdata={dat:dat,cov:cov}
return,strdata
end

function maxcoord,dat,hdr
iimap=total(dat[*,*,*],3)
pos=where(iimap eq max(iimap))
indices=array_indices(iimap,pos[0])
coord=[pix2coord(hdr,indices[0],1),pix2coord(hdr,indices[1],2)]
return,coord
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

function ab_dev, arr, value
for i=0, n_elements(arr)-1 do arr_v=abs(arr[i]-value)
return,mean(arr_v)
end

pro mkoutflow;,vres=vres;,area=area

;if ~keyword_set(vres) then print, ' Use default velocity resolution.' 

;cd, '/media/alpha/TOSHIBA EXT'
;cd, '/home/data/CO OUTFLOW'
;cd,'/home/alpha/Astrodata/NH3_100m_point/OUTFLOW'
;cd,'/home/Alpha/Astrodata/OUTFLOW'
cd,'/home/lee/W3'
;cd,'/home/alpha/Astrodata/NH3_100m_point/REDUCTION/OUTFLOW'


; 1) read outflow parameters
;readcol,'outflowpara_test',sign1,l,b,v,v_width,v_bin,v1,v2,d,LINE,bv1,bv2,rv1,rv2,area_x,area_y,spec_x,spec_y,/silent,stringskip='#',$
;  format='(A,F,F,F,F,F,F,F,F,A,F,F,F,F,F,F,F,F)'
readcol,'outflowpara',sign1,l,b,v,v_width,v_bin,v1,v2,d,LINE,bv1,bv2,rv1,rv2,area_x,area_y,spec_x,spec_y,/silent,stringskip='#',$
  format='(A,F,F,F,F,F,F,F,F,A,F,F,F,F,F,F,F,F)'

readcol,'outflowpara2',sign2,blue_l,blue_b,red_l,red_b,bszx,bszy,rszx,rszy,/silent,stringskip='#',$
  format='(A,F,F,F,F,F,F,F,F)'
readcol,'outflowlevels',sign3,dblev,blev0,drlev,rlev0,/silent,stringskip='#',format='(A,F,F,F,F)'
openw,paraout,/get_lun,'outflow_derived.txt'
openw,paracat,/get_lun,'outflow_catlog.txt'
;blue_area,bI_12,bNH2,bMH2,bPH2,bEH2,format='(a21,f8.2,f10.2,e10.2,f9.2,f13.2,e10.2)'
printf,paraout,'sign','source_lobe','area','I(CO)','N(H2)','M(H2)','P(H2)','E(H2)', 'scale','mean_v',  format='(a9,a21,a8,a10,a10,a9,a13,a10,a8,a8)'
printf,paraout,'[pc^2]','[K km/s]','[cm^-2]','[M_sun]','[M_sun km/s]','[J]','[pc]','[km/s]',format='(9x,21x,a8,a10,a10,a9,a13,a10,a8,a8)'

printf,paracat,'sign','line','source','v_c','area_b','area_r','scale_b','scale_r','delt_vb','delt_vr','I_b','I_r','N(H2)_b','N(H2)_r','M(H2)_b','M(H2)_r','M(H2)_out',$
  'p(H2)','E(H2)','time','M(H2)_out/t','Fm','Lm','M(H2)_core_18','Tex_12','scale_m','scl_limit',$
    format='(a10,a6,a15,7a10,7a11,a13,a10,a10,a15,a20,a10,a15,3a10)'
printf,paracat,'[km/s]','[pc^2]','[pc^2]','[pc]','[pc]','[km/s]','[km/s]','[K km/s]','[K km/s]','[cm^-2]','[cm^-2]','[M_sun]','[M_sun]','[M_sun]',$
  '[M_sun km/s]','[ergs]','[year]','[M_sun yr^-1]','[M_sun km/s yr^-1]','[L_sun]','[M_sun]','[K]','[pc]','[pc]',$
    format='(31x,7a10,7a11,a13,a10,a10,a15,a20,a10,a15,3a10)'

openw,fitsunit,/get_lun,'invalid_fits.txt'
openw,log,/get_lun,'log.txt'

; 2) make outflows
for i=0, n_elements(sign1)-1 do begin
  printf,log,'#######################################'
  printf,log,sign1[i],strcompress(string(l[i],format='(f7.3)')+string(b[i],format='(f+6.3)'),/remove_all),systime(),format='(a,4x,a)'
  printf,log,' '
  
  ;if bv1[i] eq -999 || bv2[i] eq -999 then print, 'No blue wings.'
  ;if rv1[i] eq -999 || rv2[i] eq -999 then print, 'No red wings.'

; 2.1  
;;;;;calculate maximum size of outflow in arcmin scale with a value of 8;;;;;;;;;;;;;;

  szmax=10.0
  printf,log,szmax,szmax/60d,format='("maxsz = ",f6.3," arcmin, ",f6.4," degree.")'
  area=double([area_x[i],area_y[i]])
  area=double(area<szmax)
  area=area>0.5
  printf,log,area[0],area[0]/60d,format='("areax = ",f6.3," arcmin, ",f6.4," degree.")'
  printf,log,area[1],area[1]/60d,format='("areay = ",f6.3," arcmin, ",f6.4," degree.")'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  print,'Begin:'+strcompress(string(i+1,format='(i)'))
 
  source=strcompress(string(l[i],format='(f7.3)')+string(b[i],format='(f+6.3)'),/remove_all)
  sourcemini='G'+strcompress(string(round(l[i]*10)/10.0,format='(f5.1)')+string(round(b[i]*10)/10.0,format='(f+4.1)'),/remove_all)
  if source eq '45.467+0.053' then sourcemini='G45.5+O.O'
  
  fitspath='./fitsfiles/'
  DLHUname=fitspath+'DLH_'+sourcemini+'_U.fits'
  DLHLname=fitspath+'DLH_'+sourcemini+'_L.fits'
  DLHL2name=fitspath+'DLH_'+sourcemini+'_L2.fits
  
  if ~file_test(fitspath+'DLH_'+sourcemini+'*.fits') then begin
  print,'  DLH_'+sourcemini+'*.fits NOT detected.'
  print,'  Making fits...'
  path='/run/media/Alpha/TOSHIBA EXT/OUTFLOWDATA/'
  ;path='/media/alpha/TOSHIBA EXT/OUTFLOWDATA/'
  ;path='./'  
  fitsname=getfitsname(l[i],b[i])
  print,fitsname
  fitsnameu=getfitsname(l[i],b[i]+0.25)
  fitsnamed=getfitsname(l[i],b[i]-0.25)
  fitsnamel=getfitsname(l[i]+0.25,b[i])
  fitsnamelu=getfitsname(l[i]+0.25,b[i]+0.25)
  fitsnameld=getfitsname(l[i]+0.25,b[i]-0.25)
  fitsnamer=getfitsname(l[i]-0.25,b[i])
  fitsnameru=getfitsname(l[i]-0.25,b[i]+0.25)
  fitsnamerd=getfitsname(l[i]-0.25,b[i]-0.25) 
  fitsnames=[fitsnamel,fitsnamelu,fitsnameld,fitsname,fitsnameu,fitsnamed,fitsnamer,fitsnameru,fitsnamerd]
  fitsnames=fitsnames[uniq(fitsnames,sort(fitsnames))]
  print,' Following '+string(n_elements(fitsnames),format='(i1)')+'*3 fits files will be needed:'
  print,fitsnames,format='(a10,"U[L/L2].fits")'
  if min(where(~file_test(path+fitsnames+'*.fits'))) ne -1 then begin
    invalid_pos=where(~file_test(path+fitsnames+'*.fits'))
    invalidfits=fitsnames[invalid_pos]
    if ~file_test(path+fitsname+'*.fits') then begin
      print,' Sorry, source of '+source+' is not available, no relevant fits file can be accessed.'
      continue
    endif
    if n_elements(invalidfits) ne 0 then begin
    print,' Following fits do not exist!'
    print,invalidfits,format='(a10,"U[L/L2].fits")'
    foreach name, invalidfits do begin 
      printf,fitsunit,name,format='(a10,"U.fits")'
      printf,fitsunit,name,format='(a10,"L.fits")'
      printf,fitsunit,name,format='(a10,"L2.fits")'
    endforeach
    endif
  endif
  if min(where(file_test(path+fitsnames+'*.fits'))) ne -1 then begin
  valid_pos=where(file_test(path+fitsnames+'*.fits'))
  validfits=fitsnames[valid_pos]
  print,validfits,format='(" Fits files,",a9,"U[L/L2].fits will be read.")'
  endif
  foreach name, validfits do begin
    cuberms,path+name+'U.fits',[v[i]-150,v[i]+150],window=[v[i]-100,v[i]+100],/silent
    cuberms,path+name+'L.fits',[v[i]-150,v[i]+150],window=[v[i]-100,v[i]+100],/silent
    cuberms,path+name+'L2.fits',[v[i]-150,v[i]+150],window=[v[i]-100,v[i]+100],/silent
  endforeach
  mosaic,l[i]-szmax/120d,l[i]+szmax/120d,b[i]-szmax/120d,b[i]+szmax/120d,v[i]-50,v[i]+50,sb='U',/silent,display=0,fitspath=path;,rmspath='./'
  spawn,'mv mosaic_U.fits '+DLHUname
  mosaic,l[i]-szmax/120d,l[i]+szmax/120d,b[i]-szmax/120d,b[i]+szmax/120d,v[i]-50,v[i]+50,sb='L',/silent,display=0,fitspath=path;,rmspath='./'
  spawn,'mv mosaic_L.fits '+DLHLname
  mosaic,l[i]-szmax/120d,l[i]+szmax/120d,b[i]-szmax/120d,b[i]+szmax/120d,v[i]-50,v[i]+50,sb='L2',/silent,display=0,fitspath=path;,rmspath='./'
  spawn,'mv mosaic_L2.fits '+DLHL2name
  print,'  Making fits done!'
  print,'  DLH_'+sourcemini+'_U[L/L2].fits will be read.'
  endif else print,'  DLH_'+sourcemini+'[M]_U[L/L2][M].fits ready to be read.'
  
  if file_test(fitspath+'DLH_'+sourcemini+'_*M.fits') then begin
    DLHUname=fitspath+'DLH_'+sourcemini+'_UM.fits'
    DLHLname=fitspath+'DLH_'+sourcemini+'_LM.fits'
    DLHL2name=fitspath+'DLH_'+sourcemini+'_L2M.fits
    print,'  DLH_'+sourcemini+'_UM[LM/L2M].fits ready to be read.'
  endif
  
  if file_test(fitspath+'DLH_'+sourcemini+'M_*.fits') then begin
    DLHUname=fitspath+'DLH_'+sourcemini+'M_U.fits'
    DLHLname=fitspath+'DLH_'+sourcemini+'M_L.fits'
    DLHL2name=fitspath+'DLH_'+sourcemini+'M_L2.fits
    print,'  DLH_'+sourcemini+'M_U[L/L2].fits ready to be read.'
  endif
  
    
;;;;;plotting configure;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  items = ['!U12!NCO(1-0)', '!U13!NCO(1-0)', '!NC!U18!NO(1-0)']
  lines = [10,10,10]
  colors = ['blue','green','red']
  p1=[0.55,0.685,0.95,0.985]
  p4=[0.55,0.35,0.95,0.625]
  p0=[0.55,0.075,0.95,0.35]
  
  p2=[0.05,0.075,0.5,0.5]
  ;p2button=[0.1,0.05,0.5,0.1]
  p3=[0.05,0.56,0.5,0.985]
  psname='f_'+source+'.eps'
  print,'  Figure '+psname+' will be created.'
  
  cgps_open,psname,font=!p.font,/portrait,/quiet,default_thickness=1.0,charsize=0.6
  cgDisplay, xsize=2000, ysize=1600
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

;;;;;;cropping feature spectra;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if spec_x[i] eq -999 then fspecxsz=area[0]/480d else fspecxsz=spec_x[i]/2d
  if spec_y[i] eq -999 then fspecysz=area[1]/480d else fspecysz=spec_y[i]/2d
  print, '  Extracting feature spectrum from a cropped area of '+string(fspecxsz*120,format='(f4.2)')+'*'+string(fspecysz*120,format='(f4.2)')+' arcmin square...'

  cropfits,DLHUname,[l[i]-fspecxsz,l[i]+fspecxsz],[b[i]-fspecysz,b[i]+fspecysz],[v[i]-v_width[i],v[i]+v_width[i]],output='stempU'
  fits_read,'stempU_C.fits',sdatU,shdrU
  
  cropfits,DLHUname,[l[i]-fspecxsz,l[i]+fspecxsz],[b[i]-fspecysz,b[i]+fspecysz],[v1[i],v2[i]],output='ctempU'
  fits_read,'ctempU_C.fits',cdatU,chdrU
  Tmb12=max(cdatU)*0.5
  Tex12=5.532/(alog(1+5.532/(0.819+Tmb12)))
  
  cv3U=double(sxpar(shdrU,'CRVAL3'))
  cp3U=double(sxpar(shdrU,'CRPIX3'))
  dl3U=double(sxpar(shdrU,'CDELT3'))
  maxcoordsU=maxcoord(cdatU,chdrU)  
  maxcoordsU_x=round(coord2pix(shdrU,maxcoordsU[0],1))
  maxcoordsU_y=round(coord2pix(shdrU,maxcoordsU[1],2))
  ;print,maxcoordsU_x,maxcoordsU_y
  specU=reform(sdatU[maxcoordsU_x,maxcoordsU_y,*])
  
  ;print,maxcoordsU,l[i],b[i]
  
  ;maxpoint=max(total(sdatU,3))
  ;sU=cubefilter(sdatU,maxpoint)
  ;specU=total(sU.dat,1)
  ;specU=total(specU,1)/sU.cov
  vU=((indgen(n_elements(specU))+1-cp3U)*dl3U+cv3U)/1000.0
  
  vU=boxsmooth(vU,v_bin[i])
  specU=boxsmooth(specU,v_bin[i])
  
  vrange=[min(vU),max(vU)]
  
  printf,log,' Predefined velocity range:'
  printf,log,[v[i]-v_width[i],v[i]+v_width[i]]
  printf,log,' Real vrange:'
  printf,log,vrange
  
  cropfits,DLHLname,[l[i]-fspecxsz,l[i]+fspecxsz],[b[i]-fspecysz,b[i]+fspecysz],[v[i]-v_width[i],v[i]+v_width[i]],output='stempL'
  fits_read,'stempL_C.fits',sdatL,shdrL
  cropfits,DLHLname,[l[i]-fspecxsz,l[i]+fspecxsz],[b[i]-fspecysz,b[i]+fspecysz],[v1[i],v2[i]],output='ctempL'
  fits_read,'ctempL_C.fits',cdatL,chdrL
  cv3L=double(sxpar(shdrL,'CRVAL3'))
  cp3L=double(sxpar(shdrL,'CRPIX3'))
  dl3L=double(sxpar(shdrL,'CDELT3'))
  maxcoordsL=maxcoord(cdatL,chdrL)
 
  maxcoordsL_x=round(coord2pix(shdrL,maxcoordsL[0],1))
  maxcoordsL_y=round(coord2pix(shdrL,maxcoordsL[1],2))
  specL=reform(sdatL[maxcoordsL_x,maxcoordsL_y,*])
  
  if line[i] eq '13CO' then begin
    maxcoordsU_x=round(coord2pix(shdrU,maxcoordsL[0],1))
    maxcoordsU_y=round(coord2pix(shdrU,maxcoordsL[1],2))
    specU=reform(sdatU[maxcoordsU_x,maxcoordsU_y,*])
    vU=((indgen(n_elements(specU))+1-cp3U)*dl3U+cv3U)/1000.0
    vU=boxsmooth(vU,v_bin[i])
    specU=boxsmooth(specU,v_bin[i])
  endif
  
  ;maxpoint=max(total(sdatL,3))
  ;sL=cubefilter(sdatL,maxpoint)
  ;specL=total(sL.dat,1)
  ;specL=total(specL,1)/sL.cov
  vL=((indgen(n_elements(specL))+1-cp3L)*dl3L+cv3L)/1000.0
  
  vL=boxsmooth(vL,v_bin[i])
  ;wspecL=specL[v1[i],v2[i]]
  specL=boxsmooth(specL,v_bin[i])
  
  v13minpos=min(where(vL ge v1[i]))
  v13maxpos=max(where(vL le v2[i]))
  v13=vL[where(specL eq max(specL[v13minpos:v13maxpos]))]
  ;print,v13
  
  cropfits,DLHL2name,[l[i]-fspecxsz,l[i]+fspecxsz],[b[i]-fspecysz,b[i]+fspecysz],[v1[i],v2[i]],output='ctempL2'
  fits_read,'ctempL2_C.fits',cdatL2,chdrL2
  
  core18=cubefilter(cdatL2,max(total(cdatL2,3))*0.5)
  T_ex=30
  c_dl3=double(sxpar(chdrL2,'CDELT3'))
  I_18=total(core18.dat)*abs(c_dl3)/(1000.0*core18.cov)
  printf,log,' Integrated intensity of C18O:',I_18,' K km/s'
  dcore_area=((0.5/60.0)*(!pi/180)*1000.0*d[i])^2
  core_area=core18.cov*dcore_area
  NH2_18=2.24e14*5e6*I_18/(1-exp(-5.27/Tex12))
  MH2_18=1.0/6.3e19*NH2_18*core_area
  
  
  cropfits,DLHL2name,[l[i]-fspecxsz,l[i]+fspecxsz],[b[i]-fspecysz,b[i]+fspecysz],[v[i]-v_width[i],v[i]+v_width[i]],output='stempL2'
  fits_read,'stempL2_C.fits',sdatL2,shdrL2
  cv3L2=double(sxpar(shdrL2,'CRVAL3'))
  cp3L2=double(sxpar(shdrL2,'CRPIX3'))
  dl3L2=double(sxpar(shdrL2,'CDELT3'))

  specL2=reform(sdatL2[maxcoordsL_x,maxcoordsL_y,*])

  vL2=((indgen(n_elements(specL2))+1-cp3L2)*dl3L2+cv3L2)/1000.0
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  maxcoordsL2=maxcoord(cdatL2,chdrL2)
  printf,log,' Peak point of C18O image:'
  printf,log,maxcoordsL2
  maxcoordsL2_x=round(coord2pix(shdrL2,maxcoordsL2[0],1))
  maxcoordsL2_y=round(coord2pix(shdrL2,maxcoordsL2[1],2))
  ;print,maxcoordsU_x,maxcoordsU_y
  specL2s=reform(sdatL2[maxcoordsL2_x,maxcoordsL2_y,*])
  vL2s=((indgen(n_elements(specL2s))+1-cp3L2)*dl3L2+cv3L2)/1000.0
  
  vL2s=boxsmooth(vL2s,v_bin[i])
  specL2s=boxsmooth(specL2s,v_bin[i])
;  VU=BOXSMOOTH(VU,V_BIN[I])
;  SPECU=BOXSMOOTH(SPECU,V_BIN[I])
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  vL2=boxsmooth(vL2,v_bin[i])
  specL2=boxsmooth(specL2,v_bin[i])
  
  v18minpos=min(where(vL2 ge v1[i]))
  v18maxpos=max(where(vL2 le v2[i]))
  v18=vL2[where(specL2 eq max(specL2[v18minpos:v18maxpos]))]
  ;print,v13,specL[v13],specL[where(specL eq max(specL))]
  ;PRINT,where(specL eq max(specL))
  ;print,MEAN(abs(specL2))*5, specL2[where(specL2 eq max(specL2))]
  
  if specL2[where(specL2 eq max(specL2))] lt MEAN(abs(specL2))*3.5 then begin
    v_c=v13
    printf,log,' Use v(13CO) as v_c:'
    printf,log,v_c
  endif else begin
    v_c=v18
    printf,log,' Use v(C18O) as v_c:'
    printf,log,v_c
  endelse

  cgplot,vU,specU,psym=10,xrange=vrange,color='blue',yrange=[0-max(specU)*0.15,max(specU)*1.15],$ ;title=source,xtitle='LSR Velocity (km s!U-1!N)',
    ytitle='T!DMB!N (K)',position=p0,xtickinterval=10,xtitle='LSR Velocity (km s!U-1!N)';,xtickformat='(a1)'
  ;cgtext,0.05*(!x.CRange[1]-!x.Crange[0])+!x.Crange[0],0.9*(!y.CRange[1]-!y.Crange[0])+!y.Crange[0],source;,charsize=0.8
  cgtext,v[i],(!y.CRange[0]*1.25)/2.0,'v(NH!D3!N)='+strcompress(string(v[i],format='(f5.1)'),/remove_all),alignment=0.5
  AL_Legend, items, LineStyle=lines, Color=colors, box=0, position=[0.575*(!x.CRange[1]-!x.Crange[0])+!x.Crange[0],0.95*(!y.CRange[1]-!y.Crange[0])+!y.Crange[0]]
  cgPlot, !x.CRange,[0,0], LineStyle=1, Color='black', /overplot,position=p0
  cgPlot, [v[i],v[i]],[0.,!y.CRange[1]], LineStyle=1, Color='black', /overplot,position=p0
  cgPlot, [v_c,v_c],[0.,!y.CRange[1]], LineStyle=1, Color='red', /overplot,position=p0
  ;print,v_c
  ;cgtext,v[i],(!y.CRange[0]+0.1)/2.0,'v(NH!D3!N)='+strcompress(string(v[i],format='(f5.1)'),/remove_all),alignment=0.5
  
  ;sfitU = GaussFit(vU, specU, coeffU, NTERMS=6)
  ;cgplot,vU,sfitU,color='black',/overplot
  ;cgplot,vL,specL,psym=10,color='green',/overplot,position=p0
  ;sfitL = GaussFit(vL, specL, coeffL, NTERMS=6)
  ;cgplot,vL,sfitL*max(specU)/max(specL),color='green',/overplot
  ;cgplot,vL2,specL2,psym=10,color='red',/overplot,position=p0
 
 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
  
 if LINE[i] eq '12CO' then begin
  baseline=make_array(n_elements(specU),value=0.)
   
  if bv1[i] eq -999 || bv2[i] eq -999 then begin 
    print, 'No blue wings.' 
  endif else begin
    bluelobe=[bv1[i],bv2[i]]
    ;bl_rg=coord2pix(shdrU,bluelobe*1000.0,3)
    bpolyrg=where(vU ge bv1[i] and vU le bv2[i])
    bx1_arr=array_band(vU[bpolyrg],/band)
    bx2_arr=array_band(Reverse(vU[bpolyrg]),/band)
    by1_arr=array_band(specU[bpolyrg])
    by2_arr=array_band(baseline[bpolyrg])
    cgColorFill,[bx1_arr, bx2_arr, bx1_arr[0]],[by1_arr,by2_arr,by1_arr[0]],Color='sky blue'
    endelse
  if rv1[i] eq -999 || rv2[i] eq -999 then begin
    print, 'No red wings.' 
  endif else begin
    redlobe=[rv1[i],rv2[i]]
    rpolyrg=where(vU ge rv1[i] and vU le rv2[i])
    ;rl_rg=coord2pix(shdrU,redlobe*1000.0,3)
    rx1_arr=array_band(vU[rpolyrg],/band)
    rx2_arr=array_band(Reverse(vU[rpolyrg]),/band)
    ry1_arr=array_band(specU[rpolyrg])
    ry2_arr=array_band(baseline[rpolyrg])
    cgColorFill,[rx1_arr, rx2_arr, rx1_arr[0]],[ry1_arr,ry2_arr,ry1_arr[0]],Color='pink'
  ;  cgColorFill, [vU[rl_rg[0]:rl_rg[1]], Reverse(vU[rl_rg[0]:rl_rg[1]]), vU[rl_rg[0]]], $
  ;     [specU[rl_rg[0]:rl_rg[1]],baseline[rl_rg[0]:rl_rg[1]],specU[rl_rg[0]]], $
  ;     Color='pink'
    endelse
  endif
  
 if LINE[i] eq '13CO' then begin
  baseline=make_array(n_elements(specL),value=0.)
  if bv1[i] eq -999 || bv2[i] eq -999 then begin 
    print, 'No blue wings.' 
  endif else begin
    bluelobe=[bv1[i],bv2[i]]
    bpolyrg=where(vL ge bv1[i] and vL le bv2[i])
    bx1_arr=array_band(vL[bpolyrg],/band)
    bx2_arr=array_band(Reverse(vL[bpolyrg]),/band)
    by1_arr=array_band(specL[bpolyrg])
    by2_arr=array_band(baseline[bpolyrg])
    cgColorFill,[bx1_arr, bx2_arr, bx1_arr[0]],[by1_arr,by2_arr,by1_arr[0]],Color='sky blue'
    endelse
  if rv1[i] eq -999 || rv2[i] eq -999 then begin
    print, 'No red wings.' 
  endif else begin
    redlobe=[rv1[i],rv2[i]]
    rpolyrg=where(vL ge rv1[i] and vL le rv2[i])
    rx1_arr=array_band(vL[rpolyrg],/band)
    rx2_arr=array_band(Reverse(vL[rpolyrg]),/band)
    ry1_arr=array_band(specL[rpolyrg])
    ry2_arr=array_band(baseline[rpolyrg])
    cgColorFill,[rx1_arr, rx2_arr, rx1_arr[0]],[ry1_arr,ry2_arr,ry1_arr[0]],Color='pink'
    endelse
   
  ;if bv1[i] eq -999 || bv2[i] eq -999 then begin 
  ;  print, 'No blue wings.' 
  ;endif else begin
  ;  bluelobe=[bv1[i],bv2[i]]
  ;  bl_rg=coord2pix(shdrL,bluelobe*1000.0,3)
  ;  cgColorFill, [vL[bl_rg[0]:bl_rg[1]], Reverse(vL[bl_rg[0]:bl_rg[1]]), vL[bl_rg[0]]], $
  ;   [specL[bl_rg[0]:bl_rg[1]],baseline[bl_rg[0]:bl_rg[1]],specL[bl_rg[0]]], $
  ;   Color='sky blue'
  ;  endelse
  ;if rv1[i] eq -999 || rv2[i] eq -999 then begin
  ;  print, 'No red wings.' 
  ;endif else begin
  ;  redlobe=[rv1[i],rv2[i]]
  ;  rl_rg=coord2pix(shdrL,redlobe*1000.0,3)
  ;  cgColorFill, [vL[rl_rg[0]:rl_rg[1]], Reverse(vL[rl_rg[0]:rl_rg[1]]), vL[rl_rg[0]]], $
  ;     [specL[rl_rg[0]:rl_rg[1]],baseline[rl_rg[0]:rl_rg[1]],specL[rl_rg[0]]], $
  ;     Color='pink'
  ;  endelse
  endif
  
  cgplot,vU,specU,psym=lines[0],/overplot,color='blue'
  cgplot,vL,specL,psym=lines[1],color='GREEN',/overplot;,position=p0
  cgplot,vL2,specL2,psym=lines[2],color='red',/overplot;,position=p0
  ;cgplot,vL2s,specL2s,psym=10,color='black',/overplot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  
;;;;;mapping and contour;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;CO DATA CHECK;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  print, '  Integrated area will be set as '+string(area[0],format='(f5.2)')+'*'+string(area[1],format='(f5.2)')+' arcmin square.'
  print, '  Integrated velocity range will be set between '+string(v1[i],format='(f+6.1)')+' to '+string(v2[i],format='(f+6.1)')+' km/s.'

  printf,log,' Cropping configuration of CO image:'
  printf,log,[l[i]+area[0]/120d,l[i]-area[0]/120d],[b[i]-area[1]/120d,b[i]+area[1]/120d]
  printf,log,' Area value:'
  printf,log,area/60d

  cropfits,DLHUname,[l[i]-area[0]/120d,l[i]+area[0]/120d],[b[i]-area[1]/120d,b[i]+area[1]/120d],[v1[i],v2[i]],output='tempU'
  fits_read,'tempU_C.fits',datU,hdrU
  iimapU=total(datU[*,*,*],3)
  imageU=bytscl(iimapU,min=0.1,max=max(iimapU))
  ;imageU=iimapU
  ;fits_write,'testmax.fits',imageU
  
  cropfits,DLHLname,[l[i]-area[0]/120d,l[i]+area[0]/120d],[b[i]-area[1]/120d,b[i]+area[1]/120d],[v1[i],v2[i]],output='tempL'
  fits_read,'tempL_C.fits',datL,hdrL
  iimapL=total(datL[*,*,*],3)
  imageL=bytscl(iimapL,min=0.1,max=max(iimapL))

  ;cropfits,DLHL2name,[l[i]-area[0]/120d,l[i]+area[0]/120d],[b[i]-area[1]/120d,b[i]+area[1]/120d],[v1[i],v2[i]],output='tempL2'
  ;fits_read,'tempL2_C.fits',datL2,hdrL2
  ;iimapL2=total(datL2[*,*,*],3)
  ;imageL2=bytscl(iimapL2,min=0.1,max=max(iimapL2)/1.2)
  
  ;img_R=reform(imageL2,1,n_elements(imageL2[*,0]),n_elements(imageL2[0,*]))
  ;img_G=reform(imageL,1,n_elements(imageL[*,0]),n_elements(imageL[0,*]))
  ;img_B=reform(imageU,1,n_elements(imageU[*,0]),n_elements(imageU[0,*]))
  ;imcolor=[img_R,img_G,img_B]
  
;;;WISE DATA CHECK;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  wiseBname='wise_4.6_'+source+'.fits'
  wiseGname='wise_12_'+source+'.fits'
  wiseRname='wise_22_'+source+'.fits'
  wiseRGBfits=[wiseBname,wiseGname,wiseRname]
  wisesurvey=['WISE 4.6','WISE 12','WISE 22']
  foreach wisename, wiseRGBfits, index do begin
    if ~file_test('./wisedata/'+wisename) then begin
      print,' Downloading '+wisename+' from SURVEY '+wisesurvey[index]
      shellcmd='./skvbatch_wget file='+wisename+" position='"+string([l[i],b[i]],format='(f7.3,",",f6.3)')+"' Survey='"+wisesurvey[index]+$
        "' Coordinates='Galactic' Projection='Car' Pixels=600 Size="+string([szmax/60.0,szmax/60.0],format='(f7.5,",",f7.5)')
      ;print,shellcmd
      spawn,shellcmd
      fileinfo=file_info(wisename)
      while fileinfo.size lt 1000000 do begin
        spawn, 'rm ./'+wisename
        spawn,shellcmd
        fileinfo=file_info(wisename)
      endwhile
      spawn,'mv wise*.fits ./wisedata/'
    endif
  endforeach
  
  whdr=headfits('./wisedata/'+wiseBname)
  
  crp1=sxpar(hdrU,'CRPIX1')
  crv1=sxpar(hdrU,'CRVAL1')
  del1=sxpar(hdrU,'CDELT1')
  l_l=(0.5-crp1)*del1+crv1
  l_r=(sxpar(hdrU,'NAXIS1')+0.5-crp1)*del1+crv1
  cl_l=l_l+sxpar(whdr,'CDELT1')
  cl_r=l_r-sxpar(whdr,'CDELT1')
  
  crp2=sxpar(hdrU,'CRPIX2')
  crv2=sxpar(hdrU,'CRVAL2')
  del2=sxpar(hdrU,'CDELT2')
  b_d=(0.5-crp2)*del2+crv2
  b_u=(sxpar(hdrU,'NAXIS2')+0.5-crp2)*del2+crv2
  cb_d=b_d+sxpar(whdr,'CDELT2')
  cb_u=b_u-sxpar(whdr,'CDELT2')
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;mkpv and lobe feature spectra;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

;readcol,'outflowpara2',sign2,blue_l,blue_b,red_l,red_b,bszx,bszy,rszx,rszy,/silent,stringskip='#',$
;  format='(A,F,F,F,F,F,F,F,F)'
;PVSLICE
j=where(sign2 eq sign1[i])
if j ne -1 then begin
  
;;;;;;;;;draw pv;;;;;;;;;;;;;;;;
  if LINE[i] eq '12CO' then begin
    if bv1[i] ne -999 && rv2[i] ne -999 then begin
    pvvmin=bv1[i]
    pvvmax=rv2[i]
    endif else if bv1[i] ne -999 && rv2[i] eq -999 then begin 
    pvvmin=bv1[i]
    pvvmax=2*v[i]-bv1[i]
    endif else if bv1[i] eq -999 && rv2[i] ne -999 then begin
    pvvmin=2*v[i]-rv2[i]
    pvvmax=rv2[i]
    endif else begin
    pvvmin=v[i]-20
    pvvmax=v[i]+20
    endelse  
    cropfits,DLHUname,[pvvmin,pvvmax],dim='v',output='pvtempU'
    ohdr=headfits('pvtempU_C.fits')
    
    ocrp1=sxpar(ohdr,'CRPIX1')
    ocrv1=sxpar(ohdr,'CRVAL1')
    odel1=sxpar(ohdr,'CDELT1')
    ol_l=(1-ocrp1)*odel1+ocrv1
    ol_r=(sxpar(ohdr,'NAXIS1')-ocrp1)*odel1+ocrv1
    
    ocrp2=sxpar(ohdr,'CRPIX2')
    ocrv2=sxpar(ohdr,'CRVAL2')
    odel2=sxpar(ohdr,'CDELT2')
    ob_d=(1-ocrp2)*odel2+ocrv2
    ob_u=(sxpar(ohdr,'NAXIS2')-ocrp2)*odel2+ocrv2
    
    l_pos0=(blue_l[j]+red_l[j])/2.0
    b_pos0=(blue_b[j]+red_b[j])/2.0
    l_pos=[(3.5*blue_l[j]-2.5*l_pos0)<ol_l,(3.5*red_l[j]-2.5*l_pos0)>ol_r]
    b_pos=[(3.5*blue_b[j]-2.5*b_pos0)>ob_d,(3.5*red_b[j]-2.5*b_pos0)<ob_u]

    mkpvslice,'pvtempU_C.fits',/gal,l_pos,b_pos
    fits_read,'pvslice.fits',pvdata,pvhdr
    ;pblue=sqrt((l_pos[0]-l_pos0)^2+(b_pos[0]-b_pos)^2)*60
    ;pred=sqrt((l_pos[1]-l_pos0)^2+(b_pos[1]-b_pos)^2)*60
    pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
    pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
    ;printf,log,' PV image vrange:'
    ;printf,log,[min(pvv),max(pvv)]
    pvvrg=[min(pvv),max(pvv)]   
    pvprg=[min(pvp),max(pvp)]*60.0
    printf,log,' PV image prange:'
    printf,log,pvprg
    cgplot,[0],[0],/nodata,/noerase,position=p1,xrange=pvvrg,yrange=pvprg,$
    ytickformat='(a1)',xtickformat='(a1)',yminor=5
    cgloadct,49,clip=[31,240]
    pvlevels=max(pvdata)*(0.05+indgen(10)*0.1)
    pvsmooth=smooth(pvdata,4,/EDGE_MIRROR)
    cgcontour,pvsmooth,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='blue'
    cgplot,[0],[0],/nodata,/noerase,position=p1,xrange=pvvrg,yrange=pvprg,$
    ytickinterval=1.0,ytitle="Position (')",yminor=5,xtitle='LSR Velocity (km s!U-1!N)'
    cgPlot, [v_c,v_c],[!y.CRange[0],!y.CRange[1]], LineStyle=1, Color='red', /overplot,position=p1
    cgPlot, [v[i],v[i]],[!y.CRange[0],!y.CRange[1]], LineStyle=1, Color='black', /overplot,position=p1    
    cgtext,!x.CRange[0]-0.05*(!x.CRange[0]-!x.CRange[1]),!y.CRange[0]+0.1*(!y.CRange[1]-!y.CRange[0]),"offset at "
    cgtext,!x.CRange[0]-0.05*(!x.CRange[0]-!x.CRange[1]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
    strcompress(string(l_pos[0],format='(f7.3)'),/remove_all)+' '+strcompress(string(b_pos[0],format='(f6.3)'),/remove_all)
    cgtext,!x.CRange[0]-0.05*(!x.CRange[0]-!x.CRange[1]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),"!U12!NCO (1-0)"
  endif
  
  if LINE[i] eq '13CO' then begin
    if bv1[i] ne -999 && rv2[i] ne -999 then begin
    pvvmin=bv1[i]
    pvvmax=rv2[i]
    endif else if bv1[i] ne -999 && rv2[i] eq -999 then begin 
    pvvmin=bv1[i]
    pvvmax=2*v[i]-bv1[i]
    endif else if bv1[i] eq -999 && rv2[i] ne -999 then begin
    pvvmin=2*v[i]-rv2[i]
    pvvmax=rv2[i]
    endif else begin
    pvvmin=v[i]-20
    pvvmax=v[i]+20
    endelse
    cropfits,DLHLname,[pvvmin,pvvmax],dim='v',output='pvtempL'
    ohdr=headfits('pvtempL_C.fits')
    
    ocrp1=sxpar(ohdr,'CRPIX1')
    ocrv1=sxpar(ohdr,'CRVAL1')
    odel1=sxpar(ohdr,'CDELT1')
    ol_l=(1-ocrp1)*odel1+ocrv1
    ol_r=(sxpar(ohdr,'NAXIS1')-ocrp1)*odel1+ocrv1
    
    ocrp2=sxpar(ohdr,'CRPIX2')
    ocrv2=sxpar(ohdr,'CRVAL2')
    odel2=sxpar(ohdr,'CDELT2')
    ob_d=(1-ocrp2)*odel2+ocrv2
    ob_u=(sxpar(ohdr,'NAXIS2')-ocrp2)*odel2+ocrv2
    
    l_pos0=(blue_l[j]+red_l[j])/2.0
    b_pos0=(blue_b[j]+red_b[j])/2.0
    l_pos=[(3.5*blue_l[j]-2.5*l_pos0)<ol_l,(3.5*red_l[j]-2.5*l_pos0)>ol_r]
    b_pos=[(3.5*blue_b[j]-2.5*b_pos0)>ob_d,(3.5*red_b[j]-2.5*b_pos0)<ob_u]

    mkpvslice,'pvtempL_C.fits',/gal,l_pos,b_pos
    fits_read,'pvslice.fits',pvdata,pvhdr
    pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
    pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
    pvvrg=[min(pvv),max(pvv)]
    pvprg=[min(pvp),max(pvp)]*60.0
    printf,log,' PV image prange:'
    printf,log,pvprg
    cgplot,[0],[0],/nodata,/noerase,position=p1,xrange=pvvrg,yrange=pvprg,$
    ytickformat='(i)',xtickformat='(a1)',yminor=5
    cgloadct,53,clip=[64,240]
    pvlevels=max(pvdata)*(0.05+indgen(10)*0.1)
    cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='green'
    cgplot,[0],[0],/nodata,/noerase,position=p1,xrange=pvvrg,yrange=pvprg,$
    ytickinterval=1.0,ytitle="Position (')",yminor=5,xtitle='LSR Velocity (km s!U-1!N)'
    cgPlot, [v_c,v_c],[!y.CRange[0],!y.CRange[1]], LineStyle=1, Color='red', /overplot,position=p1
    cgPlot, [v[i],v[i]],[!y.CRange[0],!y.CRange[1]], LineStyle=1, Color='black', /overplot,position=p1
    cgtext,!x.CRange[0]-0.05*(!x.CRange[0]-!x.CRange[1]),!y.CRange[0]+0.1*(!y.CRange[1]-!y.CRange[0]),"offset at "
    cgtext,!x.CRange[0]-0.05*(!x.CRange[0]-!x.CRange[1]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
    strcompress(string(l_pos[0],format='(f7.3)'),/remove_all)+' '+strcompress(string(b_pos[0],format='(f6.3)'),/remove_all)
    cgtext,!x.CRange[0]-0.05*(!x.CRange[0]-!x.CRange[1]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),"!U13!NCO (1-0)"
    ;cgContour,riicon,label=0,/onimage,level=levels1,c_linestyle=2,Color='red'
  endif


;;;;;;;;;;;Extract lobe spectra;;;;;;;;;;;;;;;;
;
;readcol,'outflowpara2',sign2,blue_l,blue_b,red_l,red_b,bszx,bszy,rszx,rszy,/silent,stringskip='#',$
;  format='(A,F,F,F,F,F,F,F,F)'
  
  ssz=15/3600.0
  ;print,ssz
  
  if LINE[i] eq '12CO' then begin
    
    ;bluelobe
    if bv1[i] eq -999 || bv2[i] eq -999 then begin 
     print, 'No blue wings.'
     cropfits,DLHUname,[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],[v[i]-v_width[i],v[i]+v_width[i]],output='bstempU'
     fits_read,'bstempU_C.fits',bsdatU,bshdrU
     cropfits,DLHUname,[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],[v1[i],v2[i]],output='bctempU'
     fits_read,'bctempU_C.fits',bcdatU,bchdrU
     ;[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0]
     ;bc_rg=coord2pix(bchdrU,bluelobe*1000.0,3)
     ;bluecimage=bcdatU[*,*,bc_rg[0]:bc_rg[1]]
     bluecoord=[blue_l[j],blue_b[j]]
     printf,log,' Reference point of bluelobe:
     printf,log,bluecoord
     bluemax_x=round(coord2pix(bshdrU,bluecoord[0],1))
     bluemax_y=round(coord2pix(bshdrU,bluecoord[1],2))
     bspecU=reform(bsdatU[bluemax_x,bluemax_y,*])
     bspecU=boxsmooth(bspecU,v_bin[i])
     
     cgplot,vU,bspecU,psym=10,xrange=vrange,color='blue',yrange=[0-max(specU)*0.15,max(specU)*1.15],$
      ytitle='T!DMB!N (K)',position=p4,xtickinterval=10,/noerase,xtickformat='(a1)',yminor=5
     cgtext,!x.CRange[0]-0.05*(!x.CRange[0]-!x.CRange[1]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),"!U12!NCO (1-0)"
     AL_Legend, ['blue lobe',' red lobe'], LineStyle=[0,0], Color=['blue','red'], box=0,$
      position=[0.575*(!x.CRange[1]-!x.Crange[0])+!x.Crange[0],0.95*(!y.CRange[1]-!y.Crange[0])+!y.Crange[0]]
    endif else begin
     cropfits,DLHUname,[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],[v[i]-v_width[i],v[i]+v_width[i]],output='bstempU'
     fits_read,'bstempU_C.fits',bsdatU,bshdrU
    
     cropfits,DLHUname,[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],[v1[i],v2[i]],output='bctempU'
     fits_read,'bctempU_C.fits',bcdatU,bchdrU
     
     bc_rg=coord2pix(bchdrU,bluelobe*1000.0,3)
     bluecimage=bcdatU[*,*,bc_rg[0]:bc_rg[1]]
     bluecoord=maxcoord(bluecimage,bchdrU)
     printf,log,' Peak point of bluelobe:
     printf,log,bluecoord

     bluemax_x=round(coord2pix(bshdrU,bluecoord[0],1))
     bluemax_y=round(coord2pix(bshdrU,bluecoord[1],2))
     bspecU=reform(bsdatU[bluemax_x,bluemax_y,*])
     bspecU=boxsmooth(bspecU,v_bin[i])
     cgplot,vU,bspecU,psym=10,xrange=vrange,color='blue',yrange=[0-max(bspecU)*0.15,max(bspecU)*1.15],$
       ytitle='T!DMB!N (K)',position=p4,xtickinterval=10,/noerase,xtickformat='(a1)',yminor=5
     cgtext,!x.CRange[0]-0.05*(!x.CRange[0]-!x.CRange[1]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),"!U12!NCO (1-0)"
     AL_Legend, ['blue lobe',' red lobe'], LineStyle=[0,0], Color=['blue','red'], box=0,$
       position=[0.575*(!x.CRange[1]-!x.Crange[0])+!x.Crange[0],0.95*(!y.CRange[1]-!y.Crange[0])+!y.Crange[0]]
     endelse
    
    ;redlobe
  if rv1[i] eq -999 || rv2[i] eq -999 then begin
    print, 'No red wings.'
    cropfits,DLHUname,[red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0],[red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0],[v[i]-v_width[i],v[i]+v_width[i]],output='rstempU'
    fits_read,'rstempU_C.fits',rsdatU,rshdrU
    cropfits,DLHUname,[red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0],[red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0],[v1[i],v2[i]],output='rctempU'
    fits_read,'rctempU_C.fits',rcdatU,rchdrU
    redcoord=[red_l[j],red_b[j]]
    printf,log,' Reference point of redlobe:
    printf,log,redcoord
    redmax_x=round(coord2pix(rshdrU,redcoord[0],1))
    redmax_y=round(coord2pix(rshdrU,redcoord[1],2))
    rspecU=reform(rsdatU[redmax_x,redmax_y,*])
    rspecU=boxsmooth(rspecU,v_bin[i])
    cgplot,vU,rspecU,psym=10,xrange=vrange,color='red',position=p4,xtickinterval=10,/overplot,yminor=5
  endif else begin
    cropfits,DLHUname,[red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0],[red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0],[v[i]-v_width[i],v[i]+v_width[i]],output='rstempU'
    fits_read,'rstempU_C.fits',rsdatU,rshdrU    
    cropfits,DLHUname,[red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0],[red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0],[v1[i],v2[i]],output='rctempU'
    fits_read,'rctempU_C.fits',rcdatU,rchdrU   
    rc_rg=coord2pix(rchdrU,redlobe*1000.0,3)
    redcimage=rcdatU[*,*,rc_rg[0]:rc_rg[1]]
    redcoord=maxcoord(redcimage,rchdrU)
    printf,log,' Peak point of redlobe:
    printf,log,redcoord
    redmax_x=round(coord2pix(rshdrU,redcoord[0],1))
    redmax_y=round(coord2pix(rshdrU,redcoord[1],2))
    rspecU=reform(rsdatU[redmax_x,redmax_y,*])
    rspecU=boxsmooth(rspecU,v_bin[i])
    cgplot,vU,rspecU,psym=10,xrange=vrange,color='red',position=p4,xtickinterval=10,/overplot,yminor=5
    endelse

    cgPlot, !x.CRange,[0.1,0.1], LineStyle=1, Color='black', /overplot,position=p4
    cgPlot, [v_c,v_c],[!y.CRange[0],!y.CRange[1]], LineStyle=1, Color='red', /overplot,position=p4
    cgPlot, [v[i],v[i]],[!y.CRange[0],!y.CRange[1]], LineStyle=1, Color='black', /overplot,position=p4    
    ;cgtext,v[i],(!y.CRange[0]+0.1)/2.0,'v(NH!D3!N)='+strcompress(string(v[i],format='(f5.1)'),/remove_all),alignment=0.5
  endif
  
  if LINE[i] eq '13CO' then begin
    ;bluelobe
    
    if bv1[i] eq -999 || bv2[i] eq -999 then begin 
     print, 'No blue wings.'
    cropfits,DLHLname,[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],[v[i]-v_width[i],v[i]+v_width[i]],output='bstempL'
    fits_read,'bstempL_C.fits',bsdatL,bshdrL
    cropfits,DLHLname,[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],[v1[i],v2[i]],output='bctempL'
    fits_read,'bctempL_C.fits',bcdatL,bchdrL
    bluecoord=[blue_l[j],blue_b[j]]
    printf,log,' Reference point of bluelobe:
    printf,log,bluecoord
    bluemax_x=round(coord2pix(bshdrL,bluecoord[0],1))
    bluemax_y=round(coord2pix(bshdrL,bluecoord[1],2))
    bspecL=reform(bsdatL[bluemax_x,bluemax_y,*])
    bspecL=boxsmooth(bspecL,v_bin[i])
    cgplot,vL,bspecL,psym=10,xrange=vrange,color='blue',yrange=[0-max(specL)*0.15,max(specL)*1.15],$
      ytitle='T!DMB!N (K)',position=p4,xtickinterval=10,/noerase,xtickformat='(a1)',yminor=5
    cgtext,!x.CRange[0]-0.05*(!x.CRange[0]-!x.CRange[1]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),"!U13!NCO (1-0)"
    AL_Legend, ['blue lobe',' red lobe'], LineStyle=[0,0], Color=['blue','red'], box=0,$
      position=[0.575*(!x.CRange[1]-!x.Crange[0])+!x.Crange[0],0.95*(!y.CRange[1]-!y.Crange[0])+!y.Crange[0]]
    endif else begin    
    cropfits,DLHLname,[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],[v[i]-v_width[i],v[i]+v_width[i]],output='bstempL'
    fits_read,'bstempL_C.fits',bsdatL,bshdrL
    
    cropfits,DLHLname,[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],[v1[i],v2[i]],output='bctempL'
    fits_read,'bctempL_C.fits',bcdatL,bchdrL
    
    bc_rg=coord2pix(bchdrL,bluelobe*1000.0,3)
    bluecimage=bcdatL[*,*,bc_rg[0]:bc_rg[1]]
    bluecoord=maxcoord(bluecimage,bchdrL)
    printf,log,' Peak point of bluelobe:
    printf,log,bluecoord;,format='(f9.3)'
    ;print,[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],[v1[i],v2[i]]

    bluemax_x=round(coord2pix(bshdrL,bluecoord[0],1))
    bluemax_y=round(coord2pix(bshdrL,bluecoord[1],2))
    bspecL=reform(bsdatL[bluemax_x,bluemax_y,*])
    bspecL=boxsmooth(bspecL,v_bin[i])
    cgplot,vL,bspecL,psym=10,xrange=vrange,color='blue',yrange=[0-max(bspecL)*0.15,max(bspecL)*1.15],$
      ytitle='T!DMB!N (K)',position=p4,xtickinterval=10,/noerase,xtickformat='(a1)',yminor=5
    cgtext,!x.CRange[0]-0.05*(!x.CRange[0]-!x.CRange[1]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),"!U13!NCO (1-0)"
    AL_Legend, ['blue lobe',' red lobe'], LineStyle=[0,0], Color=['blue','red'], box=0,$
      position=[0.575*(!x.CRange[1]-!x.Crange[0])+!x.Crange[0],0.95*(!y.CRange[1]-!y.Crange[0])+!y.Crange[0]]
    
    endelse
    
    ;redlobe
    if rv1[i] eq -999 || rv2[i] eq -999 then begin
    print, 'No red wings.'
    cropfits,DLHLname,[red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0],[red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0],[v[i]-v_width[i],v[i]+v_width[i]],output='rstempL'
    fits_read,'rstempL_C.fits',rsdatL,rshdrL
    cropfits,DLHLname,[red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0],[red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0],[v1[i],v2[i]],output='rctempL'
    fits_read,'rctempL_C.fits',rcdatL,rchdrL
    redcoord=[red_l[j],red_b[j]]
    printf,log,' Reference point of redlobe:
    printf,log,redcoord
    redmax_x=round(coord2pix(rshdrL,redcoord[0],1))
    redmax_y=round(coord2pix(rshdrL,redcoord[1],2))
    rspecL=reform(rsdatL[redmax_x,redmax_y,*])
    rspecL=boxsmooth(rspecL,v_bin[i])
    cgplot,vL,rspecL,psym=10,xrange=vrange,color='red',position=p4,xtickinterval=10,/overplot,yminor=5
    endif else begin
    cropfits,DLHLname,[red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0],[red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0],[v[i]-v_width[i],v[i]+v_width[i]],output='rstempL'
    fits_read,'rstempL_C.fits',rsdatL,rshdrL    
    cropfits,DLHLname,[red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0],[red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0],[v1[i],v2[i]],output='rctempL'
    fits_read,'rctempL_C.fits',rcdatL,rchdrL    
    rc_rg=coord2pix(rchdrL,redlobe*1000.0,3)
    redcimage=rcdatL[*,*,rc_rg[0]:rc_rg[1]]
    redcoord=maxcoord(redcimage,rchdrL)

    printf,log,' Peak point of redlobe:
    printf,log,redcoord

    redmax_x=round(coord2pix(rshdrL,redcoord[0],1))
    redmax_y=round(coord2pix(rshdrL,redcoord[1],2))
    rspecL=reform(rsdatL[redmax_x,redmax_y,*])
    rspecL=boxsmooth(rspecL,v_bin[i])

    cgplot,vL,rspecL,psym=10,xrange=vrange,color='red',position=p4,xtickinterval=10,/overplot,yminor=5
    endelse
    
    cgPlot, !x.CRange,[0.1,0.1], LineStyle=1, Color='black', /overplot,position=p4
    cgPlot, [v_c,v_c],[!y.CRange[0],!y.CRange[1]], LineStyle=1, Color='red', /overplot,position=p4
    cgPlot, [v[i],v[i]],[!y.CRange[0],!y.CRange[1]], LineStyle=1, Color='black', /overplot,position=p4
  endif
    

endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;draw CO and WISE;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;lon = (sxpar(hdrU,'CRVAL1')+(dindgen(sxpar(hdrU,'NAXIS1'))-sxpar(hdrU,'CRPIX1')+1)*sxpar(hdrU,'CDELT1') + 360) mod 360
  ;lat =  sxpar(hdrU,'CRVAL2')+(dindgen(sxpar(hdrU,'NAXIS2'))-sxpar(hdrU,'CRPIX2')+1)*sxpar(hdrU,'CDELT2')
  ;l_range=[max(lon),min(lon)]
  l_range=[l_l,l_r]
  ;b_range=[min(lat),max(lat)]
  b_range=[b_d,b_u]
  printf,log,' Box axis range of CO image:'
  printf,log,l_range,b_range
  
  cgplot,[0],[0],/nodata,/noerase,position=p2,aspect=area[1]/area[0],xrange=l_range,yrange=b_range,$
    xtickinterval=round((max(l_range)-min(l_range))*10*area[1]/area[0])/50d,ytickinterval=round((max(b_range)-min(b_range))*10)/50d,$
    AxisColor='black',ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')
  
  if LINE[i] eq '12CO' then begin
  ;levels0=[max(imageU)*0.5,max(imageU)*0.6,max(imageU)*0.7,max(imageU)*0.8,max(imageU)*0.9]
  cgloadct,0,clip=[0,255]
  ;iimapU=congrid(iimapU,3*n_elements(iimapU[*,0]),3*n_elements(iimapU[0,*]),cubic=-0.5,/minus_one)
  cgcontour,iimapU,/onimage,nlevels=16,label=0,/fill;,thick=3,level=levels0;c_linestyle=0
  for n=0, n_elements(iimapU[*,0])-1 do begin
    xgrid=pix2coord(hdrU,make_array(n_elements(iimapU[*,0]),value=n,/int),1);+15./3600.
    ygrid=pix2coord(hdrU,make_array(n_elements(iimapU[0,*]),/int,/index),2);+15./3600.
    cgplot,xgrid,ygrid,/overplot,psym=1,color='brown',symsize=0.5
  endfor
  
  endif
  
  if LINE[i] eq '13CO' then begin
  ;levels0=[max(imageU)*0.5,max(imageU)*0.6,max(imageU)*0.7,max(imageU)*0.8,max(imageU)*0.9]
  cgloadct,0,clip=[0,255]
  cgcontour,iimapL,/onimage,nlevels=16,label=1,thick=3,/fill;,level=levels0;c_linestyle=0
  ;;;;;;;;;;
  for n=0, n_elements(iimapL[*,0])-1 do begin
    xgrid=pix2coord(hdrL,make_array(n_elements(iimapL[*,0]),value=n,/int),1);+15./3600.
    ygrid=pix2coord(hdrL,make_array(n_elements(iimapL[0,*]),/int,/index),2);+15./3600.
    cgplot,xgrid,ygrid,/overplot,psym=1,color='brown',symsize=0.5
  endfor
  ;;;;;;;;;;
  endif
  
  cgplot,[0],[0], xrange=l_range,yrange=b_range,xtickinterval=round((max(l_range)-min(l_range))*10*area[1]/area[0])/50d,ytickinterval=round((max(b_range)-min(b_range))*10)/50d,$;, xtitle='Galactic Longlitude (!Uo!N)',$
  /nodata,AxisColor=cgcolor('brown'),/noerase,position=p2,aspect=area[1]/area[0],xtickformat='(a1)',ytickformat='(a1)';$
    ;,XTicklen=1.0, YTicklen=1.0;,xgridstyle=1,ygridstyle=1
  ;cgplot,35.2,-0.74,/overplot,psym=1,color='red'
  ;ygrid=pix2coord(hdrL,make_array(n_elements(iimapL[0,*]),/int,/index),2)+15./3600.
  ;print,ygrid
  ;cgplot,35.2,ygrid,/overplot,psym=1,color='red'
  ;;;;;;;;;;;
  ;;;;;;;;;;;
  ;sz_scale=0.5*(180.0/!pi)/(1000.0*d[i])
  sz_scale=1/60.0
  x0_pos=!x.CRange[1]+0.05*(!x.CRange[0]-!x.CRange[1])
  x1_pos=!x.CRange[1]+0.05*(!x.CRange[0]-!x.CRange[1])+sz_scale
  y0_pos=!y.CRange[1]-0.05*(!y.CRange[1]-!y.CRange[0])
  y1_pos=y0_pos
  x_txt=(x0_pos+x1_pos)/2
  y_txt=!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0])
  ;y_txt0=!y.CRange[1]-0.05*(!y.CRange[1]-!y.CRange[0])
  cgplot,[x0_pos,x1_pos],[y0_pos,y1_pos],/overplot,LineStyle=0, Color='yellow',thick=3
  cgtext,x_txt,y_txt,string(sz_scale*(!pi/180)*1000.0*d[i],format='(f4.2)')+' pc',alignment=0.5,color='yellow'
  ;cgtext,x_txt,y_txt0," 1' ",alignment=0.5,color='yellow'
  x_txt1=!x.CRange[0]-0.05*(!x.CRange[0]-!x.CRange[1])
  y_txt1=!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0])
  x_txt2=x_txt1
  y_txt2=!y.CRange[0]+0.1*(!y.CRange[1]-!y.CRange[0])
  cgtext,x_txt1,y_txt1,'from '+strcompress(string(v1[i],format='(f6.1)'),/remove_all)+' to '+strcompress(string(v2[i],format='(f6.1)'),/remove_all)+' km s!u-1!n',$
    color='yellow'
  ;cgtext,x_txt1,y_txt1,strcompress(string(v1[i],format='(f+6.1)'),/remove_all)+' to '+strcompress(string(v2[i],format='(f+6.1)'),/remove_all)+' km s!u-1!n',$
  ;  color='yellow'
  if LINE[i] eq '12CO' then cgtext,x_txt2,y_txt2,'!U12!NCO(1-0) contours',color='yellow'
  if LINE[i] eq '13CO' then cgtext,x_txt2,y_txt2,'!U13!NCO(1-0) contours',color='yellow'
  x_txt3=x_txt1
  y_txt3=!y.CRange[1]-0.075*(!y.CRange[1]-!y.CRange[0])
  cgtext,x_txt3,y_txt3,'G'+source,color='yellow'
  
  ;cgplot,[l[i]-fspecxsz,l[i]-fspecxsz],[b[i]-fspecysz,b[i]+fspecysz],/overplot,LineStyle=1,Color='white'
  ;cgplot,[l[i]+fspecxsz,l[i]+fspecxsz],[b[i]-fspecysz,b[i]+fspecysz],/overplot,LineStyle=1,Color='white'
  ;cgplot,[l[i]-fspecxsz,l[i]+fspecxsz],[b[i]-fspecysz,b[i]-fspecysz],/overplot,LineStyle=1,Color='white'
  ;cgplot,[l[i]-fspecxsz,l[i]+fspecxsz],[b[i]+fspecysz,b[i]+fspecysz],/overplot,LineStyle=1,Color='white'
  ;normalpoint1=CONVERT_COORD(l[i]-fspecxsz, b[i]-fspecysz, /DATA, /TO_NORMAL)
  ;normalpoint2=CONVERT_COORD(l[i]-fspecxsz, b[i]+fspecysz, /DATA, /TO_NORMAL)
  
  if j ne -1 then cgarrow,l_pos[0],b_pos[0],l_pos[1],b_pos[1],LineStyle=0, Color='white',/data,hsize=!D.X_SIZE/128.0
  ;print,blue_l[j],blue_b[j],bszx[j]/2.0,bszy[j]/2.0
  if j ne -1 then begin
  if bv1[i] ne -999 && bv2[i] ne -999 then begin
    cgplot,/overplot,[blue_l[j]-bszx[j]/2.0,blue_l[j]+bszx[j]/2.0,blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],$
      [blue_b[j]-bszy[j]/2.0,blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0,blue_b[j]+bszy[j]/2.0,blue_b[j]-bszy[j]/2.0],$
      color='white',thick=3
    cgplot,/overplot,[blue_l[j]-bszx[j]/2.0,blue_l[j]+bszx[j]/2.0,blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],$
      [blue_b[j]-bszy[j]/2.0,blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0,blue_b[j]+bszy[j]/2.0,blue_b[j]-bszy[j]/2.0],$
      color='blue'
  endif
  if rv1[i] ne -999 && rv2[i] ne -999 then begin
    cgplot,/overplot,[red_l[j]-rszx[j]/2.0,red_l[j]+rszx[j]/2.0,red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0,red_l[j]-rszx[j]/2.0],$
      [red_b[j]-rszy[j]/2.0,red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0,red_b[j]+rszy[j]/2.0,red_b[j]-rszy[j]/2.0],$
      color='white',thick=3
    cgplot,/overplot,[red_l[j]-rszx[j]/2.0,red_l[j]+rszx[j]/2.0,red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0,red_l[j]-rszx[j]/2.0],$
      [red_b[j]-rszy[j]/2.0,red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0,red_b[j]+rszy[j]/2.0,red_b[j]-rszy[j]/2.0],$
      color='red'
  endif
  endif
  
  if LINE[i] eq '12CO' then begin
    
    ;bluelobe
    if bv1[i] eq -999 || bv2[i] eq -999 then begin 
    print, 'No blue CONTOURS.' 
    ;scale_b=sqrt((l_pos[0]-l_pos0)^2+(b_pos[0]-b_pos0)^2)*0.75*(!pi/180)*1000.0*d[i]
    endif else begin
    bl_rg=coord2pix(hdrU,bluelobe*1000.0,3)
    blueimage=datU[*,*,bl_rg[0]:bl_rg[1]]
    ;bluecoord=maxcoord(blueimage,hdrU)
    biicon0=total(blueimage[*,*,*],3)*abs(sxpar(hdrU,'CDELT3'))/1000.0
    biisz=size(biicon0)
    ;biicon=temporary(congrid(biicon0,5*biisz[1],5*biisz[2],cubic=-0.5,/MINUS_ONE))
    biicon=temporary(congrid(biicon0,4*biisz[1],4*biisz[2],cubic=-0.3,/minus_one))
    bmax=max(biicon0[fix(biisz[1]/3):fix(2*biisz[1]/3),fix(biisz[2]/3):fix(2*biisz[2]/3)])
    levelsign=where(sign3 eq sign1[i])
    levelsign=long(mean(levelsign))
    if levelsign ne -1 then begin
      dbcon=bmax*dblev[levelsign]/10.0
      bcon0=bmax*blev0[levelsign]
    endif else begin
      dbcon=bmax*1.0/10.0
      bcon0=bmax*0.4
    endelse
    ;readcol,'outflowpara2',sign2,blue_l,blue_b,red_l,red_b,bszx,bszy,rszx,rszy,/silent,stringskip='#',$
    ;  format='(A,F,F,F,F,F,F,F,F)'
    ;calculate I(CO) M(H2) p E
    ;print,[blue_l[j]-bszx[j]/2.0,blue_l[j]+bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],bluelobe,j,sign2[j]
    if j ne -1 then begin
    cropfits,'tempU_C.fits',[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],bluelobe,output='bluelobe'
    fits_read,'bluelobe_C.fits',blue_dat,blue_hdr
    ;print,bluelobe
    ;print,n_elements(total(blue_dat,3))
    bluemap=cubefilter(blue_dat,bcon0*1000.0/abs(sxpar(blue_hdr,'CDELT3')))
    ;help,bluemap
    db_area=((0.5/60.0)*(!pi/180)*1000.0*d[i])^2
    blue_area=bluemap.cov*db_area
    ;print,bluemap.cov,db_area
    T_ex=30
    
    bcv3=double(sxpar(blue_hdr,'CRVAL3'))
    bcp3=double(sxpar(blue_hdr,'CRPIX3'))
    bdl3=double(sxpar(blue_hdr,'CDELT3'))
    
    bI_12=total(bluemap.dat)*abs(bdl3)/(1000.0*bluemap.cov)

    bspec_mean=total(bluemap.dat,1)
    bspec_mean=total(bspec_mean,1)/bluemap.cov
    ;print,n_elements(bspec_mean)
    blue_v=((indgen(n_elements(bspec_mean))+1-bcp3)*bdl3+bcv3)/1000.0
    ;print,blue_v[0],blue_v[-1],v_c
    ;print,blue_v-max(v_c)
    ;bspecv=bspec_mean
    mean_bv=total(bspec_mean*(blue_v-max(v_c)))/total(bspec_mean)
    mean_bv2=total(bspec_mean*((blue_v-max(v_c))^2))/total(bspec_mean)
    printf,log,' Mean bluelobe v:',mean_bv;,total(bspec_mean*abs(bdl3)/1000.0)
    
    bNH2=4.2e13*T_ex/(exp(-5.5/T_ex))*1e4*bI_12
    bMH2=1.0/6.3e19*bNH2*blue_area
    bPH2=bMH2*abs(mean_bv)
    bEH2=0.5*bMH2*mean_bv2*1.9818E30*1.0E6 ;[J=kg*m^2*s^-2]
    scale_b=sqrt((l_pos[0]-l[i])^2+(b_pos[0]-b[i])^2)*0.5*(!pi/180)*1000.0*d[i]
    ;scale_r=sqrt((l_pos[1]-l_pos0)^2+(b_pos[1]-b_pos0)^2)*0.5*(!pi/180)*1000.0*d[i]
    printf,paraout,sign1[i],source+'_b',blue_area,bI_12,bNH2,bMH2,bPH2,bEH2,scale_b,mean_bv,format='(a9,a21,f8.2,f10.2,e10.2,f9.2,f13.2,e10.2,f8.2,f8.2)'
    endif
    
    levels1=indgen(10)*dbcon+bcon0
    cgloadct,33,clip=[0,60],ncolors=10,/reverse
    printf,log,' Blue contours from '+strcompress(string(bcon0,format='(f7.2)'),/remove_all)+' [K km/s] with a step of '+strcompress(string(dbcon,format='(f7.2)'),/remove_all)+' [K km/s]'
    cgContour,biicon,position=p2,label=0,/onimage,level=levels1,c_colors=indgen(10);,Color='blue'
    endelse
    
    ;redlobe
    if rv1[i] eq -999 || rv2[i] eq -999 then begin 
    print, 'No red CONTOURS.' 
    ;scale_r=sqrt((l_pos[1]-l_pos0)^2+(b_pos[1]-b_pos0)^2)*0.75*(!pi/180)*1000.0*d[i]
    endif else begin
    rl_rg=coord2pix(hdrU,redlobe*1000.0,3)
    riicon0=total(datU[*,*,rl_rg[0]:rl_rg[1]],3)*abs(sxpar(hdrU,'CDELT3'))/1000.0
    riisz=size(riicon0)
    riicon=temporary(congrid(riicon0,4*riisz[1],4*riisz[2],cubic=-0.3,/MINUS_ONE))
    ;riicon=riicon0
    rmax=max(riicon0[fix(riisz[1]/3):fix(2*riisz[1]/3),fix(riisz[2]/3):fix(2*riisz[2]/3)])
    ;levels2=[rmax*0.4,rmax*0.5,rmax*0.6,rmax*0.7,rmax*0.8,rmax*0.9,rmax*0.98,rmax*1.1,rmax*1.3,rmax*1.5,rmax*1.7,rmax*2]
    levelsign=where(sign3 eq sign1[i])
    levelsign=long(mean(levelsign))
    if levelsign ne -1 then begin
      drcon=rmax*drlev[levelsign]/10.0
      rcon0=rmax*rlev0[levelsign]
    endif else begin
      drcon=rmax*1.0/10.0
      rcon0=rmax*0.4
    endelse
    
    if j ne -1 then begin
    cropfits,'tempU_C.fits',[red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0],[red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0],redlobe,output='redlobe'
    fits_read,'redlobe_C.fits',red_dat,red_hdr
    redmap=cubefilter(red_dat,rcon0*1000.0/abs(sxpar(red_hdr,'CDELT3')))
    dr_area=((0.5/60.0)*(!pi/180)*1000.0*d[i])^2
    red_area=redmap.cov*dr_area
    T_ex=30
    rcv3=double(sxpar(red_hdr,'CRVAL3'))
    rcp3=double(sxpar(red_hdr,'CRPIX3'))
    rdl3=double(sxpar(red_hdr,'CDELT3'))
    
    rI_12=total(redmap.dat)*abs(rdl3)/(1000.0*redmap.cov)
  
    rspec_mean=total(redmap.dat,1)
    rspec_mean=total(rspec_mean,1)/redmap.cov
    red_v=((indgen(n_elements(rspec_mean))+1-rcp3)*rdl3+rcv3)/1000.0
    mean_rv=total(rspec_mean*(red_v-max(v_c)))/total(rspec_mean)
    mean_rv2=total(rspec_mean*((red_v-max(v_c))^2))/total(rspec_mean)
    
    rNH2=4.2e13*T_ex/(exp(-5.5/T_ex))*1e4*rI_12
    rMH2=1.0/6.3e19*rNH2*red_area
    rPH2=rMH2*abs(mean_rv)
    rEH2=0.5*rMH2*mean_rv2*1.9818E30*1.0E6 ;[J=kg*m^2*s^-2]
    scale_r=sqrt((l_pos[1]-l[i])^2+(b_pos[1]-b[i])^2)*0.5*(!pi/180)*1000.0*d[i]
    printf,paraout,sign1[i],source+'_r',red_area,rI_12,rNH2,rMH2,rPH2,rEH2,scale_r,mean_rv,format='(a9,a21,f8.2,f10.2,e10.2,f9.2,f13.2,e10.2,f8.2,f8.2)'
    
    if bv1[i] eq -999 || bv2[i] eq -999 then begin
       time_scale=scale_r/abs(mean_rv)*3.086e13/(365*24*3600.0) ;[year]
    endif else time_scale=2.0*max([scale_b,scale_r])/(abs(mean_bv)+abs(mean_rv))*3.086e13/(365*24*3600.0) ;[year]
    printf,log,' Time scale:',time_scale,format='(a,e10.2,a)',' years'
    printf,log,'            ',time_scale*365*24*3600,format='(a,e10.2,a)',' seconds'
    printf,log,' Mean redlobe v:',mean_rv;,total(rspec_mean*abs(rdl3)/1000.0)
    
    endif
    
    levels2=indgen(10)*drcon+rcon0
    cgloadct,33,clip=[195,255],ncolors=10
    printf,log,' Red contours from '+strcompress(string(rcon0,format='(f7.2)'),/remove_all)+' [K km/s] with a step of '+strcompress(string(drcon,format='(f7.2)'),/remove_all)+' [K km/s]'
    cgContour,riicon,position=p2,label=0,/onimage,level=levels2,c_linestyle=2,c_colors=indgen(10)
    endelse
    
  endif
  
  if LINE[i] eq '13CO' then begin
        ;bluelobe
    if bv1[i] eq -999 || bv2[i] eq -999 then begin 
    print, 'No blue CONTOURS.' 
    endif else begin
    bl_rg=coord2pix(hdrL,bluelobe*1000.0,3)
    blueimage=datL[*,*,bl_rg[0]:bl_rg[1]]
    ;bluecoord=maxcoord(blueimage,hdrL)
    biicon0=total(blueimage[*,*,*],3)*abs(sxpar(hdrL,'CDELT3'))/1000.0
    biisz=size(biicon0)
    ;biicon=temporary(congrid(biicon0,5*biisz[1],5*biisz[2],cubic=-0.5,/MINUS_ONE))
    biicon=temporary(congrid(biicon0,4*biisz[1],4*biisz[2],cubic=-0.3,/minus_one))
    bmax=max(biicon0[fix(biisz[1]/3):fix(2*biisz[1]/3),fix(biisz[2]/3):fix(2*biisz[2]/3)])
    levelsign=where(sign3 eq sign1[i])
    levelsign=long(mean(levelsign))
    if levelsign ne -1 then begin
      dbcon=bmax*dblev[levelsign]/10.0
      bcon0=bmax*blev0[levelsign]
    endif else begin
      dbcon=bmax*1.0/10.0
      bcon0=bmax*0.4
    endelse

    if j ne -1 then begin
    cropfits,'tempL_C.fits',[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],bluelobe,output='bluelobe'
    fits_read,'bluelobe_C.fits',blue_dat,blue_hdr
    ;print,n_elements(total(blue_dat,3))
    
    cropfits,'tempU_C.fits',[blue_l[j]+bszx[j]/2.0,blue_l[j]-bszx[j]/2.0],[blue_b[j]-bszy[j]/2.0,blue_b[j]+bszy[j]/2.0],bluelobe,output='bluelobeU'
    fits_read,'bluelobeU_C.fits',blue_datU,blue_hdrU
    
    tau13=-alog(1-total(blue_dat,3)*abs(sxpar(blue_hdr,'CDELT3'))/(total(blue_datU,3)*abs(sxpar(blue_hdrU,'CDELT3'))))
    ;help,tau13,tau13/(1-exp(-tau13))
    ;print,tau13/(1-exp(-tau13))
    
    ;print,bluelobe
    bluemap=cubefilter(blue_dat,bcon0*1000.0/abs(sxpar(blue_hdr,'CDELT3')))
    
    ;help,bluemap
    db_area=((0.5/60.0)*(!pi/180)*1000.0*d[i])^2
    blue_area=bluemap.cov*db_area
    ;print,bluemap.cov
    ;help,bluemap.dat
    T_ex=30
    
    bcv3=double(sxpar(blue_hdr,'CRVAL3'))
    bcp3=double(sxpar(blue_hdr,'CRPIX3'))
    bdl3=double(sxpar(blue_hdr,'CDELT3'))
    
    ;bI_13=total(bluemap.dat,3)*abs(bdl3)/1000.0
    bI_13=total(bluemap.dat)*abs(bdl3)/(1000.0*bluemap.cov)
    ;bi_test=total(bluemap.dat)*abs(bdl3)/1000.0
    ;help,bI_13
    ;print,bI_13[0],bi_13[-1]

    bspec_mean=total(bluemap.dat,1)
    bspec_mean=total(bspec_mean,1)/bluemap.cov
    ;print,n_elements(bspec_mean)
    blue_v=((indgen(n_elements(bspec_mean))+1-bcp3)*bdl3+bcv3)/1000.0
    ;print,blue_v[0],blue_v[-1],v_c
    ;print,blue_v-max(v_c)
    ;bspecv=bspec_mean
    mean_bv=total(bspec_mean*(blue_v-max(v_c)))/total(bspec_mean)
    mean_bv2=total(bspec_mean*((blue_v-max(v_c))^2))/total(bspec_mean)
    printf,log,' Mean bluelobe v:',mean_bv;,total(bspec_mean*abs(bdl3)/1000.0)
    
    bNH2=1.21e20*bI_13/(1-exp(-5.29/Tex12))
    bMH2=1.0/6.3e19*bNH2*blue_area
    bPH2=bMH2*abs(mean_bv)
    bEH2=0.5*bMH2*mean_bv2*1.9818E30*1.0E6 ;[J=kg*m^2*s^-2]
    scale_b=sqrt((l_pos[0]-l[i])^2+(b_pos[0]-b[i])^2)*0.5*(!pi/180)*1000.0*d[i]
    ;scale_r=sqrt((l_pos[1]-l_pos0)^2+(b_pos[1]-b_pos0)^2)*0.5*(!pi/180)*1000.0*d[i]
    printf,paraout,sign1[i],source+'_b13',blue_area,bI_13,bNH2,bMH2,bPH2,bEH2,scale_b,mean_bv,format='(a9,a21,f8.2,f10.2,e10.2,f9.2,f13.2,e10.2,f8.2,f8.2)'
    endif
    
    levels1=indgen(10)*dbcon+bcon0
    ;print,levels1[0],bcon0
    cgloadct,33,clip=[0,60],ncolors=10,/reverse
    printf,log,' Blue contours from '+strcompress(string(bcon0,format='(f7.2)'),/remove_all)+' [K km/s] with a step of '+strcompress(string(dbcon,format='(f7.2)'),/remove_all)+' [K km/s]'
    cgContour,biicon,position=p2,label=0,/onimage,level=levels1,c_colors=indgen(10);,Color='blue'
    endelse
    
    ;redlobe
    if rv1[i] eq -999 || rv2[i] eq -999 then begin 
    print, 'No red CONTOURS.' 
    endif else begin
    rl_rg=coord2pix(hdrL,redlobe*1000.0,3)
    riicon0=total(datL[*,*,rl_rg[0]:rl_rg[1]],3)*abs(sxpar(hdrL,'CDELT3'))/1000.0
    riisz=size(riicon0)
    riicon=temporary(congrid(riicon0,4*riisz[1],4*riisz[2],cubic=-0.3,/MINUS_ONE))
    ;riicon=riicon0
    rmax=max(riicon0[fix(riisz[1]/3):fix(2*riisz[1]/3),fix(riisz[2]/3):fix(2*riisz[2]/3)])
    ;levels2=[rmax*0.4,rmax*0.5,rmax*0.6,rmax*0.7,rmax*0.8,rmax*0.9,rmax*0.98,rmax*1.1,rmax*1.3,rmax*1.5,rmax*1.7,rmax*2]
    levelsign=where(sign3 eq sign1[i])
    levelsign=long(mean(levelsign))
    if levelsign ne -1 then begin
      drcon=rmax*drlev[levelsign]/10.0
      rcon0=rmax*rlev0[levelsign]
    endif else begin
      drcon=rmax*1.0/10.0
      rcon0=rmax*0.4
    endelse
    
    if j ne -1 then begin
    cropfits,'tempL_C.fits',[red_l[j]+rszx[j]/2.0,red_l[j]-rszx[j]/2.0],[red_b[j]-rszy[j]/2.0,red_b[j]+rszy[j]/2.0],redlobe,output='redlobe'
    fits_read,'redlobe_C.fits',red_dat,red_hdr
    redmap=cubefilter(red_dat,rcon0*1000.0/abs(sxpar(red_hdr,'CDELT3')))
    dr_area=((0.5/60.0)*(!pi/180)*1000.0*d[i])^2
    red_area=redmap.cov*dr_area
    T_ex=30
    rcv3=double(sxpar(red_hdr,'CRVAL3'))
    rcp3=double(sxpar(red_hdr,'CRPIX3'))
    rdl3=double(sxpar(red_hdr,'CDELT3'))
    
    rI_13=total(redmap.dat)*abs(rdl3)/(1000.0*redmap.cov)
  
    rspec_mean=total(redmap.dat,1)
    rspec_mean=total(rspec_mean,1)/redmap.cov
    red_v=((indgen(n_elements(rspec_mean))+1-rcp3)*rdl3+rcv3)/1000.0
    mean_rv=total(rspec_mean*(red_v-max(v_c)))/total(rspec_mean)
    mean_rv2=total(rspec_mean*((red_v-max(v_c))^2))/total(rspec_mean)
    
    ;rNH2=4.2e13*T_ex/(exp(-5.5/T_ex))*1e4*rI_12
    rNH2=1.21e20*rI_13/(1-exp(-5.29/Tex12))
    rMH2=1.0/6.3e19*rNH2*red_area
    rPH2=rMH2*abs(mean_rv)
    rEH2=0.5*rMH2*mean_rv2*1.9818E30*1.0E6 ;[J=kg*m^2*s^-2]
    scale_r=sqrt((l_pos[1]-l[i])^2+(b_pos[1]-b[i])^2)*0.5*(!pi/180)*1000.0*d[i]
    printf,paraout,sign1[i],source+'_r13',red_area,rI_13,rNH2,rMH2,rPH2,rEH2,scale_r,mean_rv,format='(a9,a21,f8.2,f10.2,e10.2,f9.2,f13.2,e10.2,f8.2,f8.2)'
    
    printf,log,' Mean redlobe v:',mean_rv;,total(rspec_mean*abs(rdl3)/1000.0)
    
    endif
    
    levels2=indgen(10)*drcon+rcon0
    cgloadct,33,clip=[195,255],ncolors=10
    printf,log,' Red contours from '+strcompress(string(rcon0,format='(f7.2)'),/remove_all)+' [K km/s] with a step of '+strcompress(string(drcon,format='(f7.2)'),/remove_all)+' [K km/s]'
    cgContour,riicon,position=p2,label=0,/onimage,level=levels2,c_linestyle=2,c_colors=indgen(10)
    endelse  


  endif
  
  if j ne -1 then begin
    if LINE[i] eq '12CO' then begin
     if (bv1[i] ne -999 && bv2[i] ne -999) && (rv1[i] ne -999 && rv2[i] ne -999) then begin
     ;printf,paracat,'sign','source','v_c','area_b','area_r','scale','delt_vb','delt_vr','I_b','I_r','N(H2)_b','N(H2)_r','M(H2)_b','M(H2)_r','M(H2)_out',$
     ;'p(H2)','E(H2)','time','M(H2)_out/t','Fm','Lm','M(H2)_core_18','Tex_12',$
     ;  format='(a12,a15,6a10,7a11,a13,a10,a8,a15,a20,a8,a15,a8)'
     scale_m=(scale_b+scale_r)*0.5
     scale_limit=(55/3600.0)*(!pi/180)*1000.0*d[i]
     time_scale=2.0*max([scale_b,scale_r])/(abs(mean_bv)+abs(mean_rv))*3.086e13/(365*24*3600.0) ;[year]
     ;printf,log,line[i]
     printf,log,' Time scale:',time_scale,format='(a,e10.2,a)',' years'
     printf,log,'            ',time_scale*365*24*3600,format='(a,e10.2,a)',' seconds'
     MH2_out=bMH2+rMH2
     PH2_out=bPH2+rPH2
     EH2_out=(bEH2+rEH2)*1e7
     M_out_rate=MH2_out/time_scale
     F_mech=PH2_out/time_scale
     L_mech=EH2_out/(time_scale*365*24*3600*3.828e33)
      printf,paracat,sign1[i],LINE[i],source,v_c,blue_area,red_area,scale_b,scale_r,mean_bv,mean_rv,bI_12,rI_12,bNH2,rNH2,bMH2,rMH2,MH2_out,$
        PH2_out,EH2_out,time_scale,M_out_rate,F_mech,L_mech,MH2_18,Tex12,scale_m,scale_limit,$
        format='(a10,a6,a15,7f10.2,2f11.2,2e11.2,3f11.2,f13.2,2e10.2,e15.2,e20.2,e10.2,f15.2,f10.2,2f10.2)'
     endif
     
     if ~(bv1[i] ne -999 && bv2[i] ne -999) && (rv1[i] ne -999 && rv2[i] ne -999) then begin
     scale_m=scale_r
     scale_limit=(55/3600.0)*(!pi/180)*1000.0*d[i]
     scale_b='\nodata'
     time_scale=scale_r/abs(mean_rv)*3.086e13/(365*24*3600.0) ;[year]
     ;printf,log,line[i]
     printf,log,' Time scale:',time_scale,format='(a,e10.2,a)',' years'
     printf,log,'            ',time_scale*365*24*3600,format='(a,e10.2,a)',' seconds'
     blue_area='\nodata'
     mean_bv='\nodata'
     bI_12='\nodata'
     bNH2='\nodata'
     bMH2='\nodata'
     MH2_out=rMH2
     PH2_out=rPH2
     EH2_out=rEH2*1e7
     M_out_rate=MH2_out/time_scale
     F_mech=PH2_out/time_scale
     L_mech=EH2_out/(time_scale*365*24*3600*3.828e33)
      printf,paracat,sign1[i],LINE[i],source,v_c,blue_area,red_area,scale_b,scale_r,mean_bv,mean_rv,bI_12,rI_12,bNH2,rNH2,bMH2,rMH2,MH2_out,$
        PH2_out,EH2_out,time_scale,M_out_rate,F_mech,L_mech,MH2_18,Tex12,scale_m,scale_limit,$
        format='(a10,a6,a15,f10.2,a10,f10.2,a10,f10.2,a10,f10.2,a11,f11.2,a11,e11.2,a11,2f11.2,f13.2,2e10.2,e15.2,e20.2,e10.2,f15.2,3f10.2)'
     endif
     
     if (bv1[i] ne -999 && bv2[i] ne -999) && ~(rv1[i] ne -999 && rv2[i] ne -999) then begin
     scale_m=scale_b
     scale_limit=(55/3600.0)*(!pi/180)*1000.0*d[i]
     scale_r='\nodata'
     time_scale=scale_b/abs(mean_bv)*3.086e13/(365*24*3600.0) ;[year]
     ;printf,log,line[i]
     printf,log,' Time scale:',time_scale,format='(a,e10.2,a)',' years'
     printf,log,'            ',time_scale*365*24*3600,format='(a,e10.2,a)',' seconds'
     red_area='\nodata'
     mean_rv='\nodata'
     rI_12='\nodata'
     rNH2='\nodata'
     rMH2='\nodata'
     MH2_out=bMH2
     PH2_out=bPH2
     EH2_out=bEH2*1e7
     M_out_rate=MH2_out/time_scale
     F_mech=PH2_out/time_scale
     L_mech=EH2_out/(time_scale*365*24*3600*3.828e33)
      printf,paracat,sign1[i],LINE[i],source,v_c,blue_area,red_area,scale_b,scale_r,mean_bv,mean_rv,bI_12,rI_12,bNH2,rNH2,bMH2,rMH2,MH2_out,$
        PH2_out,EH2_out,time_scale,M_out_rate,F_mech,L_mech,MH2_18,Tex12,scale_m,scale_limit,$
        format='(a10,a6,a15,2f10.2,a10,f10.2,a10,f10.2,a10,f11.2,a11,e11.2,a11,f11.2,a11,f11.2,f13.2,2e10.2,e15.2,e20.2,e10.2,f15.2,3f10.2)'
     endif
     
    endif
    
    if LINE[i] eq '13CO' then begin
     if (bv1[i] ne -999 && bv2[i] ne -999) && (rv1[i] ne -999 && rv2[i] ne -999) then begin
     ;printf,paracat,'sign','source','v_c','area_b','area_r','scale','delt_vb','delt_vr','I_b','I_r','N(H2)_b','N(H2)_r','M(H2)_b','M(H2)_r','M(H2)_out',$
     ;'p(H2)','E(H2)','time','M(H2)_out/t','Fm','Lm','M(H2)_core_18','Tex_12',$
     ;  format='(a12,a15,6a10,7a11,a13,a10,a8,a15,a20,a8,a15,a8)'
     scale_m=(scale_b+scale_r)*0.5
     scale_limit=(55/3600.0)*(!pi/180)*1000.0*d[i]
     time_scale=2.0*max([scale_b,scale_r])/(abs(mean_bv)+abs(mean_rv))*3.086e13/(365*24*3600.0) ;[year]
     ;printf,log,line[i]
     printf,log,' Time scale:',time_scale,format='(a,e10.2,a)',' years'
     printf,log,'            ',time_scale*365*24*3600,format='(a,e10.2,a)',' seconds'
     MH2_out=bMH2+rMH2
     PH2_out=bPH2+rPH2
     EH2_out=(bEH2+rEH2)*1e7
     M_out_rate=MH2_out/time_scale
     F_mech=PH2_out/time_scale
     L_mech=EH2_out/(time_scale*365*24*3600*3.828e33)
      printf,paracat,sign1[i],LINE[i],source,v_c,blue_area,red_area,scale_b,scale_r,mean_bv,mean_rv,bI_13,rI_13,bNH2,rNH2,bMH2,rMH2,MH2_out,$
        PH2_out,EH2_out,time_scale,M_out_rate,F_mech,L_mech,MH2_18,Tex12,scale_m,scale_limit,$
        format='(a10,a6,a15,7f10.2,2f11.2,2e11.2,3f11.2,f13.2,2e10.2,e15.2,e20.2,e10.2,f15.2,3f10.2)'
     endif
     
     if ~(bv1[i] ne -999 && bv2[i] ne -999) && (rv1[i] ne -999 && rv2[i] ne -999) then begin
     scale_m=scale_r
     scale_limit=(55/3600.0)*(!pi/180)*1000.0*d[i]
     scale_b='\nodata'
     time_scale=scale_r/abs(mean_rv)*3.086e13/(365*24*3600.0) ;[year]
     ;printf,log,line[i]
     printf,log,' Time scale:',time_scale,format='(a,e10.2,a)',' years'
     printf,log,'            ',time_scale*365*24*3600,format='(a,e10.2,a)',' seconds'
     blue_area='\nodata'
     mean_bv='\nodata'
     bI_13='\nodata'
     bNH2='\nodata'
     bMH2='\nodata'
     MH2_out=rMH2
     PH2_out=rPH2
     EH2_out=rEH2*1e7
     M_out_rate=MH2_out/time_scale
     F_mech=PH2_out/time_scale
     L_mech=EH2_out/(time_scale*365*24*3600*3.828e33)
      printf,paracat,sign1[i],LINE[i],source,v_c,blue_area,red_area,scale_b,scale_r,mean_bv,mean_rv,bI_13,rI_13,bNH2,rNH2,bMH2,rMH2,MH2_out,$
        PH2_out,EH2_out,time_scale,M_out_rate,F_mech,L_mech,MH2_18,Tex12,scale_m,scale_limit,$
        format='(a10,a6,a15,f10.2,a10,f10.2,a10,f10.2,a10,f10.2,a11,f11.2,a11,e11.2,a11,2f11.2,f13.2,2e10.2,e15.2,e20.2,e10.2,f15.2,3f10.2)'
     endif
     
     if (bv1[i] ne -999 && bv2[i] ne -999) && ~(rv1[i] ne -999 && rv2[i] ne -999) then begin
     scale_m=scale_b
     scale_limit=(55/3600.0)*(!pi/180)*1000.0*d[i]
     scale_r='\nodata'
     time_scale=scale_b/abs(mean_bv)*3.086e13/(365*24*3600.0) ;[year]
     ;printf,log,line[i]
     printf,log,' Time scale:',time_scale,format='(a,e10.2,a)',' years'
     printf,log,'            ',time_scale*365*24*3600,format='(a,e10.2,a)',' seconds'
     red_area='\nodata'
     mean_rv='\nodata'
     rI_13='\nodata'
     rNH2='\nodata'
     rMH2='\nodata'
     MH2_out=bMH2
     PH2_out=bPH2
     EH2_out=bEH2*1e7
     M_out_rate=MH2_out/time_scale
     F_mech=PH2_out/time_scale
     L_mech=EH2_out/(time_scale*365*24*3600*3.828e33)
      printf,paracat,sign1[i],LINE[i],source,v_c,blue_area,red_area,scale_b,scale_r,mean_bv,mean_rv,bI_13,rI_13,bNH2,rNH2,bMH2,rMH2,MH2_out,$
        PH2_out,EH2_out,time_scale,M_out_rate,F_mech,L_mech,MH2_18,Tex12,scale_m,scale_limit,$
        format='(a10,a6,a15,2f10.2,a10,f10.2,a10,f10.2,a10,f11.2,a11,e11.2,a11,f11.2,a11,f11.2,f13.2,2e10.2,e15.2,e20.2,e10.2,f15.2,3f10.2)'
     endif    
    endif
    
  endif
  
  cgplot,l[i],b[i],psym=cgSymCat(46),/overplot,symcolor=cgcolor('brown'),symsize=1.5,/data
  ;help,maxcoordsU
  if bv1[i] eq -999 || bv2[i] eq -999 then begin 
    print, 'No blue CONTOURS.' 
    endif else if j ne -1 then cgplot,bluecoord[0],bluecoord[1],psym=cgSymCat(16),/overplot,symcolor=cgcolor('sky blue'),symsize=1,/data
 
  if rv1[i] eq -999 || rv2[i] eq -999 then begin 
    print, 'No red CONTOURS.' 
    endif else if j ne -1 then cgplot,redcoord[0],redcoord[1],psym=cgSymCat(16),/overplot,symcolor=cgcolor('pink'),symsize=1,/data
  ;if line[i] eq '12CO' then cgplot,maxcoordsU[0],maxcoordsU[1],psym=cgSymCat(36,thick=2),/overplot,symcolor=cgcolor('BLU8'),symsize=0.8,/data
  ;if line[i] eq '13CO' then cgplot,maxcoordsL[0],maxcoordsL[1],psym=cgSymCat(36,thick=2),/overplot,symcolor=cgcolor('GRN8'),symsize=0.8,/data
  if line[i] eq '12CO' then cgplot,maxcoordsU[0],maxcoordsU[1],psym=cgSymCat(36,thick=2),/overplot,symcolor=cgcolor('BLU8'),symsize=0.8,/data
  cgplot,maxcoordsL[0],maxcoordsL[1],psym=cgSymCat(36,thick=2),/overplot,symcolor=cgcolor('GRN8'),symsize=0.8,/data
  cgplot,maxcoordsL2[0],maxcoordsL2[1],psym=cgSymCat(35,thick=2),/overplot,symcolor=cgcolor('RED8'),symsize=0.8,/data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;WISE;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if no wise than spitzer?;;;
  
  cropfits,'./wisedata/'+wiseBname,[cl_r,cl_l],[cb_u,cb_d],output='tempwiseB'
  cropfits,'./wisedata/'+wiseGname,[cl_r,cl_l],[cb_u,cb_d],output='tempwiseG'
  cropfits,'./wisedata/'+wiseRname,[cl_r,cl_l],[cb_u,cb_d],output='tempwiseR'

;  fits_read,'tempwiseB_C.fits',wiseB,whdrB
;  ;wiseimg_B=reform(bytscl(wiseB,min=max(wiseB)*0.1,max=max(wiseB)),1,n_elements(wiseB[*,0]),n_elements(wiseB[0,*]))
;  wiseimg_B=reform(logscl(wiseB,min=min(wiseB)*0.5,max=min(wiseB)+(max(wiseB)-min(wiseB))*0.8),1,n_elements(wiseB[*,0]),n_elements(wiseB[0,*]))
;  
;  fits_read,'tempwiseG_C.fits',wiseG
;  ;wiseimg_G=reform(bytscl(wiseG,min=max(wiseG)*0.1,max=max(wiseG)),1,n_elements(wiseG[*,0]),n_elements(wiseG[0,*]))
;  wiseimg_G=reform(logscl(wiseG,min=min(wiseG)*0.5,max=min(wiseG)+(max(wiseG)-min(wiseG))*0.8),1,n_elements(wiseG[*,0]),n_elements(wiseG[0,*]))
;  
;  fits_read,'tempwiseR_C.fits',wiseR
;  ;wiseimg_R=reform(bytscl(wiseR,min=max(wiseR)*0.1,max=max(wiseR)),1,n_elements(wiseR[*,0]),n_elements(wiseR[0,*]))
;  wiseimg_R=reform(logscl(wiseR,min=min(wiseR)*0.5,max=min(wiseR)+(max(wiseR)-min(wiseR))*0.8),1,n_elements(wiseR[*,0]),n_elements(wiseR[0,*]))
  
  fits_read,'tempwiseB_C.fits',wiseB,whdrB
  wiseimg_B=wisescl(wiseB)
  
  fits_read,'tempwiseG_C.fits',wiseG
  wiseimg_G=wisescl(wiseG)
  
  fits_read,'tempwiseR_C.fits',wiseR
  wiseimg_R=wisescl(wiseR)
  
  wiseimg_RGB=[wiseimg_R,wiseimg_G,wiseimg_B]
  
  wcrp1=sxpar(whdrB,'CRPIX1')
  wcrv1=sxpar(whdrB,'CRVAL1')
  wdel1=sxpar(whdrB,'CDELT1')
  wl_l=(0.5-wcrp1)*wdel1+wcrv1
  wl_r=(sxpar(whdrB,'NAXIS1')+0.5-wcrp1)*wdel1+wcrv1
  wcrp2=sxpar(whdrB,'CRPIX2')
  wcrv2=sxpar(whdrB,'CRVAL2')
  wdel2=sxpar(whdrB,'CDELT2')
  wb_d=(0.5-wcrp2)*wdel2+wcrv2
  wb_u=(sxpar(whdrB,'NAXIS2')+0.5-wcrp2)*wdel2+wcrv2
  printf,log,' Real range of WISE image:'
  printf,log,[wl_l,wl_r]
  printf,log,[wb_d,wb_u]


  cgplot,[0],[0],/nodata,/noerase,position=p3,aspect=area[1]/area[0],xrange=l_range,yrange=b_range,$
    xtickinterval=round((max(l_range)-min(l_range))*10*area[1]/area[0])/50d,ytickinterval=round((max(b_range)-min(b_range))*10)/50d,$
    AxisColor='black',ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')
  ;cgimage,imcolor,position = p2,/noerase
  cgimage,wiseimg_RGB,position = p3,/noerase

  cgplot,[0],[0], xrange=l_range,yrange=b_range,xtickinterval=round((max(l_range)-min(l_range))*10*area[1]/area[0])/50d,ytickinterval=round((max(b_range)-min(b_range))*10)/50d,$;, xtitle='Galactic Longitude (!Uo!N)',$
  /nodata,AxisColor=cgcolor('brown'),$;title='Integrated Intensity Map',
  /noerase,position=p3,aspect=area[1]/area[0],xtickformat='(a1)',ytickformat='(a1)'
  cgplot,[x0_pos,x1_pos],[y0_pos,y1_pos],/overplot,LineStyle=0, Color='yellow',thick=3
  cgtext,x_txt,y_txt,string(sz_scale*(!pi/180)*1000.0*d[i],format='(f4.2)')+' pc',alignment=0.5,color='yellow'
  cgtext,x_txt1,y_txt1,'WISE 4.6 12 22 !4l!Xm',color='yellow'
  cgtext,x_txt3,y_txt3,'G'+source,color='yellow'
  
  if j ne -1 then cgarrow,l_pos[0],b_pos[0],l_pos[1],b_pos[1],LineStyle=0, Color='white',/data,hsize=!D.X_SIZE/128.0
  
  
  if LINE[i] eq '12CO' then begin
    if bv1[i] eq -999 || bv2[i] eq -999 then begin 
      print, 'No blue CONTOURS.' 
    endif else begin
      cgloadct,33,clip=[0,60],ncolors=10,/reverse
      cgContour,biicon,position=p3,label=0,/onimage,level=levels1,c_colors=indgen(10)
    endelse
    if rv1[i] eq -999 || rv2[i] eq -999 then begin 
    print, 'No red CONTOURS.' 
    endif else begin
      cgloadct,33,clip=[195,255],ncolors=10
      cgContour,riicon,label=0,/onimage,level=levels2,c_linestyle=2,c_colors=indgen(10)
    endelse
  endif
  
  if LINE[i] eq '13CO' then begin
    if bv1[i] eq -999 || bv2[i] eq -999 then begin 
    print, 'No blue CONTOURS.' 
    endif else begin
      cgloadct,33,clip=[0,60],ncolors=10,/reverse
      cgContour,biicon,position=p3,c_colors=indgen(10),label=0,/onimage,level=levels1
    endelse
    if rv1[i] eq -999 || rv2[i] eq -999 then begin 
    print, 'No red CONTOURS.' 
    endif else begin
      cgloadct,33,clip=[195,255],ncolors=10
      cgContour,riicon,position=p3,c_colors=indgen(10),label=0,/onimage,c_linestyle=2,level=levels2
    endelse
  endif
  
  ;cgplot,l[i],b[i],psym=cgSymCat(45,thick=4.0),/overplot,symcolor=cgcolor('PUR8'),symsize=1.5,/data
  cgplot,l[i],b[i],psym=cgSymCat(46),/overplot,symcolor=cgcolor('brown'),symsize=1.5,/data
  
  ;cgtext,0.275,0.01,'(a)',/normal,alignment=0.5
  ;cgtext,0.75,0.01,'(b)',/normal,alignment=0.5
  ;cgarrow,normalpoint1[0],normalpoint1[1],p0[0],p0[1],/normal,hsize=0,thick=0.5
  ;cgarrow,normalpoint2[0],normalpoint2[1],p0[0],p0[3],/normal,hsize=0,thick=0.5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  
  cgps_close
  ;cgps2pdf,psname,/delete_ps
  spawn,'rm *temp*_C.fits'
  if file_test('*mosaic*.fits') then spawn,'rm *mosaic*.fits'
  if file_test('*_rms.fits') then spawn,'rm *_rms.fits'
  if file_test('*lobe_C.fits') then spawn,'rm *lobe_C.fits'
  if ~file_test('./ps_figure') then spawn,'mkdir ps_figure'
  if file_test(psname) then spawn,'mv '+psname+' ./ps_figure'
  print,'Done!'
endfor

free_lun,fitsunit,log,paraout,paracat


end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;MKPVSLICE;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mkpvslice, fitsfile, a, d, gal=gal, step=step
fits_read,fitsfile,dat,hdr
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
sxaddhist,'PV file: '+fitsfile,pvhdr
sxaddhist,'PV path:',pvhdr
for i=0,n_elements(a)-1 do sxaddhist,string(a[i])+' '+string(d[i]),pvhdr
sxaddhist,'Position in Degree',pvhdr
sxaddhist,'Velocity in km/s',pvhdr
fits_write,'pvslice.fits',slice,pvhdr

end