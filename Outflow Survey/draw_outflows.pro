function mkpvdat, dat, hdr, a, d, gal=gal, step=step, outdat=outdat
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
sxaddhist,'PV path:',pvhdr
for i=0,n_elements(a)-1 do sxaddhist,string(a[i])+' '+string(d[i]),pvhdr
sxaddhist,'Position in Degree',pvhdr
sxaddhist,'Velocity in km/s',pvhdr
outdat={dat:slice,hdr:pvhdr}
return, outdat
end

;pvbelt

pro draw_outflows

cd,'/home/lee/W3/'

;if ~keyword_set(region) then region='BFS52'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;distance=2
;if ~keyword_set(region) then region='GGMC1'
;if ~keyword_set(Lvrange) then Lvrange=[2,9]
;distance=2
;if ~keyword_set(region) then region='GGMC2'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;distance=2
;if ~keyword_set(region) then region='GGMC3'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;distance=2
if ~keyword_set(region) then region='region_C_III'
if ~keyword_set(Lvrange) then Lvrange=[-40,-20]
distance=2.0
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


cd,'/home/lee/W3/'+region+'/candidates'
L_lim=0.6
def_rmsU=0.25
def_rmsL=0.3
outsz=(8./distance)>3
outsz=outsz<5
fits_read,'../U_C.fits',datUa,hdrUa
fits_read,'../L_C.fits',datLa,hdrLa
;fits_read,region+'_L2a_mask.fits',datL2a,hdrL2a
fits_read,'../Lpeakv0.fits',peakvmap,vhdr
fits_read,'../Ubvmap.fits',Ubvmap,bvhdr
fits_read,'../Urvmap.fits',Urvmap,rvhdr
dathdr=list(datUa,hdrUa)
dathdrL=list(datLa,hdrLa)
;dathdrL2=list(datL2a,hdrL2a)

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

readcol,'blue_out.cat',bsn,bpeak1,bpeak2,bw0,bw1,format='I,F,F,F,F',stringskip='#'
readcol,'red_out.cat',rsn,rpeak1,rpeak2,rw0,rw1,format='I,F,F,F,F',stringskip='#' 

;readcol,'bv_range.cat',bsn1,vbv1,vbv2,format='I,F,F,F',stringskip='#'
;readcol,'rv_range.cat',bsn1,vrv1,vrv2,format='I,F,F,F',stringskip='#'

;bluelobe
for i=0, n_elements(bsn)-1 do begin
  specU0=reform(datUa[coord2pix(hdrUa,bpeak1[i],1),coord2pix(hdrUa,bpeak2[i],2),*])
  specU1=reform(datUa[coord2pix(hdrUa,bpeak1[i],1)-1,coord2pix(hdrUa,bpeak2[i],2),*])
  specU2=reform(datUa[coord2pix(hdrUa,bpeak1[i],1)+1,coord2pix(hdrUa,bpeak2[i],2),*])
  specU3=reform(datUa[coord2pix(hdrUa,bpeak1[i],1),coord2pix(hdrUa,bpeak2[i],2)-1,*])
  specU4=reform(datUa[coord2pix(hdrUa,bpeak1[i],1),coord2pix(hdrUa,bpeak2[i],2)+1,*])
  spec_U=(specU0+specU1+specU2+specU3+specU4)/5.
  specU=spec_U
  specL0=reform(datLa[coord2pix(hdrLa,bpeak1[i],1),coord2pix(hdrLa,bpeak2[i],2),*])
  specL1=reform(datLa[coord2pix(hdrLa,bpeak1[i],1)-1,coord2pix(hdrLa,bpeak2[i],2),*])
  specL2=reform(datLa[coord2pix(hdrLa,bpeak1[i],1)+1,coord2pix(hdrLa,bpeak2[i],2),*])
  specL3=reform(datLa[coord2pix(hdrLa,bpeak1[i],1),coord2pix(hdrLa,bpeak2[i],2)-1,*])
  specL4=reform(datLa[coord2pix(hdrLa,bpeak1[i],1),coord2pix(hdrLa,bpeak2[i],2)+1,*])
  spec_L=(specL0+specL1+specL2+specL3+specL4)/5.
  specL=spec_L
  ;vLpeak=peakvmap[coord2pix(hdrLa,bpeak1[i],1),coord2pix(hdrLa,bpeak2[i],2)]
  vLpeak=mean(peakvmap[(coord2pix(hdrLa,bpeak1[i],1)-1):(coord2pix(hdrLa,bpeak1[i],1)+1),$
      (coord2pix(hdrLa,bpeak2[i],2)-1):(coord2pix(hdrLa,bpeak2[i],2)+1)],/nan)
  ;print,vLpeak
  if finite(vLpeak,/nan) eq 1 then begin
    vLpeak=Lvrange[0]
    print,'ALEEEEEEEEEEEEEERT!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print,'b'
    print,i
  endif
;  if i eq 0 then begin
;    vrange=[-58.0,-22.0]
;  endif else begin
  vrange=[(Lvrange[1]-Ubvmap[coord2pix(bvhdr,bpeak1[i],1),coord2pix(bvhdr,bpeak2[i],2)]-3)<(vLpeak-5),$
    (Urvmap[coord2pix(rvhdr,bpeak1[i],1),coord2pix(rvhdr,bpeak2[i],2)]+Lvrange[0]+3)>(vLpeak+5)]
;  endelse  
;    
 ; vrange=[vbv1[i],vbv2[i]]
  ;vrange=[vLpeak-8,vLpeak+8]
  baseline=make_array(n_elements(specU),value=0)
  specpsname='blueout_'+num2str(bsn[i])+'.eps'
  cgps_open,specpsname,font=!p.font,/quiet,default_thickness=1.0,charsize=1.0,/portrait
  xsize1=1250. & ysize1=1000.
  cgDisplay, xsize=round(xsize1), ysize=round(ysize1)
  pos0=[100./xsize1,700./ysize1,600./xsize1,1-50./ysize1]
  pos1=[100./xsize1,100./ysize1,600./xsize1,600./ysize1]
  pos_r=(outsz-0.5)*sqrt(2)/outsz
  pos_d=850./(2*(pos_r+1))
  pos2=[700./xsize1,100./ysize1,1.-50./xsize1,(100.+pos_d)/ysize1]
  pos3=[700./xsize1,(100.+pos_d)/ysize1,1.-50./xsize1,(100.+pos_d*(1+pos_r))/ysize1]
  pos4=[700./xsize1,(100.+pos_d*(1+pos_r))/ysize1,1.-50./xsize1,(100.+pos_d*(2+pos_r))/ysize1]
  pos5=[700./xsize1,(100.+pos_d*(2+pos_r))/ysize1,1.-50./xsize1,950./ysize1]
  cgplot,vU,specU,psym=10,/nodata,xrange=vrange,color='blue',yrange=[0-max(spec_U)*0.15,max(spec_U)*1.15],$
    ytitle='T!DMB!N (K)',position=pos0,yminor=5,xtitle='LSR Velocity (km s!U-1!N)';,$
    ;title='OUTFLOW BLUELOBE CANDIDATE '+num2str(bsn[i])+': '+num2str(bpeak1[i],format='(f7.3)')+num2str(bpeak2[i],format='(f+6.3)')
  bpolyrg=where(vU ge bw0[i] and vU le bw1[i])
  bx1_arr=array_band(vU[bpolyrg],/band)
  bx2_arr=array_band(Reverse(vU[bpolyrg]),/band)
  by1_arr=array_band(spec_U[bpolyrg])
  by2_arr=array_band(baseline[bpolyrg])
  cgColorFill,[bx1_arr, bx2_arr, bx1_arr[0]],[by1_arr,by2_arr,by1_arr[0]],Color='Sky Blue'
  cgplot,vU,spec_U,psym=10,color='blue',/overplot
  cgplot,vL,spec_L,psym=10,color='green',/overplot
  cgPlot, !x.CRange,[0.,0.], LineStyle=0, Color='black',/overplot
  cgPlot,[vLpeak,vLpeak],[0,!y.crange[1]],linestyle=1,/over
  cropfits,/dataform,dathdr,[bpeak1[i]+outsz/60.,bpeak1[i]-outsz/60.],[bpeak2[i]-outsz/60.,bpeak2[i]+outsz/60.],[bw0[i],bw1[i]],output=crop
  bdata=smooth(crop.dat,[1,1,3],/edge_mirror)
  biimap=total(bdata,3)*abs(sxpar(crop.hdr,'CDELT3'))/1000.
  biimap=congrid(biimap,3*n_elements(biimap[*,0]),3*n_elements(biimap[0,*]),cubic=-0.5,/MINUS_ONE)
  biimap=smooth(biimap,[5,5],/edge_mirror)
  levels=max(biimap[round(n_elements(biimap[*,0])/2.-n_elements(biimap[*,0])/12.-1):round(n_elements(biimap[*,0])/2.+n_elements(biimap[*,0])/12.-1),$
    round(n_elements(biimap[0,*])/2.-n_elements(biimap[0,*])/12.-1):round(n_elements(biimap[0,*])/2.+n_elements(biimap[0,*])/12.-1)])$
    *(0.18+indgen(20)*0.05)
  crp1=sxpar(crop.hdr,'CRPIX1')
  crv1=sxpar(crop.hdr,'CRVAL1')
  del1=sxpar(crop.hdr,'CDELT1')
  axs1=sxpar(crop.hdr,'NAXIS1')
  l_l=(360.+(0.5-crp1)*del1+crv1) mod 360
  l_r=(360.+(axs1+0.5-crp1)*del1+crv1) mod 360
  crp2=sxpar(crop.hdr,'CRPIX2')
  crv2=sxpar(crop.hdr,'CRVAL2')
  del2=sxpar(crop.hdr,'CDELT2')
  axs2=sxpar(crop.hdr,'NAXIS2')
  b_d=(0.5-crp2)*del2+crv2
  b_u=(axs2+0.5-crp2)*del2+crv2
  x_range=[l_l,l_r]
  y_range=[b_d,b_u]
  ratio=abs((b_u-b_d)/(l_l-l_r))
  cgplot,[0],[0],xrange=x_range,yrange=y_range,$
    xtickinterval=0.05,ytickinterval=0.05,aspect=ratio,AxisColor='black',$
    ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})'),position=pos1,/noerase
  xpos=!x.crange[0]+(!x.crange[1]-!x.crange[0])*0.10
  ypos=!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.10
  cgloadct,49,ncolors=20
  cgcontour,biimap,/onimage,levels=levels,label=1,c_colors=indgen(20)
  x1_0=bpeak1[i]
  y1_0=bpeak2[i]-(outsz)/60
  x1_1=bpeak1[i]
  y1_1=bpeak2[i]+(outsz)/60
  cgarrow,x1_0,y1_0,x1_1,y1_1,/data,color='dark green',/clip
  ;cgplot,[x1_0-1/90.,x1_1-1/90.,x1_1+1/90.,x1_0+1/90.,x1_0-1/90.],$
  ;  [y1_0,y1_1,y1_1,y1_0,y1_0],color='dark green',/over
  x2_0=bpeak1[i]+(outsz-0.5)/60.
  y2_0=bpeak2[i]-(outsz-0.5)/60.
  x2_1=bpeak1[i]-(outsz-0.5)/60.
  y2_1=bpeak2[i]+(outsz-0.5)/60.
  cgarrow,x2_0,y2_0,x2_1,y2_1,/data,color='purple',/clip
  ;cgplot,[x2_0-1/90.,x2_1-1/90.,x2_1+1/90.,x2_0+1/90.,x2_0-1/90.],$
  ;  [y2_0-1/90.,y2_1-1/90.,y2_1+1/90.,y2_0+1/90.,y2_0-1/90.],color='purple',/over
  x3_0=bpeak1[i]+(outsz)/60.;x0+2*(x0-x1)
  y3_0=bpeak2[i]
  x3_1=bpeak1[i]-(outsz)/60.;x1+2*(x1-x0)
  y3_1=bpeak2[i]
  cgarrow,x3_0,y3_0,x3_1,y3_1,/data,color='brown',/clip
  ;cgplot,[x3_0,x3_1,x3_1,x3_0,x3_0],$
  ;  [y3_0-1/90.,y3_1-1/90.,y3_1+1/90.,y3_0+1/90.,y3_0-1/90.],color='brown',/over
  x4_0=bpeak1[i]+(outsz-0.5)/60.
  y4_0=bpeak2[i]+(outsz-0.5)/60.
  x4_1=bpeak1[i]-(outsz-0.5)/60.
  y4_1=bpeak2[i]-(outsz-0.5)/60.
  cgarrow,x4_0,y4_0,x4_1,y4_1,/data,color='blue',/clip
  ;cgplot,[x4_0+1/90.,x4_1+1/90.,x4_1-1/90.,x4_0-1/90.,x4_0+1/90.],$
  ;  [y4_0-1/90.,y4_1-1/90.,y4_1+1/90.,y4_0+1/90.,y4_0-1/90.],color='blue',/over
  cropfits,dathdr,vrange,dim='v',/dataform,output=pvout
  pv=mkpvbelt(pvout.dat,pvout.hdr,[x1_0,x1_1],[y1_0,y1_1],2,/gal,outdat=outdat)
  pvdata=outdat.dat
  pvhdr=outdat.hdr
  pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
  pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
  pvvrg=[min(pvv),max(pvv)]
  pvprg=[min(pvp),max(pvp)]*60.0-(outsz-0.5);*sqrt(2)   
  cgplot,[0],[0],/nodata,/noerase,position=pos2,xrange=pvvrg,yrange=pvprg,$
    ytickformat='(i)',xtickformat='(a1)',yminor=5
  cgloadct,50,clip=[64,240]
  pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
    round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
    *(0.05+indgen(10)*0.1)
  cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='dark green'
  cgplot,[0],[0],/nodata,/noerase,position=pos2,xrange=pvvrg,yrange=pvprg,$
    ytickinterval=2.0,ytitle="Position (')",yminor=5,xtitle='LSR Velocity (km s!U-1!N)'
  cgPlot, !x.crange,[0,0], LineStyle=1, /overplot
  cgPlot, [bw0[i],bw0[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot, [bw1[i],bw1[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over     
  pv=mkpvbelt(pvout.dat,pvout.hdr,[x2_0,x2_1],[y2_0,y2_1],2,/gal,outdat=outdat)
  pvdata=outdat.dat
  pvhdr=outdat.hdr
  pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
  pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
  pvvrg=[min(pvv),max(pvv)]
  pvprg=[min(pvp),max(pvp)]*60.0-(outsz-0.5)*sqrt(2)
  cgplot,[0],[0],/nodata,/noerase,position=pos3,xrange=pvvrg,yrange=pvprg,$
    ytickformat='(i)',xtickformat='(a1)',yminor=5,ytickinterval=2.0
  cgloadct,61,clip=[64,240]
  pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
    round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
    *(0.05+indgen(10)*0.1)
  cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='purple'
  cgplot,[0],[0],/nodata,/noerase,position=pos3,xrange=pvvrg,yrange=pvprg,$
    ytickinterval=2.0,ytitle="Position (')",yminor=5,xtickformat='(a1)'
  cgPlot, !x.crange,[0,0], LineStyle=1, /overplot
  cgPlot, [bw0[i],bw0[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot, [bw1[i],bw1[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
  pv=mkpvbelt(pvout.dat,pvout.hdr,[x3_0,x3_1],[y3_0,y3_1],2,/gal,outdat=outdat)
  pvdata=outdat.dat
  pvhdr=outdat.hdr
  pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
  pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
  pvvrg=[min(pvv),max(pvv)]
  pvprg=[min(pvp),max(pvp)]*60.0-(outsz-0.5);*sqrt(2)
  cgplot,[0],[0],/nodata,/noerase,position=pos4,xrange=pvvrg,yrange=pvprg,$
    ytickformat='(i)',xtickformat='(a1)',yminor=5,ytickinterval=2.0
  cgloadct,65,clip=[64,240]
  pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
    round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
    *(0.05+indgen(10)*0.1)
  cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='brown'
  cgplot,[0],[0],/nodata,/noerase,position=pos4,xrange=pvvrg,yrange=pvprg,$
    ytickinterval=2.0,ytitle="Position (')",yminor=5,xtickformat='(a1)'
  cgPlot, !x.crange,[0,0], LineStyle=1, /overplot
  cgPlot, [bw0[i],bw0[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot, [bw1[i],bw1[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
  pv=mkpvbelt(pvout.dat,pvout.hdr,[x4_0,x4_1],[y4_0,y4_1],2,/gal,outdat=outdat)
  pvdata=outdat.dat
  pvhdr=outdat.hdr
  pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
  pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
  pvvrg=[min(pvv),max(pvv)]
  pvprg=[min(pvp),max(pvp)]*60.0-(outsz-0.5)*sqrt(2)
  cgplot,[0],[0],/nodata,/noerase,position=pos5,xrange=pvvrg,yrange=pvprg,$
    ytickformat='(i)',xtickformat='(a1)',yminor=5,ytickinterval=2.0
  cgloadct,49,clip=[64,240]
  pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
    round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
    *(0.05+indgen(10)*0.1)
  cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='blue'
  cgplot,[0],[0],/nodata,/noerase,position=pos5,xrange=pvvrg,yrange=pvprg,$
    ytickinterval=2.0,ytitle="Position (')",yminor=5,xtickformat='(a1)'
  cgPlot, !x.crange,[0,0], LineStyle=1, /overplot
  cgPlot, [bw0[i],bw0[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot, [bw1[i],bw1[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
  cgps_close
endfor

;redlobe
for i=0, n_elements(rsn)-1 do begin
  ;if i ne 3 then continue
  specU0=reform(datUa[coord2pix(hdrUa,rpeak1[i],1),coord2pix(hdrUa,rpeak2[i],2),*])
  specU1=reform(datUa[coord2pix(hdrUa,rpeak1[i],1)-1,coord2pix(hdrUa,rpeak2[i],2),*])
  specU2=reform(datUa[coord2pix(hdrUa,rpeak1[i],1)+1,coord2pix(hdrUa,rpeak2[i],2),*])
  specU3=reform(datUa[coord2pix(hdrUa,rpeak1[i],1),coord2pix(hdrUa,rpeak2[i],2)-1,*])
  specU4=reform(datUa[coord2pix(hdrUa,rpeak1[i],1),coord2pix(hdrUa,rpeak2[i],2)+1,*])
  spec_U=(specU0+specU1+specU2+specU3+specU4)/5.
  specU=spec_U
  specL0=reform(datLa[coord2pix(hdrLa,rpeak1[i],1),coord2pix(hdrLa,rpeak2[i],2),*])
  specL1=reform(datLa[coord2pix(hdrLa,rpeak1[i],1)-1,coord2pix(hdrLa,rpeak2[i],2),*])
  specL2=reform(datLa[coord2pix(hdrLa,rpeak1[i],1)+1,coord2pix(hdrLa,rpeak2[i],2),*])
  specL3=reform(datLa[coord2pix(hdrLa,rpeak1[i],1),coord2pix(hdrLa,rpeak2[i],2)-1,*])
  specL4=reform(datLa[coord2pix(hdrLa,rpeak1[i],1),coord2pix(hdrLa,rpeak2[i],2)+1,*])
  spec_L=(specL0+specL1+specL2+specL3+specL4)/5.
  specL=spec_L
  ;vLpeak=peakvmap[coord2pix(hdrLa,rpeak1[i],1),coord2pix(hdrLa,rpeak2[i],2)]
  vLpeak=mean(peakvmap[(coord2pix(hdrLa,rpeak1[i],1)-1):(coord2pix(hdrLa,rpeak1[i],1)+1),$
      (coord2pix(hdrLa,rpeak2[i],2)-1):(coord2pix(hdrLa,rpeak2[i],2)+1)],/nan)
  ;print,vLpeak
  if finite(vLpeak,/nan) eq 1 then begin
    vLpeak=Lvrange[1]
    print,'ALEEEEEEEEEEEEEERT!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print,'r'
    print,i
  endif
;  if i eq 0 then begin
;    vrange=[-58.0,-22.0]
;  endif else begin
  vrange=[(Lvrange[1]-Ubvmap[coord2pix(bvhdr,rpeak1[i],1),coord2pix(bvhdr,rpeak2[i],2)]-3)<(vLpeak-5),$
    (Urvmap[coord2pix(rvhdr,rpeak1[i],1),coord2pix(rvhdr,rpeak2[i],2)]+Lvrange[0]+3)>(vLpeak+5)]
;  endelse  
  ;vrange=[vLpeak-8,vLpeak+8]
;  vrange=[vrv1[i],vrv2[i]]
  baseline=make_array(n_elements(specU),value=0)
  specpsname='redout_'+num2str(rsn[i])+'.eps'
  cgps_open,specpsname,font=!p.font,/quiet,default_thickness=1.0,charsize=1,/portrait
  xsize1=1250. & ysize1=1000.
  cgDisplay, xsize=round(xsize1), ysize=round(ysize1)
  pos0=[100./xsize1,700./ysize1,600./xsize1,1-50./ysize1]
  pos1=[100./xsize1,100./ysize1,600./xsize1,600./ysize1]
  pos_r=(outsz-0.5)*sqrt(2)/outsz
  pos_d=850./(2*(pos_r+1))
  pos2=[700./xsize1,100./ysize1,1.-50./xsize1,(100.+pos_d)/ysize1]
  pos3=[700./xsize1,(100.+pos_d)/ysize1,1.-50./xsize1,(100.+pos_d*(1+pos_r))/ysize1]
  pos4=[700./xsize1,(100.+pos_d*(1+pos_r))/ysize1,1.-50./xsize1,(100.+pos_d*(2+pos_r))/ysize1]
  pos5=[700./xsize1,(100.+pos_d*(2+pos_r))/ysize1,1.-50./xsize1,950./ysize1]
  cgplot,vU,specU,psym=10,/nodata,xrange=vrange,color='blue',yrange=[0-max(spec_U)*0.15,max(spec_U)*1.15],$
    ytitle='T!DMB!N (K)',position=pos0,yminor=5,xtitle='LSR Velocity (km s!U-1!N)';,$
    ;title='OUTFLOW REDLOBE CANDIDATE '+num2str(rsn[i])+': '+num2str(rpeak1[i],format='(f7.3)')+num2str(rpeak2[i],format='(f+6.3)')
  rpolyrg=where(vU ge rw0[i] and vU le rw1[i])
  rx1_arr=array_band(vU[rpolyrg],/band)
  rx2_arr=array_band(Reverse(vU[rpolyrg]),/band)
  ry1_arr=array_band(spec_U[rpolyrg])
  ry2_arr=array_band(baseline[rpolyrg])
  cgColorFill,[rx1_arr, rx2_arr, rx1_arr[0]],[ry1_arr,ry2_arr,ry1_arr[0]],Color='pink'
  cgplot,vU,spec_U,psym=10,color='blue',/overplot
  cgplot,vL,spec_L,psym=10,color='green',/overplot
  cgPlot, !x.CRange,[0.,0.], LineStyle=0, /overplot
  cgPlot,[vLpeak,vLpeak],[0,!y.crange[1]],linestyle=1,/over
  cropfits,/dataform,dathdr,[rpeak1[i]+outsz/60.,rpeak1[i]-outsz/60.],[rpeak2[i]-outsz/60.,rpeak2[i]+outsz/60.],[rw0[i],rw1[i]],output=crop
  rdata=smooth(crop.dat,[1,1,3],/edge_mirror)
  riimap=total(rdata,3)*abs(sxpar(crop.hdr,'CDELT3'))/1000.
  riimap=congrid(riimap,3*n_elements(riimap[*,0]),3*n_elements(riimap[0,*]),cubic=-0.5,/MINUS_ONE)
  riimap=smooth(riimap,[5,5],/edge_mirror)
  levels=max(riimap[round(n_elements(riimap[*,0])/2.-n_elements(riimap[*,0])/12.-1):round(n_elements(riimap[*,0])/2.+n_elements(riimap[*,0])/12.-1),$
    round(n_elements(riimap[0,*])/2.-n_elements(riimap[0,*])/12.-1):round(n_elements(riimap[0,*])/2.+n_elements(riimap[0,*])/12.-1)])$
    *(0.18+indgen(20)*0.05)
  crp1=sxpar(crop.hdr,'CRPIX1')
  crv1=sxpar(crop.hdr,'CRVAL1')
  del1=sxpar(crop.hdr,'CDELT1')
  axs1=sxpar(crop.hdr,'NAXIS1')
  l_l=(360.+(0.5-crp1)*del1+crv1) mod 360
  l_r=(360.+(axs1+0.5-crp1)*del1+crv1) mod 360
  crp2=sxpar(crop.hdr,'CRPIX2')
  crv2=sxpar(crop.hdr,'CRVAL2')
  del2=sxpar(crop.hdr,'CDELT2')
  axs2=sxpar(crop.hdr,'NAXIS2')
  b_d=(0.5-crp2)*del2+crv2
  b_u=(axs2+0.5-crp2)*del2+crv2
  x_range=[l_l,l_r]
  y_range=[b_d,b_u]
  ratio=abs((b_u-b_d)/(l_l-l_r))
  cgplot,[0],[0],xrange=x_range,yrange=y_range,$
    xtickinterval=0.05,ytickinterval=0.05,aspect=ratio,AxisColor='black',$
    ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})'),position=pos1,/noerase
  xpos=!x.crange[0]+(!x.crange[1]-!x.crange[0])*0.1
  ypos=!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.1
  cgloadct,62,ncolors=20
  cgcontour,riimap,/onimage,levels=levels,label=1,c_colors=indgen(20)
  x1_0=rpeak1[i]
  y1_0=rpeak2[i]-(outsz)/60
  x1_1=rpeak1[i]
  y1_1=rpeak2[i]+(outsz)/60
  cgarrow,x1_0,y1_0,x1_1,y1_1,/data,color='dark green',/clip
  ;cgplot,[x1_0-1/90.,x1_1-1/90.,x1_1+1/90.,x1_0+1/90.,x1_0-1/90.],$
  ;  [y1_0,y1_1,y1_1,y1_0,y1_0],color='dark green',/over
  x2_0=rpeak1[i]+(outsz-0.5)/60.
  y2_0=rpeak2[i]-(outsz-0.5)/60.
  x2_1=rpeak1[i]-(outsz-0.5)/60.
  y2_1=rpeak2[i]+(outsz-0.5)/60.
  cgarrow,x2_0,y2_0,x2_1,y2_1,/data,color='purple',/clip
  ;cgplot,[x2_0-1/90.,x2_1-1/90.,x2_1+1/90.,x2_0+1/90.,x2_0-1/90.],$
  ;  [y2_0-1/90.,y2_1-1/90.,y2_1+1/90.,y2_0+1/90.,y2_0-1/90.],color='purple',/over
  x3_0=rpeak1[i]+(outsz)/60.;x0+2*(x0-x1)
  y3_0=rpeak2[i]
  x3_1=rpeak1[i]-(outsz)/60.;x1+2*(x1-x0)
  y3_1=rpeak2[i]
  cgarrow,x3_0,y3_0,x3_1,y3_1,/data,color='brown',/clip
  ;cgplot,[x3_0,x3_1,x3_1,x3_0,x3_0],$
  ;  [y3_0-1/90.,y3_1-1/90.,y3_1+1/90.,y3_0+1/90.,y3_0-1/90.],color='brown',/over
  x4_0=rpeak1[i]+(outsz-0.5)/60.
  y4_0=rpeak2[i]+(outsz-0.5)/60.
  x4_1=rpeak1[i]-(outsz-0.5)/60.
  y4_1=rpeak2[i]-(outsz-0.5)/60.
  cgarrow,x4_0,y4_0,x4_1,y4_1,/data,color='blue',/clip
  ;cgplot,[x4_0+1/90.,x4_1+1/90.,x4_1-1/90.,x4_0-1/90.,x4_0+1/90.],$
  ;  [y4_0-1/90.,y4_1-1/90.,y4_1+1/90.,y4_0+1/90.,y4_0-1/90.],color='blue',/over
  cropfits,dathdr,vrange,dim='v',/dataform,output=pvout
  pv=mkpvbelt(pvout.dat,pvout.hdr,[x1_0,x1_1],[y1_0,y1_1],2,/gal,outdat=outdat)
  pvdata=outdat.dat
  pvhdr=outdat.hdr
  pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
  pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
  pvvrg=[min(pvv),max(pvv)]
  pvprg=[min(pvp),max(pvp)]*60.0-(outsz);*sqrt(2)   
  cgplot,[0],[0],/nodata,/noerase,position=pos2,xrange=pvvrg,yrange=pvprg,$
    ytickformat='(i)',xtickformat='(a1)',yminor=5
  cgloadct,50,clip=[64,240]
  pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
    round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
    *(0.05+indgen(10)*0.1)
  cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='dark green'
  cgplot,[0],[0],/nodata,/noerase,position=pos2,xrange=pvvrg,yrange=pvprg,$
    ytickinterval=2.0,ytitle="Position (')",yminor=5,xtitle='LSR Velocity (km s!U-1!N)'
  cgPlot, !x.crange,[0,0], LineStyle=1, /overplot
  cgPlot, [rw0[i],rw0[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot, [rw1[i],rw1[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
  pv=mkpvbelt(pvout.dat,pvout.hdr,[x2_0,x2_1],[y2_0,y2_1],2,/gal,outdat=outdat)
  pvdata=outdat.dat
  pvhdr=outdat.hdr
  pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
  pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
  pvvrg=[min(pvv),max(pvv)]
  pvprg=[min(pvp),max(pvp)]*60.0-(outsz-0.5)*sqrt(2)
  cgplot,[0],[0],/nodata,/noerase,position=pos3,xrange=pvvrg,yrange=pvprg,$
    ytickformat='(i)',xtickformat='(a1)',yminor=5,ytickinterval=2.0
  cgloadct,61,clip=[64,240]
  pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
    round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
    *(0.05+indgen(10)*0.1)
  cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='purple'
  cgplot,[0],[0],/nodata,/noerase,position=pos3,xrange=pvvrg,yrange=pvprg,$
    ytickinterval=2.0,ytitle="Position (')",yminor=5,xtickformat='(a1)'
  cgPlot, !x.crange,[0,0], LineStyle=1, /overplot
  cgPlot, [rw0[i],rw0[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot, [rw1[i],rw1[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
  pv=mkpvbelt(pvout.dat,pvout.hdr,[x3_0,x3_1],[y3_0,y3_1],2,/gal,outdat=outdat)
  pvdata=outdat.dat
  pvhdr=outdat.hdr
  pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
  pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
  pvvrg=[min(pvv),max(pvv)]
  pvprg=[min(pvp),max(pvp)]*60.0-(outsz);*sqrt(2)
  cgplot,[0],[0],/nodata,/noerase,position=pos4,xrange=pvvrg,yrange=pvprg,$
    ytickformat='(i)',xtickformat='(a1)',yminor=5,ytickinterval=2.0
  cgloadct,65,clip=[64,240]
  pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
    round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
    *(0.05+indgen(10)*0.1)
  cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='brown'
  cgplot,[0],[0],/nodata,/noerase,position=pos4,xrange=pvvrg,yrange=pvprg,$
    ytickinterval=2.0,ytitle="Position (')",yminor=5,xtickformat='(a1)'
  cgPlot, !x.crange,[0,0], LineStyle=1, /overplot
  cgPlot, [rw0[i],rw0[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot, [rw1[i],rw1[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
  pv=mkpvbelt(pvout.dat,pvout.hdr,[x4_0,x4_1],[y4_0,y4_1],2,/gal,outdat=outdat)
  pvdata=outdat.dat
  pvhdr=outdat.hdr
  pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
  pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
  pvvrg=[min(pvv),max(pvv)]
  pvprg=[min(pvp),max(pvp)]*60.0-(outsz-0.5)*sqrt(2)
  cgplot,[0],[0],/nodata,/noerase,position=pos5,xrange=pvvrg,yrange=pvprg,$
    ytickformat='(i)',xtickformat='(a1)',yminor=5,ytickinterval=2.0
  cgloadct,49,clip=[64,240]
  pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
    round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
    *(0.05+indgen(10)*0.1)
  cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='blue'
  cgplot,[0],[0],/nodata,/noerase,position=pos5,xrange=pvvrg,yrange=pvprg,$
    ytickinterval=2.0,ytitle="Position (')",yminor=5,xtickformat='(a1)'
  cgPlot, !x.crange,[0,0], LineStyle=1, /overplot
  cgPlot, [rw0[i],rw0[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot, [rw1[i],rw1[i]],!y.CRange, LineStyle=1, /overplot
  cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
  cgps_close
endfor
  
if ~file_test('redspectra') then spawn,'mkdir redspectra'
if ~file_test('bluespectra') then spawn,'mkdir bluespectra'
spawn,'rm ./*spectra/*.eps'
spawn,'mv redout*.eps redspectra'
spawn,'mv blueout*.eps bluespectra'



;fits_read,region+'_Ua_C.fits',datUa,hdrUa 
;fits_read,region+'_La_C.fits',datLa,hdrLa
;;fits_read,region+'_L2a_C.fits',datL2a,hdrL2a
;;readcol,'out_12CO.cat',autosign,out_c,lobesign,peak_n,peak_l,peak_b,peak_v,blue_l,blue_r,red_l,red_r,$
;;  format='A,I,A,I,F,F,F,F,F,F,F',stringskip='#',/silent
;
;readcol,'outflowcat.txt',num,glon,glat,bw0,bw1,rw0,rw1,format='I,F,F,F,F,F,F',stringskip='#',/silent
;
;dathdr=list(datUa,hdrUa)
;for i=0, n_elements(num)-1 do begin
;  
;    psname='out_'+num2str(i+1)+'.eps'
;    cgps_open,psname,font=!p.font,/quiet,default_thickness=1.0,charsize=1.0,/portrait
;    xsize=1000. & ysize=1000
;    cgDisplay, xsize=round(xsize), ysize=round(ysize)
;    pos0=[100./xsize,100./ysize,1.-50./xsize,1-50./ysize]
;    
;    cropfits,/dataform,dathdr,[glon[i]+3./60.,glon[i]-3./60.],[glat[i]-3./60.,glat[i]+3./60.],[bw0[i],bw1[i]],output=crop
;    bdata=smooth(crop.dat,[1,1,5],/edge_mirror)
;    biimap=total(bdata,3)*abs(sxpar(crop.hdr,'CDELT3'))/1000.
;    ;biimap[where(biimap lt max(biimap))] = 0.
;    biimap=congrid(biimap,5*n_elements(biimap[*,0]),5*n_elements(biimap[0,*]),cubic=0.5,/MINUS_ONE)
;    levels=max(biimap[round(n_elements(biimap[*,0])/2.-n_elements(biimap[*,0])/6.-1):round(n_elements(biimap[*,0])/2.+n_elements(biimap[*,0])/6.-1),$
;        round(n_elements(biimap[0,*])/2.-n_elements(biimap[0,*])/6.-1):round(n_elements(biimap[0,*])/2.+n_elements(biimap[0,*])/6.-1)])$
;        *(0.18+indgen(15)*0.1)
;    
;    crp1=sxpar(crop.hdr,'CRPIX1')
;    crv1=sxpar(crop.hdr,'CRVAL1')
;    del1=sxpar(crop.hdr,'CDELT1')
;    axs1=sxpar(crop.hdr,'NAXIS1')
;    l_l=(360.+(0.5-crp1)*del1+crv1) mod 360
;    l_r=(360.+(axs1+0.5-crp1)*del1+crv1) mod 360
;   
;    crp2=sxpar(crop.hdr,'CRPIX2')
;    crv2=sxpar(crop.hdr,'CRVAL2')
;    del2=sxpar(crop.hdr,'CDELT2')
;    axs2=sxpar(crop.hdr,'NAXIS2')
;    b_d=(0.5-crp2)*del2+crv2
;    b_u=(axs2+0.5-crp2)*del2+crv2
;
;    x_range=[l_l,l_r]
;    y_range=[b_d,b_u]
;    ratio=abs((b_u-b_d)/(l_l-l_r))
;    
;    cgplot,[0],[0],xrange=x_range,yrange=y_range,$
;      xtickinterval=0.02,ytickinterval=0.02,aspect=ratio,AxisColor='black',$
;      ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})'),position=pos0
;    ;cgplot,/over,peak_l[i],peak_b[i],psym=6,symsize=5,symcolor='dark green',thick=3
;    ;cgplot,/over,peak_l[i],peak_b[i],psym=6,symsize=5,symcolor='white',thick=1
;    xpos=!x.crange[0]+(!x.crange[1]-!x.crange[0])*0.15
;    ypos=!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.15
;    symsize=(52/3600.)*85./(x_range[1]-x_range[0])
;    ;cgplot,/over,xpos,ypos,psym=cgsymcat(16),symsize=symsize,color='black'
;    ;cgimage,imRGB,/noerase,position=pos0
;    cgloadct,1,/reverse,ncolors=15
;    cgcontour,biimap,/onimage,levels=levels,label=1,c_colors=indgen(15)
;
;    cropfits,/dataform,dathdr,[glon[i]+3./60.,glon[i]-3./60.],[glat[i]-3./60.,glat[i]+3./60.],[rw0[i],rw1[i]],output=crop
;    rdata=smooth(crop.dat,[1,1,5],/edge_mirror)
;    riimap=total(rdata,3)*abs(sxpar(crop.hdr,'CDELT3'))/1000.
;    ;riimap[where(riimap lt max(riimap))] = 0.
;    riimap=congrid(riimap,5*n_elements(riimap[*,0]),5*n_elements(riimap[0,*]),cubic=0.5,/MINUS_ONE)
;    levels=max(riimap[round(n_elements(riimap[*,0])/2.-n_elements(riimap[*,0])/6.-1):round(n_elements(riimap[*,0])/2.+n_elements(riimap[*,0])/6.-1),$
;        round(n_elements(riimap[0,*])/2.-n_elements(riimap[0,*])/6.-1):round(n_elements(riimap[0,*])/2.+n_elements(riimap[0,*])/6.-1)])$
;        *(0.18+indgen(15)*0.1)
;    cgloadct,3,/reverse,ncolors=15
;    cgcontour,riimap,/onimage,levels=levels,label=1,c_colors=indgen(15)
;    x0=189.804
;    y0=0.348
;    x1=189.779
;    y1=0.342
;    xx0=x0+2*(x0-x1)
;    yy0=y0+2*(y0-y1)
;    xx1=x1+2*(x1-x0)
;    yy1=y1+2*(y1-y0)
;    cgarrow,xx0,yy0,xx1,yy1,/data,color='dark green',/clip
;
;    cgps_close   
;    
;    psname1='pv_'+num2str(i+1)+'.eps'
;    cgps_open,psname1,font=!p.font,/quiet,default_thickness=1.0,charsize=1.0,/portrait
;    ;mkpvslice,'pvtempL_C.fits',/gal,l_pos,b_pos
;    ;fits_read,'pvslice.fits',pvdata,pvhdr
;    cropfits,dathdr,[0,17],dim='v',/dataform,output=pvout
;    pv=mkpvdat(pvout.dat,pvout.hdr,[xx0,xx1],[yy0,yy1],/gal,outdat=outdat)
;    pvdata=outdat.dat
;    pvhdr=outdat.hdr
;    pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
;    pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
;    pvvrg=[min(pvv),max(pvv)]
;    pvprg=[min(pvp),max(pvp)]*60.0
;    
;    xsize=1000. & ysize=1000.
;    cgDisplay, xsize=round(xsize), ysize=round(ysize)
;    p1=[100./xsize,100./ysize,1.-50./xsize,1-50./ysize]
;    
;    cgplot,[0],[0],/nodata,/noerase,position=p1,xrange=pvvrg,yrange=pvprg,$
;    ytickformat='(i)',xtickformat='(a1)',yminor=5
;    cgloadct,49,clip=[64,240]
;    pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
;        round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
;        *(0.025+indgen(10)*0.1)
;    cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='blue'
;    cgplot,[0],[0],/nodata,/noerase,position=p1,xrange=pvvrg,yrange=pvprg,$
;    ytickinterval=1.0,ytitle="Position (')",yminor=5,xtitle='LSR Velocity (km s!U-1!N)'
;    ;cgPlot, [v_c,v_c],[!y.CRange[0],!y.CRange[1]], LineStyle=1, Color='red', /overplot
;    ;cgPlot, [v[i],v[i]],[!y.CRange[0],!y.CRange[1]], LineStyle=1, Color='black', /overplot
;    
;    cgps_close
;  
;endfor

;if ~file_test('fig_out') then spawn,'mkdir fig_out'
;spawn,'rm ./fig_out/*.eps'
;spawn,'mv *.eps fig_out'

end