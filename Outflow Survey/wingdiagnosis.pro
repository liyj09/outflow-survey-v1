function coord2pix,hdr,value,dim,fitspix=fitspix
value=float(value)
dim=strcompress(string(dim),/remove_all)
if dim eq '1' && long(sxpar(hdr,'CRPIX1')) le 0l then pix=(value-360.-sxpar(hdr,'CRVAL'+dim))/sxpar(hdr,'CDELT'+dim)+sxpar(hdr,'CRPIX'+dim) $
 else pix=round((value-sxpar(hdr,'CRVAL'+dim))/sxpar(hdr,'CDELT'+dim)+sxpar(hdr,'CRPIX'+dim),/l64)
if keyword_set(fitspix) then return,round(pix,/l64) else return,round(pix,/l64)-1
end

function pix2coord,hdr,pixel,dim,fitspix=fitspix
if keyword_set(fitspix) then pixel=round(pixel,/l64) else pixel=round(pixel,/l64)+1l
dim=strcompress(string(dim),/remove_all)
return,(pixel-sxpar(hdr,'CRPIX'+dim))*sxpar(hdr,'CDELT'+dim)+sxpar(hdr,'CRVAL'+dim)
end

function checkwing,peak3,specU,vU,FWHMU,wingU,wingL,rmsU,rmsL
B_bad=0
R_bad=0
B_wing=([wingU[0],wingL[0]])
B_wing=B_wing[sort(B_wing)]
R_wing=([wingL[1],wingU[1]])
R_wing=R_wing[sort(R_wing)]
if wingL[0]-wingU[0] le 1.0 then B_bad++
if wingU[1]-wingL[1] le 1.0 then R_bad++
Fpos0=where(vU ge FWHMU[0]-0.085 and vU le FWHMU[0]+0.085,/null) & Fpos0=Fpos0[0]
Fpos1=where(vU ge FWHMU[1]-0.085 and vU le FWHMU[1]+0.085,/null) & Fpos1=Fpos1[0]
Bpos0=where(vU ge B_wing[0]-0.085 and vU le B_wing[0]+0.085,/null) & Bpos0=Bpos0[0]
Bpos1=where(vU ge B_wing[1]-0.085 and vU le B_wing[1]+0.085,/null) & Bpos1=Bpos1[0]
Rpos0=where(vU ge R_wing[0]-0.085 and vU le R_wing[0]+0.085,/null) & Rpos0=Rpos0[0]
Rpos1=where(vU ge R_wing[1]-0.085 and vU le R_wing[1]+0.085,/null) & Rpos1=Rpos1[0]
if Fpos0 gt Fpos1 then Fpos1=-1
iicore=total(specU[Fpos0:Fpos1])
iib=total(specU[Bpos0:Bpos1])
iir=total(specU[Rpos0:Rpos1])
Fspec=specU[Fpos0:Fpos1]
Bspec=specU[Bpos0:Bpos1]
Bspec=Bspec[where(Bspec gt 0.5*rmsL)]
if n_elements(Bspec) ge 5 then Bspec=smooth(Bspec[0:-2],3,/edge_mirror)
Rspec=specU[Rpos0:Rpos1]
Rspec=Rspec[where(Rspec gt 0.5*rmsL)]
if n_elements(Rspec) ge 5 then Rspec=smooth(Rspec[1:-1],3,/edge_mirror)
if Bspec[-1] lt 1 then B_bad++
if Rspec[0] lt 1 then R_bad++
Binter=interpol([Bspec[0],Bspec[-1]],n_elements(Bspec))
Rinter=interpol([Rspec[0],Rspec[-1]],n_elements(Rspec))
dB=Binter-Bspec
dR=Rinter-Rspec
if peak3-B_wing[1] gt 3*(FWHMU[1]-FWHMU[0]) then B_bad++
if R_wing[0]-peak3 gt 3*(FWHMU[1]-FWHMU[0]) then R_bad++
if (n_elements(where(dB lt -0.5,/null)) gt 1) then B_bad++
if (n_elements(where(dR lt -0.5,/null)) gt 1) then R_bad++
;if (iib/iicore gt 1.5) || (n_elements(where(dB lt -0.15,/null)) gt 1) then B_bad++
;if (iir/iicore gt 1.5) || (n_elements(where(dR lt -0.15,/null)) gt 1) then R_bad++
;if max(Fspec)*0.75 lt max(Bspec) then B_bad++
;if max(Fspec)*0.75 lt max(Rspec) then R_bad++
return,[B_bad,R_bad]
end

pro wingdiagnosis,region=region,Lvrange=Lvrange

time0=systime(/seconds)

if ~keyword_set(region) then region='region_C_III'
if ~keyword_set(Lvrange) then Lvrange=[-40,-20]
distance=2.0
L_lim=0.6
def_rmsU=0.25
def_rmsL=0.3
cd,'/home/lee/W3/'+region
outsz=(8./distance)>3
outsz=outsz<5
fits_read,'U_C.fits',datUa,hdrUa
fits_read,'L_C.fits',datLa,hdrLa
;fits_read,region+'_L2a_mask.fits',datL2a,hdrL2a
fits_read,'Lpeakv0.fits',peakvmap,vhdr
fits_read,'Ubvmap.fits',Ubvmap,bvhdr
fits_read,'Urvmap.fits',Urvmap,rvhdr
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

readcol,'redpeaks.cat', rsn,rpeak1,rpeak2,rscen1,rscen2,rsmaj,rsmin,rsposang,rsize1,rsize2,rpeak,$
  format='I,F,F,F,F,F,F,F,F,F,F',stringskip='#',/silent
readcol,'bluepeaks.cat',bsn,bpeak1,bpeak2,bscen1,bscen2,bsmaj,bsmin,bsposang,bsize1,bsize2,bpeak,$
  format='I,F,F,F,F,F,F,F,F,F,F',stringskip='#',/silent
  
openw,rout,/get_lun,'red_out.cat'
printf,rout,'# Rn','Glon','Glat','rw[0]','rw[1]',format='(a4,a9,a7,2a6)'
openw,bout,/get_lun,'blue_out.cat'
printf,bout,'# Bn','Glon','Glat','bw[0]','bw[1]',format='(a4,a9,a7,2a6)'

;redlobe
rcount=0
for i=0, n_elements(rsn)-1 do begin

  specU0=reform(datUa[coord2pix(hdrUa,rpeak1[i],1),coord2pix(hdrUa,rpeak2[i],2),*])
  specU1=reform(datUa[coord2pix(hdrUa,rpeak1[i],1)-1,coord2pix(hdrUa,rpeak2[i],2),*])
  specU2=reform(datUa[coord2pix(hdrUa,rpeak1[i],1)+1,coord2pix(hdrUa,rpeak2[i],2),*])
  specU3=reform(datUa[coord2pix(hdrUa,rpeak1[i],1),coord2pix(hdrUa,rpeak2[i],2)-1,*])
  specU4=reform(datUa[coord2pix(hdrUa,rpeak1[i],1),coord2pix(hdrUa,rpeak2[i],2)+1,*])
  spec_U=(specU0+specU1+specU2+specU3+specU4)/5.
  ;specU=gauss_smooth(spec_U,1,/edge_mirror)
  specU=spec_U
  
  specL0=reform(datLa[coord2pix(hdrLa,rpeak1[i],1),coord2pix(hdrLa,rpeak2[i],2),*])
  specL1=reform(datLa[coord2pix(hdrLa,rpeak1[i],1)-1,coord2pix(hdrLa,rpeak2[i],2),*])
  specL2=reform(datLa[coord2pix(hdrLa,rpeak1[i],1)+1,coord2pix(hdrLa,rpeak2[i],2),*])
  specL3=reform(datLa[coord2pix(hdrLa,rpeak1[i],1),coord2pix(hdrLa,rpeak2[i],2)-1,*])
  specL4=reform(datLa[coord2pix(hdrLa,rpeak1[i],1),coord2pix(hdrLa,rpeak2[i],2)+1,*])
  spec_L=(specL0+specL1+specL2+specL3+specL4)/5.
  ;specL=gauss_smooth(spec_L,1,/edge_mirror)
  specL=spec_L
  
  vLpeak=peakvmap[coord2pix(hdrLa,rpeak1[i],1),coord2pix(hdrLa,rpeak2[i],2)]
  peakposU=where(vU ge peakvmap[coord2pix(hdrUa,rpeak1[i],1),coord2pix(hdrUa,rpeak2[i],2)]-0.085 $
    and vU le peakvmap[coord2pix(hdrUa,rpeak1[i],1),coord2pix(hdrUa,rpeak2[i],2)]+0.085)
  peakposL=where(vL ge peakvmap[coord2pix(hdrLa,rpeak1[i],1),coord2pix(hdrLa,rpeak2[i],2)]-0.085 $
    and vL le peakvmap[coord2pix(hdrLa,rpeak1[i],1),coord2pix(hdrLa,rpeak2[i],2)]+0.085);,/null)
        
  if peakposU[0] ne -1 && peakposL[0] ne -1 then begin
    rmsU=def_rmsU
    rmsL=specL[peakposL[0]]/15.>def_rmsL
    FWHMU=vU[auto_range(specU,peakposU[0],specU[peakposU[0]]/1.8)]
    ;width=FWHMU[1]-FWHMU[0]
    ;vrange=[peakvmap[coord2pix(hdrUa,rpeak1[i],1),coord2pix(hdrUa,rpeak2[i],2)]-3.5*width,$
    ;  peakvmap[coord2pix(hdrUa,rpeak1[i],1),coord2pix(hdrUa,rpeak2[i],2)]+3.5*width]
    vrange=[Lvrange[1]-Ubvmap[coord2pix(bvhdr,rpeak1[i],1),coord2pix(bvhdr,rpeak2[i],2)]-3,$
      Urvmap[coord2pix(rvhdr,rpeak1[i],1),coord2pix(rvhdr,rpeak2[i],2)]+Lvrange[0]+3]
    wingU=vU[auto_range(specU,peakposU[0],rmsU)]
    wingL=vL[auto_range(specL,peakposL[0],rmsL)]
    baseline=make_array(n_elements(specU),value=0)
    check=checkwing(peakvmap[coord2pix(hdrLa,rpeak1[i],1),coord2pix(hdrLa,rpeak2[i],2)],$
      specU,vU,FWHMU,wingU,wingL,rmsU,rmsL)
    ;check=[0,0]
    if check[1] eq 0 && specL[peakposL[0]] gt L_lim then begin
      rcount++
      printf,rout,rcount,rpeak1[i],rpeak2[i],wingL[1],wingU[1],format='(i4,f9.3,f7.3,2f6.1)'
      ;printf,routm,rcount,rpeak1[i],rpeak2[i],wingL[1],wingU[1],format='(i4,f9.3,f7.3,2f6.1)'
      specpsname='redout_'+num2str(rcount)+'.eps'
      cgps_open,specpsname,font=!p.font,/quiet,default_thickness=1.0,charsize=0.7;,/encapsulated;,/portrait

      xsize1=1450. & ysize1=1000.
      cgDisplay, xsize=round(xsize1), ysize=round(ysize1)
      pos0=[100./xsize1,750./ysize1,700./xsize1,1-50./ysize1]
      pos1=[100./xsize1,100./ysize1,700./xsize1,700./ysize1]
      pos_r=(outsz-0.5)*sqrt(2)/outsz
      pos_d=850./(2*(pos_r+1))
      pos2=[800./xsize1,100./ysize1,1.-50./xsize1,(100.+pos_d)/ysize1]
      pos3=[800./xsize1,(100.+pos_d)/ysize1,1.-50./xsize1,(100.+pos_d*(1+pos_r))/ysize1]
      pos4=[800./xsize1,(100.+pos_d*(1+pos_r))/ysize1,1.-50./xsize1,(100.+pos_d*(2+pos_r))/ysize1]
      pos5=[800./xsize1,(100.+pos_d*(2+pos_r))/ysize1,1.-50./xsize1,950./ysize1]
      
      ;cgDisplay,xsize=400,ysize=200
      cgplot,vU,specU,psym=10,/nodata,xrange=vrange,color='blue',yrange=[0-max(spec_U)*0.15,max(spec_U)*1.15],$
        ytitle='T!DMB!N (K)',position=pos0,yminor=5,xtitle='LSR Velocity (km s!U-1!N)',$
        title='OUTFLOW REDLOBE CANDIDATE '+num2str(rcount)+': '+num2str(rpeak1[i],format='(f7.3)')+num2str(rpeak2[i],format='(f+6.3)')
      rpolyrg=where(vU ge wingL[1] and vU le wingU[1])
      rx1_arr=array_band(vU[rpolyrg],/band)
      rx2_arr=array_band(Reverse(vU[rpolyrg]),/band)
      ry1_arr=array_band(spec_U[rpolyrg])
      ry2_arr=array_band(baseline[rpolyrg])
      cgColorFill,[rx1_arr, rx2_arr, rx1_arr[0]],[ry1_arr,ry2_arr,ry1_arr[0]],Color='pink'
      cgplot,vU,spec_U,psym=10,color='blue',/overplot
      cgplot,vL,spec_L,psym=10,color='green',/overplot
      cgPlot, !x.CRange,[0.,0.], LineStyle=0, /overplot
      cgPlot,[vLpeak,vLpeak],[0,!y.crange[1]],linestyle=1,/over

    
      cropfits,/dataform,dathdr,[rpeak1[i]+outsz/60.,rpeak1[i]-outsz/60.],[rpeak2[i]-outsz/60.,rpeak2[i]+outsz/60.],[wingL[1],wingU[1]],output=crop
      rdata=smooth(crop.dat,[1,1,3],/edge_mirror)
      riimap=total(rdata,3)*abs(sxpar(crop.hdr,'CDELT3'))/1000.
      riimap=congrid(riimap,3*n_elements(riimap[*,0]),3*n_elements(riimap[0,*]),cubic=-0.5,/MINUS_ONE)
      ;biimap=smooth(biimap,[5,5],/edge_mirror)
      levels=max(riimap[round(n_elements(riimap[*,0])/2.-n_elements(riimap[*,0])/6.-1):round(n_elements(riimap[*,0])/2.+n_elements(riimap[*,0])/6.-1),$
        round(n_elements(riimap[0,*])/2.-n_elements(riimap[0,*])/6.-1):round(n_elements(riimap[0,*])/2.+n_elements(riimap[0,*])/6.-1)])$
        *(0.18+indgen(15)*0.1)
    
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
        xtickinterval=0.02,ytickinterval=0.02,aspect=ratio,AxisColor='black',$
        ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})'),position=pos1,/noerase
      xpos=!x.crange[0]+(!x.crange[1]-!x.crange[0])*0.1
      ypos=!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.1
      ;symsize=(52/3600.)*85./(x_range[1]-x_range[0])
      ;cgplot,/over,xpos,ypos,psym=cgsymcat(16),symsize=symsize,color='black'
      ;cgplot,/over,xpos,ypos,psym=cgsymcat(9,thick=3),symsize=symsize*0.95,color='white'
      ;cgloadct,3,ncolors=10,/reverse
      ;cgcolor
      cgloadct,62,ncolors=15
      cgcontour,riimap,/onimage,levels=levels,label=1,c_colors=indgen(15)
      
      x1_0=rpeak1[i]
      y1_0=rpeak2[i]-(outsz)/60
      x1_1=rpeak1[i]
      y1_1=rpeak2[i]+(outsz)/60
      cgarrow,x1_0,y1_0,x1_1,y1_1,/data,color='dark green',/clip
      cgplot,[x1_0-1/90.,x1_1-1/90.,x1_1+1/90.,x1_0+1/90.,x1_0-1/90.],$
        [y1_0,y1_1,y1_1,y1_0,y1_0],color='dark green',/over
      x2_0=rpeak1[i]+(outsz-0.5)/60.
      y2_0=rpeak2[i]-(outsz-0.5)/60.
      x2_1=rpeak1[i]-(outsz-0.5)/60.
      y2_1=rpeak2[i]+(outsz-0.5)/60.
      cgarrow,x2_0,y2_0,x2_1,y2_1,/data,color='purple',/clip
      cgplot,[x2_0-1/90.,x2_1-1/90.,x2_1+1/90.,x2_0+1/90.,x2_0-1/90.],$
        [y2_0-1/90.,y2_1-1/90.,y2_1+1/90.,y2_0+1/90.,y2_0-1/90.],color='purple',/over
      x3_0=rpeak1[i]+(outsz)/60.;x0+2*(x0-x1)
      y3_0=rpeak2[i]
      x3_1=rpeak1[i]-(outsz)/60.;x1+2*(x1-x0)
      y3_1=rpeak2[i]
      cgarrow,x3_0,y3_0,x3_1,y3_1,/data,color='brown',/clip
      cgplot,[x3_0,x3_1,x3_1,x3_0,x3_0],$
        [y3_0-1/90.,y3_1-1/90.,y3_1+1/90.,y3_0+1/90.,y3_0-1/90.],color='brown',/over
      x4_0=rpeak1[i]+(outsz-0.5)/60.
      y4_0=rpeak2[i]+(outsz-0.5)/60.
      x4_1=rpeak1[i]-(outsz-0.5)/60.
      y4_1=rpeak2[i]-(outsz-0.5)/60.
      cgarrow,x4_0,y4_0,x4_1,y4_1,/data,color='blue',/clip
      cgplot,[x4_0+1/90.,x4_1+1/90.,x4_1-1/90.,x4_0-1/90.,x4_0+1/90.],$
        [y4_0-1/90.,y4_1-1/90.,y4_1+1/90.,y4_0+1/90.,y4_0-1/90.],color='blue',/over

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
      cgPlot, [wingL[1],wingL[1]],!y.CRange, LineStyle=1, /overplot
      cgPlot, [wingU[1],wingU[1]],!y.CRange, LineStyle=1, /overplot
      cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
;      cgColorFill,[wingL[1], wingU[1], wingU[1], wingL[1], wingL[1]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=45
;      cgColorFill,[wingL[1], wingU[1], wingU[1], wingL[1], wingL[1]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=-45
     
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
      cgPlot, [wingL[1],wingL[1]],!y.CRange, LineStyle=1, /overplot
      cgPlot, [wingU[1],wingU[1]],!y.CRange, LineStyle=1, /overplot
      cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
;      cgColorFill,[wingL[1], wingU[1], wingU[1], wingL[1], wingL[1]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=45
;      cgColorFill,[wingL[1], wingU[1], wingU[1], wingL[1], wingL[1]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=-45
        
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
      cgPlot, [wingL[1],wingL[1]],!y.CRange, LineStyle=1, /overplot
      cgPlot, [wingU[1],wingU[1]],!y.CRange, LineStyle=1, /overplot
      cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
;      cgColorFill,[wingL[1], wingU[1], wingU[1], wingL[1], wingL[1]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=45
;      cgColorFill,[wingL[1], wingU[1], wingU[1], wingL[1], wingL[1]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=-45
        
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
      cgPlot, [wingL[1],wingL[1]],!y.CRange, LineStyle=1, /overplot
      cgPlot, [wingU[1],wingU[1]],!y.CRange, LineStyle=1, /overplot
      cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
;      cgColorFill,[wingL[1], wingU[1], wingU[1], wingL[1], wingL[1]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=45
;      cgColorFill,[wingL[1], wingU[1], wingU[1], wingL[1], wingL[1]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=-45
     
      
      cgps_close
      
      
    endif
  endif
  
endfor

;bluelobe
bcount=0
for i=0, n_elements(bsn)-1 do begin
  
  specU0=reform(datUa[coord2pix(hdrUa,bpeak1[i],1),coord2pix(hdrUa,bpeak2[i],2),*])
  specU1=reform(datUa[coord2pix(hdrUa,bpeak1[i],1)-1,coord2pix(hdrUa,bpeak2[i],2),*])
  specU2=reform(datUa[coord2pix(hdrUa,bpeak1[i],1)+1,coord2pix(hdrUa,bpeak2[i],2),*])
  specU3=reform(datUa[coord2pix(hdrUa,bpeak1[i],1),coord2pix(hdrUa,bpeak2[i],2)-1,*])
  specU4=reform(datUa[coord2pix(hdrUa,bpeak1[i],1),coord2pix(hdrUa,bpeak2[i],2)+1,*])
  spec_U=(specU0+specU1+specU2+specU3+specU4)/5.
  ;specU=gauss_smooth(spec_U,1,/edge_mirror)
  specU=spec_U
  
  specL0=reform(datLa[coord2pix(hdrLa,bpeak1[i],1),coord2pix(hdrLa,bpeak2[i],2),*])
  specL1=reform(datLa[coord2pix(hdrLa,bpeak1[i],1)-1,coord2pix(hdrLa,bpeak2[i],2),*])
  specL2=reform(datLa[coord2pix(hdrLa,bpeak1[i],1)+1,coord2pix(hdrLa,bpeak2[i],2),*])
  specL3=reform(datLa[coord2pix(hdrLa,bpeak1[i],1),coord2pix(hdrLa,bpeak2[i],2)-1,*])
  specL4=reform(datLa[coord2pix(hdrLa,bpeak1[i],1),coord2pix(hdrLa,bpeak2[i],2)+1,*])
  spec_L=(specL0+specL1+specL2+specL3+specL4)/5.
  ;specL=gauss_smooth(spec_L,1,/edge_mirror)
  specL=spec_L
   
  vLpeak=peakvmap[coord2pix(hdrLa,bpeak1[i],1),coord2pix(hdrLa,bpeak2[i],2)]   
  peakposU=where(vU ge peakvmap[coord2pix(hdrUa,bpeak1[i],1),coord2pix(hdrUa,bpeak2[i],2)]-0.085 $
    and vU le peakvmap[coord2pix(hdrUa,bpeak1[i],1),coord2pix(hdrUa,bpeak2[i],2)]+0.085)
  peakposL=where(vL ge peakvmap[coord2pix(hdrLa,bpeak1[i],1),coord2pix(hdrLa,bpeak2[i],2)]-0.085 $
    and vL le peakvmap[coord2pix(hdrLa,bpeak1[i],1),coord2pix(hdrLa,bpeak2[i],2)]+0.085);,/null)
      
  if peakposU[0] ne -1 && peakposL[0] ne -1 then begin
    rmsU=def_rmsU
    rmsL=specL[peakposL[0]]/15.>def_rmsL
    FWHMU=vU[auto_range(specU,peakposU[0],specU[peakposU[0]]/1.8)]
    ;width=FWHMU[1]-FWHMU[0]
    ;vrange=[peakvmap[coord2pix(hdrUa,bpeak1[i],1),coord2pix(hdrUa,bpeak2[i],2)]-3.5*width,$
    ;  peakvmap[coord2pix(hdrUa,bpeak1[i],1),coord2pix(hdrUa,bpeak2[i],2)]+3.5*width]
    vrange=[Lvrange[1]-Ubvmap[coord2pix(bvhdr,bpeak1[i],1),coord2pix(bvhdr,bpeak2[i],2)]-5,$
      Urvmap[coord2pix(rvhdr,bpeak1[i],1),coord2pix(rvhdr,bpeak2[i],2)]+Lvrange[0]+5]
    wingU=vU[auto_range(specU,peakposU[0],rmsU)]
    wingL=vL[auto_range(specL,peakposL[0],rmsL)]
    baseline=make_array(n_elements(specU),value=0)
    check=checkwing(peakvmap[coord2pix(hdrLa,bpeak1[i],1),coord2pix(hdrLa,bpeak2[i],2)],$
      specU,vU,FWHMU,wingU,wingL,rmsU,rmsL)
    if check[0] eq 0 && specL[peakposL[0]] gt L_lim then begin
      bcount++
      printf,bout,bcount,bpeak1[i],bpeak2[i],wingU[0],wingL[0],format='(i4,f9.3,f7.3,2f6.1)'
      ;printf,boutm,bcount,bpeak1[i],bpeak2[i],wingU[0],wingL[0],format='(i4,f9.3,f7.3,2f6.1)'
      specpsname='blueout_'+num2str(bcount)+'.eps'
      cgps_open,specpsname,font=!p.font,/quiet,default_thickness=1.0,charsize=0.6;,/encapsulated;,/portrait
      ;cgDisplay,xsize=400,ysize=200
      xsize1=1450. & ysize1=1000.
      cgDisplay, xsize=round(xsize1), ysize=round(ysize1)
      pos0=[100./xsize1,750./ysize1,700./xsize1,1-50./ysize1]
      pos1=[100./xsize1,100./ysize1,700./xsize1,700./ysize1]
      pos_r=(outsz-0.5)*sqrt(2)/outsz
      pos_d=850./(2*(pos_r+1))
      pos2=[800./xsize1,100./ysize1,1.-50./xsize1,(100.+pos_d)/ysize1]
      pos3=[800./xsize1,(100.+pos_d)/ysize1,1.-50./xsize1,(100.+pos_d*(1+pos_r))/ysize1]
      pos4=[800./xsize1,(100.+pos_d*(1+pos_r))/ysize1,1.-50./xsize1,(100.+pos_d*(2+pos_r))/ysize1]
      pos5=[800./xsize1,(100.+pos_d*(2+pos_r))/ysize1,1.-50./xsize1,950./ysize1]
      
      cgplot,vU,specU,psym=10,/nodata,xrange=vrange,color='blue',yrange=[0-max(spec_U)*0.15,max(spec_U)*1.15],$
        ytitle='T!DMB!N (K)',position=pos0,yminor=5,xtitle='LSR Velocity (km s!U-1!N)',$
        title='OUTFLOW BLUELOBE CANDIDATE '+num2str(bcount)+': '+num2str(bpeak1[i],format='(f7.3)')+num2str(bpeak2[i],format='(f+6.3)')
      bpolyrg=where(vU ge wingU[0] and vU le wingL[0])
      bx1_arr=array_band(vU[bpolyrg],/band)
      bx2_arr=array_band(Reverse(vU[bpolyrg]),/band)
      by1_arr=array_band(spec_U[bpolyrg])
      by2_arr=array_band(baseline[bpolyrg])
      cgColorFill,[bx1_arr, bx2_arr, bx1_arr[0]],[by1_arr,by2_arr,by1_arr[0]],Color='Sky Blue'
      cgplot,vU,spec_U,psym=10,color='blue',/overplot
      cgplot,vL,spec_L,psym=10,color='green',/overplot
      cgPlot, !x.CRange,[0.,0.], LineStyle=0, Color='black',/overplot
      cgPlot,[vLpeak,vLpeak],[0,!y.crange[1]],linestyle=1,/over
    
      cropfits,/dataform,dathdr,[bpeak1[i]+outsz/60.,bpeak1[i]-outsz/60.],[bpeak2[i]-outsz/60.,bpeak2[i]+outsz/60.],[wingU[0],wingL[0]],output=crop
      bdata=smooth(crop.dat,[1,1,3],/edge_mirror)
      biimap=total(bdata,3)*abs(sxpar(crop.hdr,'CDELT3'))/1000.
      biimap=congrid(biimap,3*n_elements(biimap[*,0]),3*n_elements(biimap[0,*]),cubic=-0.5,/MINUS_ONE)
      ;biimap=smooth(biimap,[5,5],/edge_mirror)
      levels=max(biimap[round(n_elements(biimap[*,0])/2.-n_elements(biimap[*,0])/6.-1):round(n_elements(biimap[*,0])/2.+n_elements(biimap[*,0])/6.-1),$
        round(n_elements(biimap[0,*])/2.-n_elements(biimap[0,*])/6.-1):round(n_elements(biimap[0,*])/2.+n_elements(biimap[0,*])/6.-1)])$
        *(0.18+indgen(15)*0.1)
    
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
        xtickinterval=0.02,ytickinterval=0.02,aspect=ratio,AxisColor='black',$
        ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})'),position=pos1,/noerase
      xpos=!x.crange[0]+(!x.crange[1]-!x.crange[0])*0.10
      ypos=!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.10
      ;symsize=(52/3600.)*85./(x_range[1]-x_range[0])
      ;cgplot,/over,xpos,ypos,psym=cgsymcat(16),symsize=symsize,color='black'
      ;cgplot,/over,xpos,ypos,psym=cgsymcat(9,thick=3),symsize=symsize*0.95,color='white'
      ;cgloadct,3,ncolors=10,/reverse
      cgloadct,49,ncolors=15
      cgcontour,biimap,/onimage,levels=levels,label=1,c_colors=indgen(15)
      
;      xx0=bpeak1[i]+(outsz-0.5)/60.;x0+2*(x0-x1)
;      yy0=bpeak2[i]-(outsz-0.5)/60;y0+2*(y0-y1)
;      xx1=bpeak1[i]-(outsz-0.5)/60.;x1+2*(x1-x0)
;      yy1=bpeak2[i]+(outsz-0.5)/60;y1+2*(y1-y0)
;      cgarrow,xx0,yy0,xx1,yy1,/data,color='dark green',/clip
;      cgplot,[xx0-1/90.,xx1-1/90.,xx1+1/90.,xx0+1/90.,xx0-1/90.],$
;        [yy0-1/90.,yy1-1/90.,yy1+1/90.,yy0+1/90.,yy0-1/90.],color='dark green',/over
;
;      xx00=bpeak1[i]+(outsz-0.5)/60.;x0+2*(x0-x1)
;      yy00=bpeak2[i]+(outsz-0.5)/60;y0+2*(y0-y1)
;      xx11=bpeak1[i]-(outsz-0.5)/60.;x1+2*(x1-x0)
;      yy11=bpeak2[i]-(outsz-0.5)/60;y1+2*(y1-y0)
;      cgarrow,xx00,yy00,xx11,yy11,/data,color='purple';,/clip
;      cgplot,[xx00+1/90.,xx11+1/90.,xx11-1/90.,xx00-1/90.,xx00+1/90.],$
;        [yy00-1/90.,yy11-1/90.,yy11+1/90.,yy00+1/90.,yy00-1/90.],color='purple',/over

      x1_0=bpeak1[i]
      y1_0=bpeak2[i]-(outsz)/60
      x1_1=bpeak1[i]
      y1_1=bpeak2[i]+(outsz)/60
      cgarrow,x1_0,y1_0,x1_1,y1_1,/data,color='dark green',/clip
      cgplot,[x1_0-1/90.,x1_1-1/90.,x1_1+1/90.,x1_0+1/90.,x1_0-1/90.],$
        [y1_0,y1_1,y1_1,y1_0,y1_0],color='dark green',/over
      x2_0=bpeak1[i]+(outsz-0.5)/60.
      y2_0=bpeak2[i]-(outsz-0.5)/60.
      x2_1=bpeak1[i]-(outsz-0.5)/60.
      y2_1=bpeak2[i]+(outsz-0.5)/60.
      cgarrow,x2_0,y2_0,x2_1,y2_1,/data,color='purple',/clip
      cgplot,[x2_0-1/90.,x2_1-1/90.,x2_1+1/90.,x2_0+1/90.,x2_0-1/90.],$
        [y2_0-1/90.,y2_1-1/90.,y2_1+1/90.,y2_0+1/90.,y2_0-1/90.],color='purple',/over
      x3_0=bpeak1[i]+(outsz)/60.;x0+2*(x0-x1)
      y3_0=bpeak2[i]
      x3_1=bpeak1[i]-(outsz)/60.;x1+2*(x1-x0)
      y3_1=bpeak2[i]
      cgarrow,x3_0,y3_0,x3_1,y3_1,/data,color='brown',/clip
      cgplot,[x3_0,x3_1,x3_1,x3_0,x3_0],$
        [y3_0-1/90.,y3_1-1/90.,y3_1+1/90.,y3_0+1/90.,y3_0-1/90.],color='brown',/over
      x4_0=bpeak1[i]+(outsz-0.5)/60.
      y4_0=bpeak2[i]+(outsz-0.5)/60.
      x4_1=bpeak1[i]-(outsz-0.5)/60.
      y4_1=bpeak2[i]-(outsz-0.5)/60.
      cgarrow,x4_0,y4_0,x4_1,y4_1,/data,color='blue',/clip
      cgplot,[x4_0+1/90.,x4_1+1/90.,x4_1-1/90.,x4_0-1/90.,x4_0+1/90.],$
        [y4_0-1/90.,y4_1-1/90.,y4_1+1/90.,y4_0+1/90.,y4_0-1/90.],color='blue',/over

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
      cgPlot, [wingL[0],wingL[0]],!y.CRange, LineStyle=1, /overplot
      cgPlot, [wingU[0],wingU[0]],!y.CRange, LineStyle=1, /overplot
      cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
;      cgColorFill,[wingU[0], wingL[0], wingL[0], wingU[0], wingU[0]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=45
;      cgColorFill,[wingU[0], wingL[0], wingL[0], wingU[0], wingU[0]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=-45
     
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
      cgPlot, [wingL[0],wingL[0]],!y.CRange, LineStyle=1, /overplot
      cgPlot, [wingU[0],wingU[0]],!y.CRange, LineStyle=1, /overplot
      cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
;      cgColorFill,[wingU[0], wingL[0], wingL[0], wingU[0], wingU[0]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=45
;      cgColorFill,[wingU[0], wingL[0], wingL[0], wingU[0], wingU[0]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=-45
        
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
      cgPlot, [wingL[0],wingL[0]],!y.CRange, LineStyle=1, /overplot
      cgPlot, [wingU[0],wingU[0]],!y.CRange, LineStyle=1, /overplot
      cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
;      cgColorFill,[wingU[0], wingL[0], wingL[0], wingU[0], wingU[0]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=45
;      cgColorFill,[wingU[0], wingL[0], wingL[0], wingU[0], wingU[0]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=-45
        
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
      cgPlot, [wingL[0],wingL[0]],!y.CRange, LineStyle=1, /overplot
      cgPlot, [wingU[0],wingU[0]],!y.CRange, LineStyle=1, /overplot
      cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
;      cgColorFill,[wingU[0], wingL[0], wingL[0], wingU[0], wingU[0]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=45
;      cgColorFill,[wingU[0], wingL[0], wingL[0], wingU[0], wingU[0]],$
;        [!y.crange[0],!y.crange,!y.crange[1],!y.crange[0]],Color='black',/line_fill,ORIENTATION=-45
      
      cgps_close
      
    endif
  endif
  
endfor

free_lun,rout,bout;,routm,boutm
  
if ~file_test('redspectra') then spawn,'mkdir redspectra'
;if ~file_test('redcontours') then spawn,'mkdir redcontours'
if ~file_test('bluespectra') then spawn,'mkdir bluespectra'
;if ~file_test('bluecontours') then spawn,'mkdir bluecontours'
;spawn,'rm ./*spectra/*.eps'
;;spawn,'rm ./*contours/*.eps'
;spawn,'mv redout*.eps redspectra'
;;spawn,'mv redcon*.eps redcontours'
;spawn,'mv blueout*.eps bluespectra'
;;spawn,'mv bluecon*.eps bluecontours'

;crp1=sxpar(hdrUa,'CRPIX1')
;crv1=sxpar(hdrUa,'CRVAL1')
;del1=sxpar(hdrUa,'CDELT1')
;l_l=(360.+(0.5-crp1)*del1+crv1) mod 360
;l_r=(360.+(sxpar(hdrUa,'NAXIS1')+0.5-crp1)*del1+crv1) mod 360
;   
;crp2=sxpar(hdrUa,'CRPIX2')
;crv2=sxpar(hdrUa,'CRVAL2')
;del2=sxpar(hdrUa,'CDELT2')
;b_d=(0.5-crp2)*del2+crv2
;b_u=(sxpar(hdrUa,'NAXIS2')+0.5-crp2)*del2+crv2
;
;x_range=[l_l,l_r]
;y_range=[b_d,b_u]
;ratio=abs((b_u-b_d)/(l_l-l_r))
;
;;rcount,rpeak1[i],rpeak2[i],rcen1[i],rcen2[i],rmajor[i]<5./60,rminor[i],rposangle[i],rsize1[i],rsize2[i],rpeak[i]
;readcol,'redpeaks.cat',rsn,rpeak1,rpeak2,rscen1,rscen2,rsmaj,rsmin,rsposang,rsize1,rsize2,rpeak,$
;  format='I,F,F,F,F,F,F,F,F,F,F',stringskip='#',/silent
;readcol,'bluepeaks.cat',bsn,bpeak1,bpeak2,bscen1,bscen2,bsmaj,bsmin,bsposang,bsize1,bsize2,bpeak,$
;  format='I,F,F,F,F,F,F,F,F,F,F',stringskip='#',/silent
;
;psname='r_out.eps'
;cgps_open,psname,font=!p.font,/quiet,default_thickness=1.0,charsize=1.0;,/encapsulated;,/portrait
;xsize=1000. & ysize=1000.*ratio
;cgDisplay, xsize=round(xsize), ysize=round(ysize)
;pos0=[100./xsize,100.*ratio/ysize,1.-50./xsize,1-50./ysize]
;cgplot,[0],[0],xrange=x_range,yrange=y_range,$
;  xtickinterval=0.5,ytickinterval=0.5,aspect=ratio,AxisColor='black',xthick=5,ythick=5,$
;  ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})'),position=pos0
;cgloadct,0,/reverse;,ncolors=10
;
;cropfits,dathdrL,Lvrange,dim='v',/dataform,output=crL
;Liimap=total(smooth(crL.dat,[1,1,3]),3)
;cgimage,bytscl(Liimap<max(liimap)/1.5),/overplot;,/noerase,position=pos0
;
;;cropfits,dathdrL2,Lvrange,dim='v',/dataform,output=crL2
;;L2iimap=total(crL2.dat,3)
;;L2iimap=smooth(L2iimap,[3,3],/edge_mirror)
;;levels=max(L2iimap)*(0.15+indgen(10)*0.1)
;;cgcontour,/onimage,L2iimap,label=0,color='red',level=levels
;;levels=levels,label=1,color='blue'
;
;cgplot,[0],[0],xrange=x_range,yrange=y_range,position=pos0,$
;  xtickinterval=0.5,ytickinterval=0.5,aspect=ratio,$
;  AxisColor='black',/noerase,xtickformat='(a1)',ytickformat='(a1)'
;for i=0, n_elements(rsn)-1 do begin
;  ;cgPlotS,draw_ellipse(rscen1[i],rscen2[i],rsmaj[i],rsmin[i],pa=-90-rsposang[i]), Color='white',NOCLIP=0,thick=2
;  ;cgplotS,draw_ellipse(rscen1[i],rscen2[i],rsmaj[i],rsmin[i],pa=-90-rsposang[i]), Color='red',NOCLIP=0,thick=1
;  ;cgplot,/over,rpeak1[i],rpeak2[i],psym=7,symsize=0.75,symcolor='white',thick=3
;  ;cgplot,/over,rpeak1[i],rpeak2[i],psym=7,symsize=0.75,symcolor='red',thick=3
;  ;cgsymcat
;endfor
;for i=0, n_elements(bsn)-1 do begin
;  ;cgPlotS,draw_ellipse(bscen1[i],bscen2[i],bsmaj[i],bsmin[i],pa=-90-bsposang[i]), Color='white',NOCLIP=0,thick=2
;  ;cgplotS,draw_ellipse(bscen1[i],bscen2[i],bsmaj[i],bsmin[i],pa=-90-bsposang[i]), Color='blue',NOCLIP=0,thick=1
;  ;cgplot,/over,bpeak1[i],bpeak2[i],psym=1,symsize=0.75,symcolor='white',thick=1
;  ;cgplot,/over,bpeak1[i],bpeak2[i],psym=1,symsize=0.75,symcolor='blue',thick=3
;endfor
;readcol,'red_out.cat',nnum,npeak1,npeak2,nwing1,nwing2,format='I,F,F,F,F',/silent
;for i=0, n_elements(npeak1)-1 do begin
;  ;cgplot,/over,npeak1[i],npeak2[i],psym=9,symsize=0.75,symcolor='white',thick=3
;  cgplot,/over,npeak1[i],npeak2[i],psym=9,symsize=0.75,symcolor='red',thick=3
;  ;cgplot,/over,npeak1[i],npeak2[i],psym=7,symsize=0.5,symcolor='white',thick=3
;  ;cgplot,/over,npeak1[i],npeak2[i],psym=7,symsize=0.5,symcolor='red',thick=1
;  cgtext,npeak1[i]-0.03,npeak2[i]+0.015,color='white',/data,num2str(nnum[i]),charsize=0.8,alignment=0.5,charthick=3
;  cgtext,npeak1[i]-0.03,npeak2[i]+0.015,color='red',/data,num2str(nnum[i]),charsize=0.8,alignment=0.5,charthick=1
;endfor
;readcol,'blue_out.cat',bnnum,bnpeak1,bnpeak2,bnwing1,bnwing2,format='I,F,F,F,F',/silent
;for i=0, n_elements(bnpeak1)-1 do begin
;  ;cgplot,/over,bnpeak1[i],bnpeak2[i],psym=9,symsize=0.75,symcolor='white',thick=3
;  cgplot,/over,bnpeak1[i],bnpeak2[i],psym=9,symsize=0.75,symcolor='blue',thick=3
;  ;cgplot,/over,bnpeak1[i],bnpeak2[i],psym=1,symsize=0.75,symcolor='white',thick=3
;  ;cgplot,/over,bnpeak1[i],bnpeak2[i],psym=1,symsize=0.75,symcolor='blue',thick=1
;  cgtext,bnpeak1[i]+0.03,bnpeak2[i]+0.015,color='white',/data,num2str(bnnum[i]),charsize=0.8,alignment=0.5,charthick=3
;  cgtext,bnpeak1[i]+0.03,bnpeak2[i]+0.015,color='blue',/data,num2str(bnnum[i]),charsize=0.8,alignment=0.5,charthick=1
;endfor
;cgps_close

time1=systime(/seconds)

print,' TIME COMSUMPTION: '+num2str(time1-time0,format='(i)')+' SECONDS.'
end
