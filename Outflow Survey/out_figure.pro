function wise_scl,dat
sz=size(dat)
datc=dat[round(sz[1]/2.-sz[1]/3-1):round(sz[1]/2.+sz[1]/3-1),$
  round(sz[2]/2.-sz[2]/3-1):round(sz[2]/2.+sz[2]/3-1)]
;datc=dat
min=median(datc)-(max(datc)-median(datc))*0.2
;max=max(datc)-(max(datc)-median(datc))*0.4
;min=median(datc);-max(datc)*0.2
max=max(datc)*0.5+median(datc)*0.505
re_dat=reform(logscl(dat,min=min,max=max),1,sz[1],sz[2])
return, re_dat
end

pro out_figure

;cd,'/media/alpha/W3_supplement/'
cd,'/media/alpha/W3/'
;cd,'/media/lee/data1/W3/'
;cd,'/media/dell/data1/W3_supplement/'

;if ~keyword_set(region) then region='region_B'
;if ~keyword_set(Lvrange) then Lvrange=[-85,-55]
;distance=0.599
;if ~keyword_set(region) then region='region_C_I'
;if ~keyword_set(Lvrange) then Lvrange=[-70,-30]
;distance=2
;if ~keyword_set(region) then region='region_C_II'
;if ~keyword_set(Lvrange) then Lvrange=[-40,-34]
;distance=2
if ~keyword_set(region) then region='region_C_III'
if ~keyword_set(Lvrange) then Lvrange=[-34,-28]
distance=2
;if ~keyword_set(region) then region='region_D'
;if ~keyword_set(Lvrange) then Lvrange=[-29,-20]
;distance=1.3
;if ~keyword_set(region) then region='region_E_I'
;if ~keyword_set(Lvrange) then Lvrange=[-25,-5]
;distance=0.6
;if ~keyword_set(region) then region='region_E_II'
;if ~keyword_set(Lvrange) then Lvrange=[-20,0]
;distance=0.6
;if ~keyword_set(region) then region='region_F'
;if ~keyword_set(Lvrange) then Lvrange=[-15,10]
;distance=0.6
;if ~keyword_set(region) then region='part1'
;if ~keyword_set(Lvrange) then Lvrange=[-60,-20]
;distance=2.0
;if ~keyword_set(region) then region='part2'
;if ~keyword_set(Lvrange) then Lvrange=[-25,10]
;distance=0.8
;if ~keyword_set(region) then region='part3'
;if ~keyword_set(Lvrange) then Lvrange=[-60,0]
;distance=2.0
;if ~keyword_set(region) then region='part5'
;if ~keyword_set(Lvrange) then Lvrange=[-25,10]
;distance=0.2
;if ~keyword_set(region) then region='part6'
;if ~keyword_set(Lvrange) then Lvrange=[-85,-50]
;distance=6.3
;if ~keyword_set(region) then region='part7'
;if ~keyword_set(Lvrange) then Lvrange=[-70,-35]
;distance=2.0
;if ~keyword_set(region) then region='part8'
;if ~keyword_set(Lvrange) then Lvrange=[-60,-25]
;distance=2.0
;if ~keyword_set(region) then region='part9'
;if ~keyword_set(Lvrange) then Lvrange=[-55,-25]
;distance=2.0
;if ~keyword_set(region) then region='part10'
;if ~keyword_set(Lvrange) then Lvrange=[-25,-5]
;distance=0.8
;if ~keyword_set(region) then region='part11'
;if ~keyword_set(Lvrange) then Lvrange=[-25,-5]
;distance=0.8
;if ~keyword_set(region) then region='part12'
;if ~keyword_set(Lvrange) then Lvrange=[-15,10]
;distance=0.8
;if ~keyword_set(region) then region='part13'
;if ~keyword_set(Lvrange) then Lvrange=[-15,10]
;distance=0.8

;cd,'/media/alpha/W3_supplement/'+region+'/candidates'
cd,'/media/alpha/W3/'+region+'/candidates'
;cd,'/media/lee/data1/W3/'+region+'/candidates'
;cd,'/media/dell/data1/W3_supplement/'+region+'/candidates'
L_lim=0.6
def_rmsU=0.25
def_rmsL=0.3
outsz=(6./distance)>3
outsz=outsz<5
;fits_read,'../'+region+'_U_C.fits',datUa,hdrUa
;fits_read,'../'+region+'_L_C.fits',datLa,hdrLa
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
baseline=make_array(n_elements(specU),value=0)

readcol,'blue_out.cat',bsn,bpeak1,bpeak2,bw0,bw1,format='I,F,F,F,F',stringskip='#'
readcol,'red_out.cat',rsn,rpeak1,rpeak2,rw0,rw1,format='I,F,F,F,F',stringskip='#'
readcol,'outflowcat.cat',num,glon,glat,bluev0,bluev1,redv0,redv1,tag,bnum,rnum,index,format='I,F,F,F,F,F,F,A,I,I,I',stringskip='#'
readcol,'con_start.cat',blow,bmax,rlow,rmax,format='F,F,F,F',stringskip='#'
readcol,'distance.cat',id,distance1,rank,format='I,F,A',stringskip='#'
;readcol,'vl_range.cat',vvv1,vvv2,iddd,format='F,F,I' ;only for C_I
readcol,'vl_range.cat',vvv1,vvv2,format='F,F'

xsize1=1450 & ysize1=1200
p0=[100./xsize1,662.5/ysize1,800./xsize1,1-25./ysize1]
p1=[100./xsize1,75./ysize1,800./xsize1,587.5/ysize1]
p2=[900./xsize1,700./ysize1,1.-50./xsize1,1-75./ysize1]
p3=[900./xsize1,125./ysize1,1.-50./xsize1,550./ysize1]

;cgcolorbar
for i=0, n_elements(num)-1 do begin
 ; if i ne 46 then continue
 ; print,glon,glat
;  print,vvv1
  if (index[i] eq 204) or (index[i] eq 317) then index[i]=0
  if index[i] gt 204 then index[i]=index[i]-1
  if index[i] gt 316 then index[i]=index[i]-1
  distance=distance1[i]/1000.d
  outsz=(6./distance)>3
  outsz=outsz<4.5
  if tag[i] eq 'D' then begin
    j=where(bsn eq bnum[i])
    k=where(rsn eq rnum[i])
    ;12CO spectra
    bspecU0=reform(datUa[coord2pix(hdrUa,bpeak1[j],1),  coord2pix(hdrUa,bpeak2[j],2),  *])
    bspecU1=reform(datUa[coord2pix(hdrUa,bpeak1[j],1)-1,coord2pix(hdrUa,bpeak2[j],2),  *])
    bspecU2=reform(datUa[coord2pix(hdrUa,bpeak1[j],1)+1,coord2pix(hdrUa,bpeak2[j],2),  *])
    bspecU3=reform(datUa[coord2pix(hdrUa,bpeak1[j],1),  coord2pix(hdrUa,bpeak2[j],2)-1,*])
    bspecU4=reform(datUa[coord2pix(hdrUa,bpeak1[j],1),  coord2pix(hdrUa,bpeak2[j],2)+1,*])
    bspec_U=(bspecU0+bspecU1+bspecU2+bspecU3+bspecU4)/5.
    rspecU0=reform(datUa[coord2pix(hdrUa,rpeak1[k],1),  coord2pix(hdrUa,rpeak2[k],2),  *])
    rspecU1=reform(datUa[coord2pix(hdrUa,rpeak1[k],1)-1,coord2pix(hdrUa,rpeak2[k],2),  *])
    rspecU2=reform(datUa[coord2pix(hdrUa,rpeak1[k],1)+1,coord2pix(hdrUa,rpeak2[k],2),  *])
    rspecU3=reform(datUa[coord2pix(hdrUa,rpeak1[k],1),  coord2pix(hdrUa,rpeak2[k],2)-1,*])
    rspecU4=reform(datUa[coord2pix(hdrUa,rpeak1[k],1),  coord2pix(hdrUa,rpeak2[k],2)+1,*])
    rspec_U=(rspecU0+rspecU1+rspecU2+rspecU3+rspecU4)/5.
    ;13CO spectrum
    specL0=reform(datLa[coord2pix(hdrLa,bpeak1[j],1),  coord2pix(hdrLa,bpeak2[j],2),  *])
    specL1=reform(datLa[coord2pix(hdrLa,bpeak1[j],1)-1,coord2pix(hdrLa,bpeak2[j],2),  *])
    specL2=reform(datLa[coord2pix(hdrLa,bpeak1[j],1)+1,coord2pix(hdrLa,bpeak2[j],2),  *])
    specL3=reform(datLa[coord2pix(hdrLa,bpeak1[j],1),  coord2pix(hdrLa,bpeak2[j],2)-1,*])
    specL4=reform(datLa[coord2pix(hdrLa,bpeak1[j],1),  coord2pix(hdrLa,bpeak2[j],2)+1,*])
    specL5=reform(datLa[coord2pix(hdrLa,rpeak1[k],1),  coord2pix(hdrLa,rpeak2[k],2),  *])
    specL6=reform(datLa[coord2pix(hdrLa,rpeak1[k],1)-1,coord2pix(hdrLa,rpeak2[k],2),  *])
    specL7=reform(datLa[coord2pix(hdrLa,rpeak1[k],1)+1,coord2pix(hdrLa,rpeak2[k],2),  *])
    specL8=reform(datLa[coord2pix(hdrLa,rpeak1[k],1),  coord2pix(hdrLa,rpeak2[k],2)-1,*])
    specL9=reform(datLa[coord2pix(hdrLa,rpeak1[k],1),  coord2pix(hdrLa,rpeak2[k],2)+1,*])
    spec_L=(specL0+specL1+specL2+specL3+specL4+specL5+specL6+specL7+specL8+specL9)/10.
    ;vLpeak
    vLpeak0=mean(peakvmap[(coord2pix(hdrLa,bpeak1[j],1)-1):(coord2pix(hdrLa,bpeak1[j],1)+1),$
      (coord2pix(hdrLa,bpeak2[j],2)-1):(coord2pix(hdrLa,bpeak2[j],2)+1)],/nan)
    if finite(vLpeak0,/nan) eq 1 then vLpeak0=Lvrange[0]
    vLpeak1=mean(peakvmap[(coord2pix(hdrLa,rpeak1[k],1)-1):(coord2pix(hdrLa,rpeak1[k],1)+1),$
      (coord2pix(hdrLa,rpeak2[k],2)-1):(coord2pix(hdrLa,rpeak2[k],2)+1)],/nan)
    if finite(vLpeak1,/nan) eq 1 then vLpeak1=Lvrange[1]
    vLpeak=ceil((vLpeak0+vLpeak1)*5)*0.1
;    if i eq 4 then begin vrange=[-58,-22]
;    endif else begin
    vrange=[(Lvrange[1]-Ubvmap[coord2pix(bvhdr,bpeak1[j],1),coord2pix(bvhdr,bpeak2[j],2)]-3)<(vLpeak-5),$
      (Urvmap[coord2pix(rvhdr,rpeak1[k],1),coord2pix(rvhdr,rpeak2[k],2)]+Lvrange[0]+3)>(vLpeak+5)]
;    endelse  
     vrange=[vvv1[i],vvv2[i]]
    
    ;plot spectra
    psname='out_'+num2str(num[i])+'_'+num2str(index[i])+'.eps'
    cgps_open,psname,font=!p.font,/quiet,default_thickness=1.0,charsize=0.9,/portrait
    cgDisplay, xsize=round(xsize1), ysize=round(ysize1)
    cgplot,vU,bspec_U,psym=10,/nodata,xrange=vrange,color='blue',yrange=[0-max([bspec_U,rspec_U])*0.15,max([bspec_U,rspec_U])*1.15],$
      ytitle='T!DMB!N (K)',position=p3,yminor=5,xtitle='LSR Velocity (km s!U-1!N)'
      ;title='OUTFLOW BLUELOBE CANDIDATE '+num2str(bsn[i])+': '+num2str(bpeak1[i],format='(f7.3)')+num2str(bpeak2[i],format='(f+6.3)')
    bpolyrg=where(vU ge bluev0[i] and vU le bluev1[i])
    bx1_arr=array_band(vU[bpolyrg],/band)
    bx2_arr=array_band(Reverse(vU[bpolyrg]),/band)
    by1_arr=array_band(bspec_U[bpolyrg])
    by2_arr=array_band(baseline[bpolyrg])
    cgColorFill,[bx1_arr, bx2_arr, bx1_arr[0]],[by1_arr,by2_arr,by1_arr[0]],Color='Sky Blue'
    rpolyrg=where(vU ge redv0[i] and vU le redv1[i])
    rx1_arr=array_band(vU[rpolyrg],/band)
    rx2_arr=array_band(Reverse(vU[rpolyrg]),/band)
    ry1_arr=array_band(rspec_U[rpolyrg])
    ry2_arr=array_band(baseline[rpolyrg])
    cgColorFill,[rx1_arr, rx2_arr, rx1_arr[0]],[ry1_arr,ry2_arr,ry1_arr[0]],Color='pink'
    cgplot,vU,bspec_U,psym=10,color='blue',/overplot
    cgplot,vU,rspec_U,psym=10,color='red', /overplot
    cgplot,vL,spec_L,psym=10,color='green',/overplot
    cgplot, !x.CRange,[0.,0.], LineStyle=0, Color='black',/overplot
    cgplot,[vLpeak,vLpeak],[0,!y.crange[1]],linestyle=1,/over
           
    cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(d)',color='black'
    cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]-0.15*(!y.CRange[1]-!y.CRange[0]),$
      rank[i],color='black'  

    
    ;plot CO outflow contours ?and WISE?
    cropfits,/dataform,dathdrL,[max([bpeak1[j],rpeak1[k]])+outsz/60.,min([bpeak1[j],rpeak1[k]])-outsz/60.],$
      [min([bpeak2[j],rpeak2[k]])-outsz/60.,max([bpeak2[j],rpeak2[k]])+outsz/60.],[bw1[j],rw1[k]],output=lcr
    ldata=smooth(lcr.dat,[1,1,3],/edge_mirror)
    liimap=total(ldata,3)*abs(sxpar(lcr.hdr,'CDELT3'))/1000.
    liimap=congrid(liimap,3*n_elements(liimap[*,0]),3*n_elements(liimap[0,*]),cubic=-0.5,/MINUS_ONE)
    liimap=smooth(liimap,[3,3],/edge_mirror)
    ;lmax=max(liimap[round(n_elements(liimap[*,0])/2.-n_elements(liimap[*,0])/12.-1):round(n_elements(liimap[*,0])/2.+n_elements(liimap[*,0])/12.-1),$
    ;  round(n_elements(liimap[0,*])/2.-n_elements(liimap[0,*])/12.-1):round(n_elements(liimap[0,*])/2.+n_elements(liimap[0,*])/12.-1)])
    ;llevels=lmax*(0.1-0.01+indgen(15)*(1-0.9)*0.2)
    ;print,lmax
        
    cropfits,/dataform,dathdr,[max([bpeak1[j],rpeak1[k]])+outsz/60.,min([bpeak1[j],rpeak1[k]])-outsz/60.],$
      [min([bpeak2[j],rpeak2[k]])-outsz/60.,max([bpeak2[j],rpeak2[k]])+outsz/60.],[bw0[j],bw1[j]],output=bcr
    cropfits,/dataform,dathdr,[max([bpeak1[j],rpeak1[k]])+outsz/60.,min([bpeak1[j],rpeak1[k]])-outsz/60.],$
      [min([bpeak2[j],rpeak2[k]])-outsz/60.,max([bpeak2[j],rpeak2[k]])+outsz/60.],[rw0[k],rw1[k]],output=rcr
    bdata=smooth(bcr.dat,[1,1,3],/edge_mirror)
    rdata=smooth(rcr.dat,[1,1,3],/edge_mirror)
    biimap=total(bdata,3)*abs(sxpar(bcr.hdr,'CDELT3'))/1000.
    biimap=congrid(biimap,3*n_elements(biimap[*,0]),3*n_elements(biimap[0,*]),cubic=-0.5,/MINUS_ONE)
    biimap=smooth(biimap,[3,3],/edge_mirror)
    blevels=bmax[i]*((blow[i]<0.8)-0.01+indgen(15)*(1-(blow[i]<0.8))*0.2)
    riimap=total(rdata,3)*abs(sxpar(rcr.hdr,'CDELT3'))/1000.
    riimap=congrid(riimap,3*n_elements(riimap[*,0]),3*n_elements(riimap[0,*]),cubic=-0.5,/MINUS_ONE)
    riimap=smooth(riimap,[3,3],/edge_mirror)
    rlevels=rmax[i]*((rlow[i]<0.8)-0.01+indgen(15)*(1-(rlow[i]<0.8))*0.2)  
    crp1=sxpar(bcr.hdr,'CRPIX1')
    crv1=sxpar(bcr.hdr,'CRVAL1')
    del1=sxpar(bcr.hdr,'CDELT1')
    axs1=sxpar(bcr.hdr,'NAXIS1')
    l_l=(360.+(0.5-crp1)*del1+crv1) mod 360
    l_r=(360.+(axs1+0.5-crp1)*del1+crv1) mod 360
;    l_l=l_l+0.02
;    l_r=l_r-0.02
    crp2=sxpar(bcr.hdr,'CRPIX2')
    crv2=sxpar(bcr.hdr,'CRVAL2')
    del2=sxpar(bcr.hdr,'CDELT2')
    axs2=sxpar(bcr.hdr,'NAXIS2')
    b_d=(0.5-crp2)*del2+crv2
    b_u=(axs2+0.5-crp2)*del2+crv2
    x_range=[l_l,l_r]
    y_range=[b_d,b_u]
    ratio=abs((b_u-b_d)/(l_l-l_r))
    ;pos1=[100./xsize1,75./ysize1,800./xsize1,587.5/ysize1]
    pos1=p1
    ;print,pos1
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase,xtickformat='(a1)',ytickformat='(a1)'
    cgloadct,53,clip=[0,255-64],ncolors=15
    ;cgimage,bytscl(liimap),position=pos1,/noerase
    cgcontour,liimap,/onimage,nlevels=15,label=0,c_colors=indgen(15),/fill
    
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')
    cgloadct,3,ncolors=15,/reverse,clip=[64,255-64]
    cgcontour,riimap,/onimage,levels=rlevels,label=0,c_colors=indgen(15),thick=3
    cgloadct,1,ncolors=15,/reverse,clip=[64,255-64]
    cgcontour,biimap,/onimage,levels=blevels,label=0,c_colors=indgen(15),thick=3
    ;arrow
    if bpeak1[j] eq rpeak1[k] then k_br=999999999.0 else k_br=(rpeak2[k]-bpeak2[j])/(rpeak1[k]-bpeak1[j])
    k_theta=atan(k_br)
    x_0=bpeak1[j]-(outsz/60)*abs(cos(k_theta))*(-1)^((rpeak1[k]-bpeak1[j]) lt 0)
    y_0=bpeak2[j]-(outsz/60)*abs(sin(k_theta))*(-1)^((rpeak2[k]-bpeak2[j]) lt 0)
    x_1=rpeak1[k]+(outsz/60)*abs(cos(k_theta))*(-1)^((rpeak1[k]-bpeak1[j]) lt 0)
    y_1=rpeak2[k]+(outsz/60)*abs(sin(k_theta))*(-1)^((rpeak2[k]-bpeak2[j]) lt 0)
    cgplot,/over,glon[i],glat[i],psym=cgsymcat(16),color='purple'
    cgarrow,x_0,y_0,x_1,y_1,/data,color='black',/clip,hsize=!D.X_SIZE/128.,hthick=2.0
;    if STRMID(region,0,1) eq 'G' then IDstr= STRMID(region,0,1)+STRMID(region,0,1,/reverse)+'-'+num2str(num[i]) $
;      else IDstr= STRUPCASE(STRMID(region,7,5))+'-'+num2str(num[i])
       IDstr=num2str(index[i])
    cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),$
      IDstr,color='black'
      
    points=(2*!PI/99.0)*findgen(100)
    beam_x=l_l-0.01+24.5/3600.0*cos(points)
    beam_y=b_d+0.01+24.5/3600.0*sin(points)
 
    cgplot,beam_x,beam_y, psym=0,SYMSIZE=1,/data,color='black',xrange=x_range,yrange=y_range,Aspect=ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')
      
    cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(c)',color='black'
    beam_xx=min(beam_x)
    scalex_c=[beam_xx-0.01,beam_xx-0.03]
    scaley_c=[b_d+0.01,b_d+0.01]
    distance_phy=round((0.02*3600.d*distance1[i]*1000.d/206265)*100)/100000.d
    dist_phy=num2str(distance_phy)
    dist_phy=strmid(dist_phy,0,4) 
    cgplot,scalex_c,scaley_c, psym=0,SYMSIZE=1,/data,color='black',xrange=x_range,yrange=y_range,Aspect=ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})') 
    cgtext, beam_xx-0.01 ,b_d+0.015,dist_phy+' pc',color='black',/data  
    
    wiseBname='wise_4.6_'+num2str(bnum[i])+'-'+num2str(rnum[i])+'.fits'
    wiseGname='wise_12_'+num2str(bnum[i])+'-'+num2str(rnum[i])+'.fits'
    wiseRname='wise_22_'+num2str(bnum[i])+'-'+num2str(rnum[i])+'.fits'
    fits_read,'./wisedata/'+wiseBname,datwB,hdrwB
    fits_read,'./wisedata/'+wiseGname,datwG,hdrwG
    fits_read,'./wisedata/'+wiseRname,datwR,hdrwR
    dathdrwB=list(datwB,hdrwB)
    dathdrwG=list(datwG,hdrwG)
    dathdrwR=list(datwR,hdrwR)
    whdr=headfits('./wisedata/'+wiseBname)
    cl_l=l_l+sxpar(whdr,'CDELT1')
    cl_r=l_r-sxpar(whdr,'CDELT1')
    cb_d=b_d+sxpar(whdr,'CDELT2')
    cb_u=b_u-sxpar(whdr,'CDELT2')
    cropfits,dathdrwB,[cl_r,cl_l],[cb_u,cb_d],output=wB,/dataform
    cropfits,dathdrwG,[cl_r,cl_l],[cb_u,cb_d],output=wG,/dataform
    cropfits,dathdrwR,[cl_r,cl_l],[cb_u,cb_d],output=wR,/dataform
    wiseimg_B=wise_scl(wB.dat)
    wiseimg_G=wise_scl(wG.dat)
    wiseimg_R=wise_scl(wR.dat)  
    wiseimg_RGB=[wiseimg_R,wiseimg_G,wiseimg_B]
    pos0=p0
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos0,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')
    cgimage,wiseimg_RGB,position=pos0,/noerase
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='white',position=pos0,/noerase,xtickformat='(a1)',ytickformat='(a1)'
    cgloadct,3,ncolors=15,/reverse,clip=[64,255-64]
    cgcontour,riimap,/onimage,levels=rlevels,label=0,c_colors=indgen(15)
    cgloadct,1,ncolors=15,/reverse,clip=[64,255-64]
    cgcontour,biimap,/onimage,levels=blevels,label=0,c_colors=indgen(15)
    cgplot,/over,glon[i],glat[i],psym=cgsymcat(16),color='purple'
    
    cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(a)',color='white'
    cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),$
      textoidl('WISE 4.6 12 22 \mum'),color='white'  
    scalex_w=[l_l-0.01,l_l-0.03]
    scaley_w=[b_d+0.01,b_d+0.01]
    distance_phy=round((0.02*3600.d*distance1[i]*1000.d/206265)*100)/100000.d
    dist_phy=num2str(distance_phy)
    dist_phy=strmid(dist_phy,0,4)
    cgplot,scalex_w,scaley_w, psym=0,SYMSIZE=1,/data,color='white',xrange=x_range,yrange=y_range,Aspect=ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='white',position=pos0,/noerase,xtickformat='(a1)',ytickformat='(a1)';ytitle='(a1)',xtitle=textoidl('Galactic Longitude (^{o})')
    cgtext, l_l-0.01 ,b_d+0.015,dist_phy+' pc',color='white',/data

    ;pv diagram
    cropfits,dathdr,vrange,dim='v',/dataform,output=pvout
    pv=mkpvbelt(pvout.dat,pvout.hdr,[x_0,x_1],[y_0,y_1],2,/gal)
    pvdata=pv.dat
    pvdata=smooth(pvdata,[3,3],/edge_mirror)
    pvhdr=pv.hdr
    pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
    pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
    pvvrg=[min(pvv),max(pvv)]
    pvprg=[min(pvp),max(pvp)]*60.0-mean([min(pvp),max(pvp)])*60;*sqrt(2)   
    cgplot,[0],[0],/nodata,/noerase,position=p2,xrange=pvvrg,yrange=pvprg,$
      ytickformat='(i)',xtickformat='(a1)',yminor=5
    cgloadct,50,clip=[64,240]
    pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
      round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
      *(0.05+indgen(10)*0.1)
    cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline;',outcolor='green'
    cgplot,[0],[0],/nodata,/noerase,position=p2,xrange=pvvrg,yrange=pvprg,$
      ytickinterval=2.0,ytitle="Position (')",yminor=5,xtitle='LSR Velocity (km s!U-1!N)'
    cgPlot, !x.crange,[0,0], LineStyle=1, /overplot
    cgPlot, [bw0[j],bw0[j]],!y.CRange, LineStyle=1, /overplot
    cgPlot, [bw1[j],bw1[j]],!y.CRange, LineStyle=1, /overplot
    cgPlot, [rw0[k],rw0[k]],!y.CRange, LineStyle=1, /overplot
    cgPlot, [rw1[k],rw1[k]],!y.CRange, LineStyle=1, /overplot
    cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
    
    pvpx=mean([x_0,x_1])
    pvpxs=strmid(num2str(pvpx),0,7)
    pvpy=mean([y_0,y_1])
    pa=round(atan((y_1-y_0)/(x_0-x_1))*180/!pi)
    if (pa lt 0) then begin 
      if (y_1-y_0 gt 0) then pa=pa+180
    endif else begin 
      if (x_0-x_1 lt 0) then pa=pa+180
    endelse
    pam=round(pa)
    if strmid(num2str(pvpy),0,1) eq '-' then pvpys=strmid(num2str(pvpy),0,6) else pvpys=strmid(num2str(pvpy),0,5)
    pvstring='offset from '+pvpxs+' '+pvpys+' at PA '+num2str(pam)+textoidl('^{o}')
    cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]+0.05*(!y.CRange[1]-!y.CRange[0]), pvstring, color='black'
    cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(b)',color='black' 
    
    cgps_close
  endif
  
  if tag[i] eq 'B' then begin
    j=where(bsn eq bnum[i])
    ;12CO spectra
    bspecU0=reform(datUa[coord2pix(hdrUa,bpeak1[j],1),  coord2pix(hdrUa,bpeak2[j],2),  *])
    bspecU1=reform(datUa[coord2pix(hdrUa,bpeak1[j],1)-1,coord2pix(hdrUa,bpeak2[j],2),  *])
    bspecU2=reform(datUa[coord2pix(hdrUa,bpeak1[j],1)+1,coord2pix(hdrUa,bpeak2[j],2),  *])
    bspecU3=reform(datUa[coord2pix(hdrUa,bpeak1[j],1),  coord2pix(hdrUa,bpeak2[j],2)-1,*])
    bspecU4=reform(datUa[coord2pix(hdrUa,bpeak1[j],1),  coord2pix(hdrUa,bpeak2[j],2)+1,*])
    bspec_U=(bspecU0+bspecU1+bspecU2+bspecU3+bspecU4)/5.
    ;13CO spectrum
    bspecL0=reform(datLa[coord2pix(hdrLa,bpeak1[j],1),  coord2pix(hdrLa,bpeak2[j],2),  *])
    bspecL1=reform(datLa[coord2pix(hdrLa,bpeak1[j],1)-1,coord2pix(hdrLa,bpeak2[j],2),  *])
    bspecL2=reform(datLa[coord2pix(hdrLa,bpeak1[j],1)+1,coord2pix(hdrLa,bpeak2[j],2),  *])
    bspecL3=reform(datLa[coord2pix(hdrLa,bpeak1[j],1),  coord2pix(hdrLa,bpeak2[j],2)-1,*])
    bspecL4=reform(datLa[coord2pix(hdrLa,bpeak1[j],1),  coord2pix(hdrLa,bpeak2[j],2)+1,*])
    bspec_L=(bspecL0+bspecL1+bspecL2+bspecL3+bspecL4)/5.
    ;vLpeak
    vLpeak0=mean(peakvmap[(coord2pix(hdrLa,bpeak1[j],1)-1):(coord2pix(hdrLa,bpeak1[j],1)+1),$
      (coord2pix(hdrLa,bpeak2[j],2)-1):(coord2pix(hdrLa,bpeak2[j],2)+1)],/nan)
    if finite(vLpeak0,/nan) eq 1 then vLpeak0=Lvrange[0]
    vLpeak=ceil(vLpeak0*5)*0.2
    vrange=[(Lvrange[1]-Ubvmap[coord2pix(bvhdr,bpeak1[j],1),coord2pix(bvhdr,bpeak2[j],2)]-3)<(vLpeak-5),$
      (Urvmap[coord2pix(rvhdr,bpeak1[j],1),coord2pix(rvhdr,bpeak2[j],2)]+Lvrange[0]+3)>(vLpeak+5)]
     vrange=[vvv1[i],vvv2[i]]
    
    ;plot spectra
    psname='out_'+num2str(num[i])+'_'+num2str(index[i])+'.eps'
    cgps_open,psname,font=!p.font,/quiet,default_thickness=1.0,charsize=0.9,/portrait
    cgDisplay, xsize=round(xsize1), ysize=round(ysize1)
    cgplot,vU,bspec_U,psym=10,/nodata,xrange=vrange,color='blue',yrange=[0-max(bspec_U)*0.15,max(bspec_U)*1.15],$
      ytitle='T!DMB!N (K)',position=p3,yminor=5,xtitle='LSR Velocity (km s!U-1!N)'
      ;title='OUTFLOW BLUELOBE CANDIDATE '+num2str(bsn[i])+': '+num2str(bpeak1[i],format='(f7.3)')+num2str(bpeak2[i],format='(f+6.3)')
    bpolyrg=where(vU ge bluev0[i] and vU le bluev1[i])
    bx1_arr=array_band(vU[bpolyrg],/band)
    bx2_arr=array_band(Reverse(vU[bpolyrg]),/band)
    by1_arr=array_band(bspec_U[bpolyrg])
    by2_arr=array_band(baseline[bpolyrg])
    cgColorFill,[bx1_arr, bx2_arr, bx1_arr[0]],[by1_arr,by2_arr,by1_arr[0]],Color='Sky Blue'
    cgplot,vU,bspec_U,psym=10,color='blue',/overplot
    cgplot,vL,bspec_L,psym=10,color='green',/overplot
    cgplot, !x.CRange,[0.,0.], LineStyle=0, Color='black',/overplot
    cgplot,[vLpeak,vLpeak],[0,!y.crange[1]],linestyle=1,/over
    
    cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(d)',color='black'
    cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]-0.15*(!y.CRange[1]-!y.CRange[0]),$
      rank[i],color='black'  
    
    ;plot CO outflow contours ?and WISE?
    cropfits,/dataform,dathdrL,[bpeak1[j]+outsz/60.,bpeak1[j]-outsz/60.],$
      [bpeak2[j]-outsz/60.,bpeak2[j]+outsz/60.],[bw1[j],2*vLpeak-bw1[j]],output=lcr
    ldata=smooth(lcr.dat,[1,1,3],/edge_mirror)
    liimap=total(ldata,3)*abs(sxpar(lcr.hdr,'CDELT3'))/1000.
    liimap=congrid(liimap,3*n_elements(liimap[*,0]),3*n_elements(liimap[0,*]),cubic=-0.5,/MINUS_ONE)
    liimap=smooth(liimap,[3,3],/edge_mirror)
    ;lmax=max(liimap[round(n_elements(liimap[*,0])/2.-n_elements(liimap[*,0])/12.-1):round(n_elements(liimap[*,0])/2.+n_elements(liimap[*,0])/12.-1),$
    ;  round(n_elements(liimap[0,*])/2.-n_elements(liimap[0,*])/12.-1):round(n_elements(liimap[0,*])/2.+n_elements(liimap[0,*])/12.-1)])
    ;llevels=lmax[i]*(0.1-0.01+indgen(15)*(1-0.9))*0.2)
    ;print,lmax
        
    cropfits,/dataform,dathdr,[bpeak1[j]+outsz/60.,bpeak1[j]-outsz/60.],$
      [bpeak2[j]-outsz/60.,bpeak2[j]+outsz/60.],[bw0[j],bw1[j]],output=bcr
    bdata=smooth(bcr.dat,[1,1,3],/edge_mirror)
    biimap=total(bdata,3)*abs(sxpar(bcr.hdr,'CDELT3'))/1000.
    biimap=congrid(biimap,3*n_elements(biimap[*,0]),3*n_elements(biimap[0,*]),cubic=-0.5,/MINUS_ONE)
    biimap=smooth(biimap,[3,3],/edge_mirror)
    blevels=bmax[i]*((blow[i]<0.8)-0.01+indgen(15)*(1-(blow[i]<0.8))*0.2) 
    crp1=sxpar(bcr.hdr,'CRPIX1')
    crv1=sxpar(bcr.hdr,'CRVAL1')
    del1=sxpar(bcr.hdr,'CDELT1')
    axs1=sxpar(bcr.hdr,'NAXIS1')
    l_l=(360.+(0.5-crp1)*del1+crv1) mod 360
    l_r=(360.+(axs1+0.5-crp1)*del1+crv1) mod 360
    crp2=sxpar(bcr.hdr,'CRPIX2')
    crv2=sxpar(bcr.hdr,'CRVAL2')
    del2=sxpar(bcr.hdr,'CDELT2')
    axs2=sxpar(bcr.hdr,'NAXIS2')
    b_d=(0.5-crp2)*del2+crv2
    b_u=(axs2+0.5-crp2)*del2+crv2
    x_range=[l_l,l_r]
    y_range=[b_d,b_u]
    ratio=abs((b_u-b_d)/(l_l-l_r))
    ;pos1=[100./xsize1,75./ysize1,800./xsize1,587.5/ysize1]
    pos1=p1
    ;print,pos1
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase
    cgloadct,53,clip=[0,255-64],ncolors=15
    ;cgimage,bytscl(liimap),position=pos1,/noerase
    cgcontour,liimap,/onimage,nlevels=15,label=0,c_colors=indgen(15),/fill
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')
    cgloadct,1,ncolors=15,/reverse,clip=[64,255-64]
    cgcontour,biimap,/onimage,levels=blevels,label=0,c_colors=indgen(15),thick=3
    ;arrow
    x_0=bpeak1[j]+(outsz/45)*sqrt(0.5)
    y_0=bpeak2[j]-(outsz/45)*sqrt(0.5)
    x_1=bpeak1[j]-(outsz/45)*sqrt(0.5)
    y_1=bpeak2[j]+(outsz/45)*sqrt(0.5)
    cgplot,/over,glon[i],glat[i],psym=cgsymcat(16),color='purple'
    cgarrow,x_0,y_0,x_1,y_1,/data,color='black',/clip,hsize=!D.X_SIZE/128.,hthick=2.0
;    if STRMID(region,0,1) eq 'G' then IDstr= STRMID(region,0,1)+STRMID(region,0,1,/reverse)+'-'+num2str(num[i]) $
;      else IDstr= STRUPCASE(STRMID(region,7,5))+'-'+num2str(num[i])
    IDstr=num2str(index[i])
    cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),$
      IDstr,color='black'
      
    cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(c)',color='black'    
      
    points=(2*!PI/99.0)*findgen(100)
    beam_x=l_l-0.01+24.5/3600.0*cos(points)
    beam_y=b_d+0.01+24.5/3600.0*sin(points)
 
    cgplot,beam_x,beam_y, psym=0,SYMSIZE=1,/data,color='black',xrange=x_range,yrange=y_range,Aspect=ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')    
      
    beam_xx=min(beam_x)
    scalex_c=[beam_xx-0.01,beam_xx-0.03]
    scaley_c=[b_d+0.01,b_d+0.01]
    distance_phy=round((0.02*3600.d*distance1[i]*1000.d/206265)*100)/100000.d
    dist_phy=num2str(distance_phy)
    dist_phy=strmid(dist_phy,0,4) 
    cgplot,scalex_c,scaley_c, psym=0,SYMSIZE=1,/data,color='black',xrange=x_range,yrange=y_range,Aspect=ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})') 
    cgtext, beam_xx-0.01 ,b_d+0.015,dist_phy+' pc',color='black',/data  
          
    
    wiseBname='wise_4.6_'+num2str(bnum[i])+'-'+num2str(rnum[i])+'.fits'
    wiseGname='wise_12_'+num2str(bnum[i])+'-'+num2str(rnum[i])+'.fits'
    wiseRname='wise_22_'+num2str(bnum[i])+'-'+num2str(rnum[i])+'.fits'
    fits_read,'./wisedata/'+wiseBname,datwB,hdrwB
    fits_read,'./wisedata/'+wiseGname,datwG,hdrwG
    fits_read,'./wisedata/'+wiseRname,datwR,hdrwR
    dathdrwB=list(datwB,hdrwB)
    dathdrwG=list(datwG,hdrwG)
    dathdrwR=list(datwR,hdrwR)
    whdr=headfits('./wisedata/'+wiseBname)
    cl_l=l_l+sxpar(whdr,'CDELT1')
    cl_r=l_r-sxpar(whdr,'CDELT1')
    cb_d=b_d+sxpar(whdr,'CDELT2')
    cb_u=b_u-sxpar(whdr,'CDELT2')
    cropfits,dathdrwB,[cl_r,cl_l],[cb_u,cb_d],output=wB,/dataform
    cropfits,dathdrwG,[cl_r,cl_l],[cb_u,cb_d],output=wG,/dataform
    cropfits,dathdrwR,[cl_r,cl_l],[cb_u,cb_d],output=wR,/dataform
    wiseimg_B=wise_scl(wB.dat)
    wiseimg_G=wise_scl(wG.dat)
    wiseimg_R=wise_scl(wR.dat)  
    wiseimg_RGB=[wiseimg_R,wiseimg_G,wiseimg_B]
    pos0=p0
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos0,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')
    cgimage,wiseimg_RGB,position=pos0,/noerase
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='white',position=pos0,/noerase,xtickformat='(a1)',ytickformat='(a1)'
    cgloadct,1,ncolors=15,/reverse,clip=[64,255-64]
    cgcontour,biimap,/onimage,levels=blevels,label=0,c_colors=indgen(15)
    cgplot,/over,glon[i],glat[i],psym=cgsymcat(16),color='purple'
    
     cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(a)',color='white'
     cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),$
      textoidl('WISE 4.6 12 22 \mum'),color='white'  
    scalex_w=[l_l-0.01,l_l-0.03]
    scaley_w=[b_d+0.01,b_d+0.01]
    distance_phy=round((0.02*3600.d*distance1[i]*1000.d/206265)*100)/100000.d
    dist_phy=num2str(distance_phy)
    dist_phy=strmid(dist_phy,0,4)
    cgplot,scalex_w,scaley_w, psym=0,SYMSIZE=1,/data,color='white',xrange=x_range,yrange=y_range,Aspect=ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='white',position=pos0,/noerase,xtickformat='(a1)',ytickformat='(a1)';ytitle='(a1)',xtitle=textoidl('Galactic Longitude (^{o})')
    cgtext, l_l-0.01 ,b_d+0.015,dist_phy+' pc',color='white',/data     

    ;pv diagram
    cropfits,dathdr,vrange,dim='v',/dataform,output=pvout
    pv=mkpvbelt(pvout.dat,pvout.hdr,[x_0,x_1],[y_0,y_1],2,/gal)
    pvdata=pv.dat
    pvdata=smooth(pvdata,[3,3],/edge_mirror)
    pvhdr=pv.hdr
    pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
    pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
    pvvrg=[min(pvv),max(pvv)]
    pvprg=[min(pvp),max(pvp)]*60.0-mean([min(pvp),max(pvp)])*60;*sqrt(2)   
    cgplot,[0],[0],/nodata,/noerase,position=p2,xrange=pvvrg,yrange=pvprg,$
      ytickformat='(i)',xtickformat='(a1)',yminor=5
    cgloadct,50,clip=[64,240]
    pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
      round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
      *(0.05+indgen(10)*0.1)
    cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline;',outcolor='green'
    cgplot,[0],[0],/nodata,/noerase,position=p2,xrange=pvvrg,yrange=pvprg,$
      ytickinterval=2.0,ytitle="Position (')",yminor=5,xtitle='LSR Velocity (km s!U-1!N)'
    cgPlot, !x.crange,[0,0], LineStyle=1, /overplot
    cgPlot, [bw0[j],bw0[j]],!y.CRange, LineStyle=1, /overplot
    cgPlot, [bw1[j],bw1[j]],!y.CRange, LineStyle=1, /overplot
    cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
    
        pvpx=mean([x_0,x_1])
    pvpxs=strmid(num2str(pvpx),0,7)
    pvpy=mean([y_0,y_1])
    if strmid(num2str(pvpy),0,1) eq '-' then pvpys=strmid(num2str(pvpy),0,6) else pvpys=strmid(num2str(pvpy),0,5)
    pvstring='offset from '+pvpxs+' '+pvpys+' at PA 45'+textoidl('^{o}')
    cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]+0.05*(!y.CRange[1]-!y.CRange[0]), pvstring, color='black'
    cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(b)',color='black'        
    cgps_close
  endif
  
    if tag[i] eq 'R' then begin
    k=where(rsn eq rnum[i])
    ;12CO spectra
    rspecU0=reform(datUa[coord2pix(hdrUa,rpeak1[k],1),  coord2pix(hdrUa,rpeak2[k],2),  *])
    rspecU1=reform(datUa[coord2pix(hdrUa,rpeak1[k],1)-1,coord2pix(hdrUa,rpeak2[k],2),  *])
    rspecU2=reform(datUa[coord2pix(hdrUa,rpeak1[k],1)+1,coord2pix(hdrUa,rpeak2[k],2),  *])
    rspecU3=reform(datUa[coord2pix(hdrUa,rpeak1[k],1),  coord2pix(hdrUa,rpeak2[k],2)-1,*])
    rspecU4=reform(datUa[coord2pix(hdrUa,rpeak1[k],1),  coord2pix(hdrUa,rpeak2[k],2)+1,*])
    rspec_U=(rspecU0+rspecU1+rspecU2+rspecU3+rspecU4)/5.
    ;13CO spectrum
    rspecL5=reform(datLa[coord2pix(hdrLa,rpeak1[k],1),  coord2pix(hdrLa,rpeak2[k],2),  *])
    rspecL6=reform(datLa[coord2pix(hdrLa,rpeak1[k],1)-1,coord2pix(hdrLa,rpeak2[k],2),  *])
    rspecL7=reform(datLa[coord2pix(hdrLa,rpeak1[k],1)+1,coord2pix(hdrLa,rpeak2[k],2),  *])
    rspecL8=reform(datLa[coord2pix(hdrLa,rpeak1[k],1),  coord2pix(hdrLa,rpeak2[k],2)-1,*])
    rspecL9=reform(datLa[coord2pix(hdrLa,rpeak1[k],1),  coord2pix(hdrLa,rpeak2[k],2)+1,*])
    rspec_L=(rspecL5+rspecL6+rspecL7+rspecL8+rspecL9)/5.
    ;vLpeak
;    if k eq 1 then begin
;       vLpeak=-2.5
;    endif else begin 
    vLpeak1=mean(peakvmap[(coord2pix(hdrLa,rpeak1[k],1)-1):(coord2pix(hdrLa,rpeak1[k],1)+1),$
      (coord2pix(hdrLa,rpeak2[k],2)-1):(coord2pix(hdrLa,rpeak2[k],2)+1)],/nan)
;    endelse 
    if finite(vLpeak1,/nan) eq 1 then vLpeak1=Lvrange[1]  
    vLpeak=ceil(vLpeak1*5)*0.2
;    if num[i] eq 36 then vLpeak=4.2
;        if i eq 35 then begin vrange=[-58,-42]
;    endif else begin
    vrange=[(Lvrange[1]-Ubvmap[coord2pix(bvhdr,rpeak1[k],1),coord2pix(bvhdr,rpeak2[k],2)]-3)<(vLpeak-5),$
      (Urvmap[coord2pix(rvhdr,rpeak1[k],1),coord2pix(rvhdr,rpeak2[k],2)]+Lvrange[0]+3)>(vLpeak+5)]
;    endelse
     vrange=[vvv1[i],vvv2[i]]
    
    ;plot spectra
    psname='out_'+num2str(num[i])+'_'+num2str(index[i])+'.eps'
    cgps_open,psname,font=!p.font,/quiet,default_thickness=1.0,charsize=0.9,/portrait
    cgDisplay, xsize=round(xsize1), ysize=round(ysize1)
    cgplot,vU,rspec_U,psym=10,/nodata,xrange=vrange,color='blue',yrange=[0-max(rspec_U)*0.15,max(rspec_U)*1.15],$
      ytitle='T!DMB!N (K)',position=p3,yminor=5,xtitle='LSR Velocity (km s!U-1!N)'
      ;title='OUTFLOW BLUELOBE CANDIDATE '+num2str(bsn[i])+': '+num2str(bpeak1[i],format='(f7.3)')+num2str(bpeak2[i],format='(f+6.3)')
    rpolyrg=where(vU ge redv0[i] and vU le redv1[i])
    rx1_arr=array_band(vU[rpolyrg],/band)
    rx2_arr=array_band(Reverse(vU[rpolyrg]),/band)
    ry1_arr=array_band(rspec_U[rpolyrg])
    ry2_arr=array_band(baseline[rpolyrg])
    cgColorFill,[rx1_arr, rx2_arr, rx1_arr[0]],[ry1_arr,ry2_arr,ry1_arr[0]],Color='pink'
    ;cgplot,vU,bspec_U,psym=10,color='blue',/overplot
    cgplot,vU,rspec_U,psym=10,color='red', /overplot
    cgplot,vL,rspec_L,psym=10,color='green',/overplot
    cgplot, !x.CRange,[0.,0.], LineStyle=0, Color='black',/overplot
    cgplot,[vLpeak,vLpeak],[0,!y.crange[1]],linestyle=1,/over
    
    cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(d)',color='black' 
    cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]-0.15*(!y.CRange[1]-!y.CRange[0]),$
      rank[i],color='black'  
    
    ;plot CO outflow contours ?and WISE?
    cropfits,/dataform,dathdrL,[rpeak1[k]+outsz/60.,rpeak1[k]-outsz/60.],$
      [rpeak2[k]-outsz/60.,rpeak2[k]+outsz/60.],[rw0[k],2*vLpeak-rw0[k]],output=lcr
    ldata=smooth(lcr.dat,[1,1,3],/edge_mirror)
    liimap=total(ldata,3)*abs(sxpar(lcr.hdr,'CDELT3'))/1000.
    liimap=congrid(liimap,3*n_elements(liimap[*,0]),3*n_elements(liimap[0,*]),cubic=-0.5,/MINUS_ONE)
    liimap=smooth(liimap,[3,3],/edge_mirror)
    ;lmax=max(liimap[round(n_elements(liimap[*,0])/2.-n_elements(liimap[*,0])/12.-1):round(n_elements(liimap[*,0])/2.+n_elements(liimap[*,0])/12.-1),$
    ;  round(n_elements(liimap[0,*])/2.-n_elements(liimap[0,*])/12.-1):round(n_elements(liimap[0,*])/2.+n_elements(liimap[0,*])/12.-1)])
    cropfits,/dataform,dathdr,[rpeak1[k]+outsz/60.,rpeak1[k]-outsz/60.],$
      [rpeak2[k]-outsz/60.,rpeak2[k]+outsz/60.],[rw0[k],rw1[k]],output=rcr
    rdata=smooth(rcr.dat,[1,1,3],/edge_mirror)
    riimap=total(rdata,3)*abs(sxpar(rcr.hdr,'CDELT3'))/1000.
    riimap=congrid(riimap,3*n_elements(riimap[*,0]),3*n_elements(riimap[0,*]),cubic=-0.5,/MINUS_ONE)
    riimap=smooth(riimap,[3,3],/edge_mirror)
    rlevels=rmax[i]*((rlow[i]<0.8)-0.01+indgen(15)*(1-(rlow[i]<0.8))*0.2)  
    crp1=sxpar(rcr.hdr,'CRPIX1')
    crv1=sxpar(rcr.hdr,'CRVAL1')
    del1=sxpar(rcr.hdr,'CDELT1')
    axs1=sxpar(rcr.hdr,'NAXIS1')
    l_l=(360.+(0.5-crp1)*del1+crv1) mod 360
    l_r=(360.+(axs1+0.5-crp1)*del1+crv1) mod 360
    crp2=sxpar(rcr.hdr,'CRPIX2')
    crv2=sxpar(rcr.hdr,'CRVAL2')
    del2=sxpar(rcr.hdr,'CDELT2')
    axs2=sxpar(rcr.hdr,'NAXIS2')
    b_d=(0.5-crp2)*del2+crv2
    b_u=(axs2+0.5-crp2)*del2+crv2
    x_range=[l_l,l_r]
    y_range=[b_d,b_u]
    ratio=abs((b_u-b_d)/(l_l-l_r))
    ;pos1=[100./xsize1,75./ysize1,800./xsize1,587.5/ysize1]
    pos1=p1
    ;print,pos1
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase
    cgloadct,53,clip=[0,255-64],ncolors=15
    ;cgimage,bytscl(liimap),position=pos1,/noerase
    cgcontour,liimap,/onimage,nlevels=15,label=0,c_colors=indgen(15),/fill
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')
    cgloadct,3,ncolors=15,/reverse,clip=[64,255-64]
    cgcontour,riimap,/onimage,levels=rlevels,label=0,c_colors=indgen(15),thick=3
    ;arrow
    x_0=rpeak1[k]+(outsz/45)*sqrt(0.5)
    y_0=rpeak2[k]-(outsz/45)*sqrt(0.5)
    x_1=rpeak1[k]-(outsz/45)*sqrt(0.5)
    y_1=rpeak2[k]+(outsz/45)*sqrt(0.5)
    cgplot,/over,glon[i],glat[i],psym=cgsymcat(16),color='purple'
    cgarrow,x_0,y_0,x_1,y_1,/data,color='black',/clip,hsize=!D.X_SIZE/128.,hthick=2.0
;    if STRMID(region,0,1) eq 'G' then IDstr= STRMID(region,0,1)+STRMID(region,0,1,/reverse)+'-'+num2str(num[i]) $
;      else IDstr= STRUPCASE(STRMID(region,7,5))+'-'+num2str(num[i])
      IDstr=num2str(index[i])
    cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),$
      IDstr,color='black'
      
    points=(2*!PI/99.0)*findgen(100)
    beam_x=l_l-0.01+24.5/3600.0*cos(points)
    beam_y=b_d+0.01+24.5/3600.0*sin(points)
   ; print,l_l,scalex_c,scaley_c,distance_phy,dist_phy
 
    cgplot,beam_x,beam_y, psym=0,SYMSIZE=1,/data,color='black',xrange=x_range,yrange=y_range,Aspect=ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')

    cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(c)',color='black'
      
    beam_xx=min(beam_x)
    scalex_c=[beam_xx-0.01,beam_xx-0.03]
    scaley_c=[b_d+0.01,b_d+0.01]
    distance_phy=round((0.02*3600.d*distance1[i]*1000.d/206265)*100)/100000.d
    dist_phy=num2str(distance_phy)
    dist_phy=strmid(dist_phy,0,4) 
    cgplot,scalex_c,scaley_c, psym=0,SYMSIZE=1,/data,color='black',xrange=x_range,yrange=y_range,Aspect=ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos1,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})') 
    cgtext, beam_xx-0.01 ,b_d+0.015,dist_phy+' pc',color='black',/data
    
    wiseBname='wise_4.6_'+num2str(bnum[i])+'-'+num2str(rnum[i])+'.fits'
    wiseGname='wise_12_'+num2str(bnum[i])+'-'+num2str(rnum[i])+'.fits'
    wiseRname='wise_22_'+num2str(bnum[i])+'-'+num2str(rnum[i])+'.fits'
    fits_read,'./wisedata/'+wiseBname,datwB,hdrwB
    fits_read,'./wisedata/'+wiseGname,datwG,hdrwG
    fits_read,'./wisedata/'+wiseRname,datwR,hdrwR
    dathdrwB=list(datwB,hdrwB)
    dathdrwG=list(datwG,hdrwG)
    dathdrwR=list(datwR,hdrwR)
    whdr=headfits('./wisedata/'+wiseBname)
    cl_l=l_l+sxpar(whdr,'CDELT1')
    cl_r=l_r-sxpar(whdr,'CDELT1')
    cb_d=b_d+sxpar(whdr,'CDELT2')
    cb_u=b_u-sxpar(whdr,'CDELT2')
    cropfits,dathdrwB,[cl_r,cl_l],[cb_u,cb_d],output=wB,/dataform
    cropfits,dathdrwG,[cl_r,cl_l],[cb_u,cb_d],output=wG,/dataform
    cropfits,dathdrwR,[cl_r,cl_l],[cb_u,cb_d],output=wR,/dataform
    wiseimg_B=wise_scl(wB.dat)
    wiseimg_G=wise_scl(wG.dat)
    wiseimg_R=wise_scl(wR.dat)  
    wiseimg_RGB=[wiseimg_R,wiseimg_G,wiseimg_B]
    pos0=p0
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      AxisColor='black',position=pos0,/noerase,ytitle='Galactic Latitude (!Uo!N)',xtitle=textoidl('Galactic Longitude (^{o})')
    cgimage,wiseimg_RGB,position=pos0,/noerase
    cgplot,[0],[0],xrange=x_range,yrange=y_range,aspect = ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      position=pos0,/noerase,xtickformat='(a1)',ytickformat='(a1)'
    cgloadct,3,ncolors=15,/reverse,clip=[64,255-64]
    cgcontour,riimap,/onimage,levels=rlevels,label=0,c_colors=indgen(15)
    cgplot,/over,glon[i],glat[i],psym=cgsymcat(16),color='purple'
    
    cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(a)',color='white'
    cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]-0.1*(!y.CRange[1]-!y.CRange[0]),$
      textoidl('WISE 4.6 12 22 \mum'),color='white'  
    scalex_w=[l_l-0.01,l_l-0.03]
    scaley_w=[b_d+0.01,b_d+0.01]
    distance_phy=round((0.02*3600.d*distance1[i]*1000.d/206265)*100)/100000.d
    dist_phy=num2str(distance_phy)
    dist_phy=strmid(dist_phy,0,4)
    cgplot,scalex_w,scaley_w, psym=0,SYMSIZE=1,/data,color='white',xrange=x_range,yrange=y_range,Aspect=ratio,xminor=4,yminor=4,$
      xtickinterval=round((max(x_range)-min(x_range))*10*ratio)/25d,ytickinterval=round((max(y_range)-min(y_range))*10)/25d,$
      Axiscolor='white',position=pos0,/noerase,xtickformat='(a1)',ytickformat='(a1)';ytitle='(a1)',xtitle=textoidl('Galactic Longitude (^{o})')
    cgtext, l_l-0.01 ,b_d+0.015,dist_phy+' pc',color='white',/data

    ;pv diagram
    cropfits,dathdr,vrange,dim='v',/dataform,output=pvout
    pv=mkpvbelt(pvout.dat,pvout.hdr,[x_0,x_1],[y_0,y_1],2,/gal)
    pvdata=pv.dat
    pvdata=smooth(pvdata,[3,3],/edge_mirror)
    pvhdr=pv.hdr
    pvp=sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
    pvv=sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
    pvvrg=[min(pvv),max(pvv)]
    pvprg=[min(pvp),max(pvp)]*60.0-mean([min(pvp),max(pvp)])*60;*sqrt(2)   
    cgplot,[0],[0],/nodata,/noerase,position=p2,xrange=pvvrg,yrange=pvprg,$
      ytickformat='(i)',xtickformat='(a1)',yminor=5
    cgloadct,50,clip=[64,240]
;    cgimage,bytscl(pvdata),/overplot;,ncolors=numc;,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline,outcolor='blue'
;    
;    cropfits,dathdrL,vrange,dim='v',/dataform,output=pvoutL
;    pvL=mkpvbelt(pvoutL.dat,pvoutL.hdr,[x_0,x_1],[y_0,y_1],2,/gal)
;    pvdataL=pvL.dat
;    pvdataL=smooth(pvdataL,[3,3],/edge_mirror)
;    pvhdrL=pvL.hdr
;    pvpL=sxpar(pvhdrL,'CRVAL2')+(dindgen(sxpar(pvhdrL,'NAXIS2'))-sxpar(pvhdrL,'CRPIX2')+1)*sxpar(pvhdrL,'CDELT2')
;    pvvL=sxpar(pvhdrL,'CRVAL1')+(dindgen(sxpar(pvhdrL,'NAXIS1'))-sxpar(pvhdrL,'CRPIX1')+1)*sxpar(pvhdrL,'CDELT1')
;    pvvrgL=[min(pvvL),max(pvvL)]
;    pvprgL=[min(pvpL),max(pvpL)]*60.0-mean([min(pvpL),max(pvpL)])*60;*sqrt(2)
;    pvlevels=max(pvdataL[round(n_elements(pvdataL[*,0])/2.-n_elements(pvdataL[*,0])/6.-1):round(n_elements(pvdataL[*,0])/2.+n_elements(pvdataL[*,0])/6.-1),$
;      round(n_elements(pvdataL[0,*])/2.-n_elements(pvdataL[0,*])/6.-1):round(n_elements(pvdataL[0,*])/2.+n_elements(pvdataL[0,*])/6.-1)])$
;      *(0.05+indgen(10)*0.1)
    pvlevels=max(pvdata[round(n_elements(pvdata[*,0])/2.-n_elements(pvdata[*,0])/6.-1):round(n_elements(pvdata[*,0])/2.+n_elements(pvdata[*,0])/6.-1),$
      round(n_elements(pvdata[0,*])/2.-n_elements(pvdata[0,*])/6.-1):round(n_elements(pvdata[0,*])/2.+n_elements(pvdata[0,*])/6.-1)])$
      *(0.05+indgen(10)*0.1)
    cgcontour,pvdata,label=0,/onimage,/fill,level=pvlevels,c_linestyle=0,/outline;',outcolor='green'
    cgplot,[0],[0],/nodata,/noerase,position=p2,xrange=pvvrg,yrange=pvprg,$
      ytickinterval=2.0,ytitle="Position (')",yminor=5,xtitle='LSR Velocity (km s!U-1!N)'
    cgPlot, !x.crange,[0,0], LineStyle=1, /overplot
    cgPlot, [rw0[k],rw0[k]],!y.CRange, LineStyle=1, /overplot
    cgPlot, [rw1[k],rw1[k]],!y.CRange, LineStyle=1, /overplot
    cgPlot,[vLpeak,vLpeak],!y.crange,linestyle=0,/over
    
    pvpx=mean([x_0,x_1])
    pvpxs=strmid(num2str(pvpx),0,7)
    pvpy=mean([y_0,y_1])
    if strmid(num2str(pvpy),0,1) eq '-' then pvpys=strmid(num2str(pvpy),0,6) else pvpys=strmid(num2str(pvpy),0,5)
    pvstring='offset from '+pvpxs+' '+pvpys+' at PA 45'+textoidl('^{o}')
    cgtext,!x.CRange[0]+0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[1]+0.05*(!y.CRange[1]-!y.CRange[0]), pvstring, color='black'
    cgtext,!x.CRange[1]-0.1*(!x.CRange[1]-!x.CRange[0]),!y.CRange[0]+0.05*(!y.CRange[1]-!y.CRange[0]),$
      '(b)',color='black'    
    cgps_close
  endif
  
endfor

  
if ~file_test('outfigures') then spawn,'mkdir outfigures'
;spawn,'rm ./outfigures/*.eps'
spawn,'mv out_*.eps ./outfigures/'

end