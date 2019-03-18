pro cal_phy

;read outflowcat.cat
;read blue/red_out.cat
;read corresponding fits
;draw figures
;calculations

;cd,'/home/lee/W3/'
cd,'/media/alpha/W3/'

;if ~keyword_set(region) then region='region_B'
;if ~keyword_set(Lvrange) then Lvrange=[-85,-55]
;distance=6
;if ~keyword_set(region) then region='region_C_I'
;if ~keyword_set(Lvrange) then Lvrange=[-60,-40]
;distance=2
;if ~keyword_set(region) then region='region_C_II'
;if ~keyword_set(Lvrange) then Lvrange=[-40,-34]
;distance=2
;if ~keyword_set(region) then region='region_C_III'
;if ~keyword_set(Lvrange) then Lvrange=[-34,-28]
;distance=2
;if ~keyword_set(region) then region='region_D'
;if ~keyword_set(Lvrange) then Lvrange=[-29,-20]
;distance=1.3
;if ~keyword_set(region) then region='region_E_I'
;if ~keyword_set(Lvrange) then Lvrange=[-20,-13]
;distance=0.6
;if ~keyword_set(region) then region='region_E_II'
;if ~keyword_set(Lvrange) then Lvrange=[-13,-5]
;distance=0.6
if ~keyword_set(region) then region='region_F'
if ~keyword_set(Lvrange) then Lvrange=[-5,10]
distance=0.6

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
;if ~keyword_set(region) then region='GGMC4'
;if ~keyword_set(Lvrange) then Lvrange=[-3,5]
;distance=2
;if ~keyword_set(region) then region='lynds'
;if ~keyword_set(Lvrange) then Lvrange=[-3,3]
;distance=0.4
;if ~keyword_set(region) then region='region_C_III'
;if ~keyword_set(Lvrange) then Lvrange=[-40,-20]
;distance=2.0
;if ~keyword_set(region) then region='swallow'
;if ~keyword_set(Lvrange) then Lvrange=[12,18]
;distance=3.8
;if ~keyword_set(region) then region='horn'
;if ~keyword_set(Lvrange) then Lvrange=[12,18]
;distance=3.8
;if ~keyword_set(region) then region='remote'
;if ~keyword_set(Lvrange) then Lvrange=[18,28]
;distance=8.5

;cd,'/home/lee/W3/'+region+'/candidates'
cd,'/media/alpha/W3/'+region+'/candidates'
L_lim=0.6
def_rmsU=0.25
def_rmsL=0.3
outsz=(8./distance)>3
outsz=outsz<5

readcol,'blue_out.cat',bsn,bpeak1,bpeak2,bw0,bw1,format='I,F,F,F,F',stringskip='#'
readcol,'red_out.cat',rsn,rpeak1,rpeak2,rw0,rw1,format='I,F,F,F,F',stringskip='#'
readcol,'outflowcat.cat',num,glon,glat,bluev0,bluev1,redv0,redv1,tag,bnum,rnum,index,format='I,F,F,F,F,F,F,A,I,I,I',stringskip='#'
readcol,'distance.cat',id,distance1,format='I,F',stringskip='#'

;openw,lobeb,/get_lun,'lobe_blue.cat'
;openw,lober,/get_lun,'lobe_red.cat'
;printf,lobeb,'# bn','Glon','Glat','bw[0]','bw[1]','spec','pv','con',$
;  format='(a4,a9,a7,2a6,3a5)'
;printf,lober,'# rn','Glon','Glat','rw[0]','rw[1]','spec','pv','con',$
;  format='(a4,a9,a7,2a6,3a5)'
;for i=0, n_elements(bsn)-1 do printf,lobeb,bsn[i],bpeak1[i],bpeak2[i],bw0[i],bw1[i],'-','-','-',$
;  format='(i4,f9.3,f7.3,2f6.1,3a5)'
;for i=0, n_elements(rsn)-1 do printf,lober,rsn[i],rpeak1[i],rpeak2[i],rw0[i],rw1[i],'-','-','-',$
;  format='(i4,f9.3,f7.3,2f6.1,3a5)'
;free_lun,lobeb,lober
;
;openw,para,/get_lun,'lobe_info.cat'
;printf,para,'#osn','tag','lsn','Glon','Glat','w0','w1','spec','pv','con','rank',$
;  format='(3a4,a9,a7,2a8,4a6)'
;;printf,para,'#####################################################################'
;for i=0, n_elements(num)-1 do begin
;  j=where(bsn eq bnum[i])
;  if j ne -1 then printf,para,num[i],'B',bsn[j],bpeak1[j],bpeak2[j],bw0[j],bw1[j],'-','-','-','-',$
;    format='(i4,a4,i4,f9.3,f7.3,2f8.1,4a6)' $
;  else printf,para,num[i],'B','-','-','-','-','-','-','-','-','-',format='(i4,2a4,a9,a7,2a8,4a6)'
;  k=where(rsn eq rnum[i])
;  if k ne -1 then printf,para,'    ','R',rsn[k],rpeak1[k],rpeak2[k],rw0[k],rw1[k],'-','-','-','-',$
;    format='(2a4,i4,f9.3,f7.3,2f8.1,4a6)' $
;  else printf,para,'    ','R','-','-','-','-','-','-','-','-','-',format='(3a4,a9,a7,2a8,4a6)'
;  ;printf,para,'#####################################################################'
;endfor
;free_lun,para

fits_read,'../U_C.fits',datUa,hdrUa
fits_read,'../L_C.fits',datLa,hdrLa
;fits_read,region+'_L2a_mask.fits',datL2a,hdrL2a
fits_read,'../Lpeakv0.fits',peakvmap,vhdr
fits_read,'../Ubvmap.fits',Ubvmap,bvhdr
fits_read,'../Urvmap.fits',Urvmap,rvhdr
dathdr=list(datUa,hdrUa)
dathdrL=list(datLa,hdrLa)

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
  
openw,phy,/get_lun,'derived_phy2.cat'
printf,phy,'#osn','index','tag','Glon','Glat','v_c','del_v','dist','len','mass','moment','E','t','L_flow',$
  format='(a4,2a6,a7,a7,a7,a8,a6,a4,5a11)'
;openw,con,/get_lun,'con_start.cat'
;printf,con,'#num','tag','bsn','rsn','blow','rlow',format='(4a4,2a7)'
;printf,con,'#  blow','bmax','rlow','rmax',format='(4a7)'

;openw,phy,/get_lun,'derived_phy.cat'
;printf,phy,'#osn','tag','lsn','Glon','Glat','area','acm2','low','len','v_c','del_v','mass','moment','E','t','L_flow',$
;  format='(3a4,a9,a7,2a7,2a7,2a7,5a11)'
;openw,con,/get_lun,'con_start.cat'
;;printf,con,'#num','tag','bsn','rsn','blow','rlow',format='(4a4,2a7)'
;printf,con,'#  blow','bmax','rlow','rmax',format='(4a7)'


;measure physical properties
T_ex=30

for i=0, n_elements(num)-1 do begin
 ; print,i
  distance=distance1[i]/1000.0
  ;if (bluev1(i) ge -5.0) or (redv0(i) gt -5.0) then distance=0.2 else distance=0.8
;  if ((bluev1(i) ge -5.0) and (bluev1(i) lt 20)) then distance=0.2 else distance=0.8
;  if ((redv0(i) gt -5.0) and (redv0(i) lt 20)) then distance=0.2 else distance=0.8
  j=where(bsn eq bnum[i])
  k=where(rsn eq rnum[i])
  if j ne -1 then begin
    cropfits,/dataform,dathdr,[bpeak1[j]+outsz/60.,bpeak1[j]-outsz/60.],[bpeak2[j]-outsz/60.,bpeak2[j]+outsz/60.],[bw0[j],bw1[j]],output=bcr
    bcv3=double(sxpar(bcr.hdr,'CRVAL3'))
    bcp3=double(sxpar(bcr.hdr,'CRPIX3'))
    bdl3=double(sxpar(bcr.hdr,'CDELT3'))
    bdata=smooth(bcr.dat,[1,1,3],/edge_mirror)
    
    specUB=reform(bdata[0,0,*])
    vUB=((indgen(n_elements(specUB))+1-bcp3)*bdl3+bcv3)/1000.0
    
    biimap=total(bdata,3)*abs(bdl3)/1000.
    biimap=congrid(biimap,3*n_elements(biimap[*,0]),3*n_elements(biimap[0,*]),cubic=-0.5,/MINUS_ONE)
    biimap=smooth(biimap,[3,3],/edge_mirror)
    bmax=max(biimap[round(n_elements(biimap[*,0])/2.-n_elements(biimap[*,0])/12.-1):round(n_elements(biimap[*,0])/2.+n_elements(biimap[*,0])/12.-1),$
    round(n_elements(biimap[0,*])/2.-n_elements(biimap[0,*])/12.-1):round(n_elements(biimap[0,*])/2.+n_elements(biimap[0,*])/12.-1)])
    ;area,length
    db_area=((0.5/180.0)*(!pi/180)*1000.0*distance)^2
    cal_area,biimap,lmt=0.45,area=barea,mask=bmask,low=blow
    blue_area=barea*db_area
    blue_size=barea*(0.5/3.)^2
    blen=1.77*sqrt(blue_area)
    ;mean integrated intensity
    b_I=total(biimap*bmask)/total(bmask)
    ;delta velocity
    vLpeak=mean(peakvmap[(coord2pix(hdrLa,bpeak1[j],1)-1):(coord2pix(hdrLa,bpeak1[j],1)+1),$
      (coord2pix(hdrLa,bpeak2[j],2)-1):(coord2pix(hdrLa,bpeak2[j],2)+1)],/nan)
    if finite(vLpeak,/nan) eq 1 then begin
      vLpeak=Lvrange[0]
      print,'ALEEEEEEEEEEEEEERT!!!!!!!!!!!!!!!!!!!!!!!!!!'
    endif
    blue_Dv=bw0[j]-vLpeak
    blue_t=abs(blen/blue_Dv)*3.086e13/(365*24*3600.0)
    ;bcv3=double(sxpar(bcr.hdr,'CRVAL3'))
    ;bcp3=double(sxpar(bcr.hdr,'CRPIX3'))
    ;bdl3=double(sxpar(bcr.hdr,'CDELT3'))
    ;bspec_mean=total(bdata[round((size(bdata))[1]/2.-1):round((size(bdata))[1]/2.+1),round((size(bdata))[2]/2.-1):round((size(bdata))[2]/2.+1)],1)
    ;bspec_mean=total(bspec_mean,1)/4.
    ;blue_v=((indgen(n_elements(bspec_mean))+1-bcp3)*bdl3+bcv3)/1000.0
    ;mean_bv=total(bspec_mean*(blue_v-vLpeak))/total(bspec_mean)
    ;mean_bv2=total(bspec_mean*((blue_v-vLpeak)^2))/total(bspec_mean)
    mean_bv=bw0[j]-vLpeak+(bw1[j]-bw0[j])/2.5
    mean_bv2=mean_bv^2
    
    ;physical parameters
    bNH2=4.2e17*b_I*T_ex/exp(-5.5/T_ex)
    
    n1=n_elements(bdata[*,0,0])
    n2=n_elements(bdata[0,*,0])
    n3=n_elements(bdata[0,0,*])
    cropfits,/dataform,dathdrL,[bpeak1[j]+outsz/60.,bpeak1[j]-outsz/60.],[bpeak2[j]-outsz/60.,bpeak2[j]+outsz/60.],[bw0[j],bw1[j]],output=bcrL
    bdataL=smooth(bcrL.dat,[1,1,3],/edge_mirror)
    
    bcv3L=double(sxpar(bcrL.hdr,'CRVAL3'))
    bcp3L=double(sxpar(bcrL.hdr,'CRPIX3'))
    bdl3L=double(sxpar(bcrL.hdr,'CDELT3'))
    
    specLB=reform(bdataL[0,0,*])
    vLB=((indgen(n_elements(specLB))+1-bcp3L)*bdl3L+bcv3L)/1000.0
    
    for l=0,n1-1 do begin
      for m=0,n2-1 do begin
        if bmask[l,m] eq 1 then begin
          for n=0,n3-1 do begin
            bindexv=where((vLB gt vUB[n]-0.1) and (vLB lt vUB[n]+0.1))
            T13B=bdataL[l,m,bindexv]
            T13B=max(T13B)
            if T13B le 0.25 then begin
               bNH2=bNH2+4.2e17*bdata[l,m,n]*T_ex/exp(-5.5/T_ex)*abs(bdl3)/1000. 
            endif else begin
              print,i+1,'B'
               tau12=70*(-alog(1-T13B/26.565))
               bNH2=bNH2+4.2e17*bdata[l,m,n]*T_ex/exp(-5.5/T_ex)*abs(bdl3)/1000.*tau12/(1-exp(-1.*tau12))
            endelse
          endfor
        endif
      endfor
     endfor
     
    bMH2=1.0/6.3e19*bNH2*blue_area
    bPH2=bMH2*abs(mean_bv)
    bEH2=0.5*bMH2*(mean_bv2)*1.9818E30*1.0E6*1e7 ;[J=kg*m^2*s^-2]
    bLm=bEH2/(blue_t*365*24*3600)
    
    printf,phy,'',index[i],'Blue',bpeak1[j],bpeak2[j],vLpeak,abs(mean_bv),distance1[i],blen,bMH2,bPH2,bEH2,blue_t,bLm,$
      format='(a5,i5,a5,f9.3,f7.3,2f7.2,i5,f7.2,5e11.2)'
 ;     print,distance
  endif ;else printf,phy,'',index[i],'Blue','-','-','-','-','-','-','-','-','-','-','-',format='(a5,i5,12a)'
  
;      printf,phy,num[i],'B',bsn[j],bpeak1[j],bpeak2[j],blue_area,blue_size,blow,blen,vLpeak,mean_bv,bMH2,bPH2,bEH2,blue_t,bLm,$
;      format='(i4,a4,i4,f9.3,f7.3,2f7.2,2f7.2,2f7.2,5e11.2)'
;  endif else printf,phy,num[i],'B','-','-','-','-','-','-','-','-','-','-','-','-','-','-',format='(i4,2a4,a9,a7,2a7,4a7,5a11)'
  
  
  if k ne -1 then begin
    cropfits,/dataform,dathdr,[rpeak1[k]+outsz/60.,rpeak1[k]-outsz/60.],[rpeak2[k]-outsz/60.,rpeak2[k]+outsz/60.],[rw0[k],rw1[k]],output=rcr
    rdata=smooth(rcr.dat,[1,1,3],/edge_mirror)
    
    rcv3=double(sxpar(rcr.hdr,'CRVAL3'))
    rcp3=double(sxpar(rcr.hdr,'CRPIX3'))
    rdl3=double(sxpar(rcr.hdr,'CDELT3'))
    
    specUR=reform(rdata[0,0,*])
    vUR=((indgen(n_elements(specUR))+1-rcp3)*rdl3+rcv3)/1000.0
    
    riimap=total(rdata,3)*abs(sxpar(rcr.hdr,'CDELT3'))/1000.
    riimap=congrid(riimap,3*n_elements(riimap[*,0]),3*n_elements(riimap[0,*]),cubic=-0.5,/MINUS_ONE)
    riimap=smooth(riimap,[3,3],/edge_mirror)
    rmax=max(riimap[round(n_elements(riimap[*,0])/2.-n_elements(riimap[*,0])/12.-1):round(n_elements(riimap[*,0])/2.+n_elements(riimap[*,0])/12.-1),$
    round(n_elements(riimap[0,*])/2.-n_elements(riimap[0,*])/12.-1):round(n_elements(riimap[0,*])/2.+n_elements(riimap[0,*])/12.-1)])
    dr_area=((0.5/180.0)*(!pi/180)*1000.0*distance)^2
    cal_area,riimap,lmt=0.45,area=rarea,mask=rmask,low=rlow
    red_area=rarea*dr_area
    red_size=rarea*(0.5/3.)^2
    rlen=1.77*sqrt(red_area)
    ;print,red_area
    R_I=total(Riimap*Rmask)/total(Rmask)
    ;delta velocity
;    if k eq 1 then begin
;      vLpeak=-2.5
;    endif else begin  
    vLpeak=mean(peakvmap[(coord2pix(hdrLa,rpeak1[k],1)-1):(coord2pix(hdrLa,rpeak1[k],1)+1),$
      (coord2pix(hdrLa,rpeak2[k],2)-1):(coord2pix(hdrLa,rpeak2[k],2)+1)],/nan)
;    endelse  
    if finite(vLpeak,/nan) eq 1 then begin
      vLpeak=Lvrange[1]
      print,'ALEEEEEEEEEEEEEERT!!!!!!!!!!!!!!!!!!!!!!!!!!'
    endif
    red_Dv=rw1[k]-vLpeak
    red_t=abs(rlen/red_Dv)*3.086e13/(365*24*3600.0)
    ;mean lobe velocity
    ;rcv3=double(sxpar(rcr.hdr,'CRVAL3'))
    ;rcp3=double(sxpar(rcr.hdr,'CRPIX3'))
    ;rdl3=double(sxpar(rcr.hdr,'CDELT3'))
    ;rspec_mean=total(rdata[round((size(rdata))[1]/2.-1):round((size(rdata))[1]/2.+1),round((size(rdata))[2]/2.-1):round((size(rdata))[2]/2.+1)],1)
    ;rspec_mean=total(rspec_mean,1)/4.
    ;red_v=((indgen(n_elements(rspec_mean))+1-rcp3)*rdl3+rcv3)/1000.0
    ;mean_rv=total(rspec_mean*(red_v-vLpeak))/total(rspec_mean)
    ;mean_rv2=total(rspec_mean*((red_v-vLpeak)^2))/total(rspec_mean)
    mean_rv=rw0[k]-vLpeak+(rw1[k]-rw0[k])/2.5
    mean_rv2=mean_rv^2
    
    ;physical parameters
    rNH2=4.2e17*r_I*T_ex/exp(-5.5/T_ex)
    
    n1=n_elements(rdata[*,0,0])
    n2=n_elements(rdata[0,*,0])
    n3=n_elements(rdata[0,0,*])
    cropfits,/dataform,dathdrL,[rpeak1[j]+outsz/60.,rpeak1[j]-outsz/60.],[rpeak2[j]-outsz/60.,rpeak2[j]+outsz/60.],[rw0[j],rw1[j]],output=rcrL
    rdataL=smooth(rcrL.dat,[1,1,3],/edge_mirror)
    
    rcv3L=double(sxpar(rcrL.hdr,'CRVAL3'))
    rcp3L=double(sxpar(rcrL.hdr,'CRPIX3'))
    rdl3L=double(sxpar(rcrL.hdr,'CDELT3'))
    
    specLR=reform(rdataL[0,0,*])
    vLR=((indgen(n_elements(specLR))+1-rcp3L)*rdl3L+rcv3L)/1000.0
    
    for l=0,n1-1 do begin
      for m=0,n2-1 do begin
        if rmask[l,m] eq 1 then begin
          for n=0,n3-1 do begin
            rindexv=where((vLR gt vUR[n]-0.1) and (vLR lt vUR[n]+0.1))
            T13B=rdataL[l,m,rindexv]
            if T13B le 0.25 then begin
              rNH2=rNH2+4.2e17*rdata[l,m,n]*T_ex/exp(-5.5/T_ex)*abs(rdl3)/1000.
            endif else begin
              print,i+1,'R'
            tau12=70*(-alog(1-T13B/26.565))
            rNH2=rNH2+4.2e17*rdata[l,m,n]*T_ex/exp(-5.5/T_ex)*abs(rdl3)/1000.*tau12/(1-exp(-1.*tau12))
            endelse
          endfor
        endif
      endfor
     endfor
     
    rMH2=1.0/6.3e19*rNH2*red_area
    rPH2=rMH2*abs(mean_rv)
    rEH2=0.5*rMH2*(mean_rv2)*1.9818E30*1.0E6*1e7 ;[J=kg*m^2*s^-2]
    rLm=rEH2/(red_t*365*24*3600)
    
    l=where(bsn eq bnum[i])
    if l eq -1 then begin
    printf,phy,'',index[i],'Red',rpeak1[k],rpeak2[k],vLpeak,abs(mean_rv),distance1[i],rlen,rMH2,rPH2,rEH2,red_t,rLm,$
      format='(a5,i5,a5,f9.3,f7.3,2f7.2,i5,f7.2,5e11.2)'
    endif else printf, phy,'','','Red',rpeak1[k],rpeak2[k],vLpeak,abs(mean_rv),distance1[i],rlen,rMH2,rPH2,rEH2,red_t,rLm,$
      format='(a5,a5,a5,f9.3,f7.3,2f7.2,i5,f7.2,5e11.2)'
 ;     print,distance
    endif ;else printf,phy,'',index[i],'Red','-','-','-','-','-','-','-','-','-','-','-',format='(a5,i5,12a)'
    
;        printf,phy,'','R',rsn[k],rpeak1[k],rpeak2[k],red_area,red_size,rlow,rlen,vLpeak,mean_rv,rMH2,rPH2,rEH2,red_t,rLm,$
;      format='(a4,a4,i4,f9.3,f7.3,2f7.2,2f7.2,2f7.2,5e11.2)'
;    endif else printf,phy,'','R','-','-','-','-','-','-','-','-','-','-','-','-','-','-',format='(a4,2a4,a9,a7,2a7,4a7,5a11)'
    
    if j eq -1 then begin
      blow=0 & bmax=0
    endif
    if k eq -1 then begin
      rlow=0 & rmax=0
    endif
    ;printf,con,num[i],tag[i],bnum[i],rnum[i],blow,rlow,format='(i4,a4,2i4,2f7.2)'
  ;  printf,con,blow,bmax,rlow,rmax,format='(4f7.2)'
endfor

free_lun,phy;,con

;draw outflows


end