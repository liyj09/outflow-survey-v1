function checkrep,listarr
record=list()
poslist=list()
for i=0, n_elements(listarr)-1 do begin
  pos=where(listarr[i] eq listarr)
  if n_elements(pos) ge 2 && where(listarr[i] eq record) eq -1 then begin
    record.add,listarr[i]
    poslist.add,pos
  endif
endfor
return,{rec:record,pos:poslist}
end

pro bipolarcmp,region;,bcat,rcat,Lpeakvmap

;if ~keyword_set(region) then region='BFS52'
;distance=2
;if ~keyword_set(region) then region='GGMC1'
;distance=2
;if ~keyword_set(region) then region='GGMC2'
;distance=2
;if ~keyword_set(region) then region='GGMC3'
;distance=2
if ~keyword_set(region) then region='region_C_III'
distance=2.0
;if ~keyword_set(region) then region='lynds'
;distance=0.4
;if ~keyword_set(region) then region='west'
;distance=0.6
;if ~keyword_set(region) then region='swallow'
;distance=3.8
;if ~keyword_set(region) then region='horn'
;distance=3.8
;if ~keyword_set(region) then region='remote'
;distance=8.5

cd,'/home/lee/W3/'+region+'/candidates/'

;Lpeakvmap='Lpeakv0.fits'
bcat='blue_out.cat'
rcat='red_out.cat'

;fits_read,Lpeakvmap,dat,hdr
readcol,bcat,bn,bgl,bgb,bw0,bw1,format='I,F,F,F,F',stringskip='#'
readcol,rcat,rn,rgl,rgb,rw0,rw1,format='I,F,F,F,F',stringskip='#'

;print,bn[1],bgl[1],bgb[1],bw0[1],bw1[1]

;stop

openw,ocat,/get_lun,'outflowcat.txt'
printf,ocat,'#num','Glon','Glat','bw[0]','bw[1]','rw[0]','rw[1]','tag','b_n','r_n',$
  format='(a4,a9,a7,4a6,a3,2a4)'
outcount=0
bipolar=list()
norepeat=list()
rest=list()
used=list()
r_used=list()
climit=(0.1/distance)<0.1
climit=climit>0.06
for i=0, n_elements(bn)-1 do begin
  dis=sqrt((bgl[i]-rgl)^2+(bgb[i]-rgb)^2)
  if (where(dis lt climit))[0] ne -1 then begin
    pos=where(dis lt climit)
    bipolar.add,pos,/extract
  endif    
  ;mindis=min(sqrt((bgl[i]-rgl)^2+(bgb[i]-rgb)^2),minpos)
  ;if mindis lt climit then bipolar.add,minpos
endfor
  ;rep=checkrep(bipolar)
for i=0, n_elements(bn)-1 do begin
  if (where(i eq used))[0] eq -1 then begin
    mindis=min(sqrt((bgl[i]-rgl)^2+(bgb[i]-rgb)^2),minpos)
    if mindis lt climit && (where(minpos eq norepeat))[0] eq -1 then begin
      if n_elements(where(minpos eq bipolar)) ge 2 then begin
        outcount++
        mp=min(sqrt((bgl-rgl[minpos])^2+(bgb-rgb[minpos])^2),posb)
        newgl=(bgl[posb]+rgl[minpos])/2.
        newgb=(bgb[posb]+rgb[minpos])/2.
        printf,ocat,outcount,newgl,newgb,bw0[posb],bw1[posb],rw0[minpos],rw1[minpos],'D',bn[posb],rn[minpos],$
          format='(i4,f9.3,f7.3,4f6.1,a3,2i4)'
        norepeat.add,minpos
        r_used.add,minpos
        if posb ne i then rest.add,i
        used.add,posb
        ;print,norepeat
      endif else begin
        r_used.add,minpos
        outcount++
        newgl=(bgl[i]+rgl[minpos])/2.
        newgb=(bgb[i]+rgb[minpos])/2.
        printf,ocat,outcount,newgl,newgb,bw0[i],bw1[i],rw0[minpos],rw1[minpos],'D',bn[i],rn[minpos],$
          format='(i4,f9.3,f7.3,4f6.1,a3,2i4)'
      endelse
    endif else begin
      outcount++
      printf,ocat,outcount,bgl[i],bgb[i],bw0[i],bw1[i],99999,99999,'B',bn[i],0,$
         format='(i4,f9.3,f7.3,2f6.1,2i6,a3,2i4)'
    endelse
  endif
endfor

for i=0, n_elements(rest)-1 do begin
  n=rest[i]
  outcount++
  printf,ocat,outcount,bgl[n],bgb[n],bw0[n],bw1[n],99999,99999,'B',bn[n],0,$
    format='(i4,f9.3,f7.3,2f6.1,2i6,a3,2i4)'
endfor

for i=0, n_elements(rn)-1 do begin
  if (where(r_used eq i))[0] eq -1 then begin
    outcount++
    printf,ocat,outcount,rgl[i],rgb[i],99999,99999,rw0[i],rw1[i],'R',0,rn[i],$
      format='(i4,f9.3,f7.3,2i6,2f6.1,A3,2i4)'
  endif
endfor

free_lun,ocat

readcol,'outflowcat.txt',num,glon,glat,wb0,wb1,wr0,wr1,tag,bnum,rnum,format='I,F,F,F,F,F,F,A,I,I',stringskip='#'
openw,cat,/get_lun,'outflowcat.cat'
printf,cat,'#num','Glon','Glat','bw[0]','bw[1]','rw[0]','rw[1]','tag','bnum','rnum',$
  format='(a4,a9,a7,4a8,a4,2a5)'

seq=sort(glon)
glon0=glon[seq]
glat0=glat[seq]
rep=checkrep(glon0)
for i=0, n_elements(rep.rec)-1 do begin
  pos=(rep.pos)[i]
  barr=glat0[pos]
  seqb=sort(barr)+pos[0]
  seq[pos]=seq[seqb]  
endfor
glon=glon[seq]
glat=glat[seq]
wb0=wb0[seq]
wb1=wb1[seq]
wr0=wr0[seq]
wr1=wr1[seq]
tag=tag[seq]
bnum=bnum[seq]
rnum=rnum[seq]

for i=0,n_elements(num)-1 do begin
  printf,cat,i+1,glon[i],glat[i],wb0[i],wb1[i],wr0[i],wr1[i],tag[i],bnum[i],rnum[i],$
    format='(i4,f9.3,f7.3,4f8.1,a4,2i5)'
endfor

free_lun,cat



end